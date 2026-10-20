#!/usr/bin/env bash
# Author: Nathaniel Sun
# Date: 2026-06-17
# Description:
#   Linux project setup: create .venv and install the 1D-fluid model dependencies.
#
#   NumPy is installed from a standard wheel because custom OpenBLAS source builds provide no
#   measurable benefit to the fluid/optics hot path (np.fft is pocketfft, SciPy's SuperLU links its
#   own BLAS, fluid/optics kernels are Numba/elementwise). Thread/CPU configuration
#   is handled at runtime by the src/backend/ module, which also auto-selects the
#   best available Numba threading layer (tbb on Linux/x86 if installed, else
#   workqueue) -- see the optional ``tbb`` entry in requirements.txt.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INSTALL_ROOT="${INSTALL_ROOT:-$SCRIPT_DIR}"
VENV_DIR="$INSTALL_ROOT/.venv"
REQUIREMENTS="$INSTALL_ROOT/requirements.txt"

echo "=== 1. Ensuring Python venv support is present (per distro) ==="
if [[ -f /etc/os-release ]]; then
    # shellcheck disable=SC1091
    source /etc/os-release
fi

case "${ID:-}" in
    debian|ubuntu)
        sudo apt-get update
        sudo apt-get install -y python3-venv python3-dev git
        ;;
    fedora)
        sudo dnf install -y python3 python3-devel git
        ;;
    rhel|centos|rocky|almalinux)
        PKG_MANAGER="$(command -v dnf || command -v yum || true)"
        [[ -n "$PKG_MANAGER" ]] && sudo "$PKG_MANAGER" install -y python3 python3-devel git
        ;;
    arch)
        sudo pacman -Syu --noconfirm --needed python git
        ;;
    *)
        echo "Unknown distro '${ID:-?}'. Proceeding; ensure python3-venv and git are installed."
        ;;
esac

echo "=== 2. Checking Python version ==="
python3 - <<'PY'
import sys

minimum = (3, 11)
if sys.version_info < minimum:
    detected = ".".join(map(str, sys.version_info[:3]))
    required = ".".join(map(str, minimum))
    raise SystemExit(f"Python {required}+ is required; found {detected}.")
PY

echo "=== 3. Creating Python virtual environment (.venv) ==="
if [[ ! -d "$VENV_DIR" ]]; then
    python3 -m venv "$VENV_DIR"
fi
PYTHON="$VENV_DIR/bin/python"
export PATH="$VENV_DIR/bin:$PATH"

"$PYTHON" -m pip install --upgrade pip

echo "=== 4. Installing 1D-fluid model dependencies (NumPy, SciPy, Numba, ...) ==="
if [[ ! -f "$REQUIREMENTS" ]]; then
    echo "ERROR: ${REQUIREMENTS} not found." >&2
    exit 1
fi
"$PYTHON" -m pip install -r "$REQUIREMENTS"
"$PYTHON" -m pip install --no-deps -e "$INSTALL_ROOT"

echo "=== 5. Verifying install ==="
"$PYTHON" -c 'import numpy, scipy, numba; print(f"NumPy {numpy.__version__}, SciPy {scipy.__version__}, Numba {numba.__version__}")'

echo "===================================================="
echo "Setup complete. Activate with:"
printf '  source %q\n' "$VENV_DIR/bin/activate"
echo "Thread count is auto-selected at runtime (src/backend/);"
echo "the compute backend (src/backend/) auto-selects the Numba threading layer."
echo
echo "===================================================="
