#!/usr/bin/env bash
# Author: Nathaniel Sun
# Date: 2026-06-17
# Description:
#   macOS project setup: create .venv and install the 1D-fluid model dependencies.
#
#   NumPy is installed from a standard wheel because custom OpenBLAS source builds provide no
#   measurable benefit to the fluid/optics hot path: np.fft is pocketfft (never BLAS), SciPy's
#   sparse SuperLU links its own BLAS independently of NumPy, and the fluid/optics
#   kernels are Numba/elementwise. Thread and CPU-topology configuration is
#   handled at runtime by the src/backend/ module (auto-detects performance vs.
#   efficiency cores and the best Numba threading layer).

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INSTALL_ROOT="${INSTALL_ROOT:-$SCRIPT_DIR}"
VENV_DIR="$INSTALL_ROOT/.venv"
REQUIREMENTS="$INSTALL_ROOT/requirements.txt"

echo "=== 0. Checking Xcode Command Line Tools (git) ==="
if ! xcode-select -p >/dev/null 2>&1; then
    echo "Xcode Command Line Tools are required (git)."
    echo "A dialog should open; complete the install, then re-run this script."
    xcode-select --install
    exit 1
fi

echo "=== 1. Checking Python version ==="
python3 - <<'PY'
import sys

minimum = (3, 11)
if sys.version_info < minimum:
    detected = ".".join(map(str, sys.version_info[:3]))
    required = ".".join(map(str, minimum))
    raise SystemExit(f"Python {required}+ is required; found {detected}.")
PY

echo "=== 2. Creating Python virtual environment (.venv) ==="
if [[ ! -d "$VENV_DIR" ]]; then
    python3 -m venv "$VENV_DIR"
fi
PYTHON="$VENV_DIR/bin/python"
export PATH="$VENV_DIR/bin:$PATH"

"$PYTHON" -m pip install --upgrade pip

echo "=== 3. Installing 1D-fluid model dependencies (NumPy, SciPy, Numba, ...) ==="
if [[ ! -f "$REQUIREMENTS" ]]; then
    echo "ERROR: ${REQUIREMENTS} not found." >&2
    exit 1
fi
"$PYTHON" -m pip install -r "$REQUIREMENTS"
"$PYTHON" -m pip install --no-deps -e "$INSTALL_ROOT"

echo "=== 4. Verifying install ==="
"$PYTHON" -c 'import numpy, scipy, numba; print(f"NumPy {numpy.__version__}, SciPy {scipy.__version__}, Numba {numba.__version__}")'

echo "===================================================="
echo "Setup complete. Activate with:"
printf '  source %q\n' "$VENV_DIR/bin/activate"
echo "Thread count is auto-selected at runtime (src/backend/);"
echo "the compute backend (src/backend/) auto-selects the Numba threading layer."
echo
echo "===================================================="
