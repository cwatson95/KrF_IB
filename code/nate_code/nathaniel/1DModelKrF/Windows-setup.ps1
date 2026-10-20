#!/usr/bin/env pwsh
# Author: Nathaniel Sun
# Date: 2026-06-17
# Description:
#   Windows project setup: create .venv and install the 1D-fluid model dependencies.
#
#   NumPy is installed from a standard wheel because custom MSYS2/mingw OpenBLAS source builds
#   provide no measurable benefit to the fluid/optics hot path (np.fft is pocketfft, SciPy's
#   SuperLU links its own BLAS, and fluid/optics kernels are Numba/elementwise). Thread/CPU
#   configuration is handled at
#   runtime by the src/backend/ module.
# NOTE: You may need: Set-ExecutionPolicy RemoteSigned -Scope CurrentUser

$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir

$InstallRoot = if ($env:INSTALL_ROOT) { $env:INSTALL_ROOT } else { $ScriptDir }
$VenvDir = Join-Path $InstallRoot ".venv"
$Requirements = Join-Path $InstallRoot "requirements.txt"

function Test-Command {
    param([string]$Name)
    return [bool](Get-Command $Name -ErrorAction SilentlyContinue)
}

try {
    Write-Host "=== 1. Checking prerequisites ==="
    if (-not (Test-Command python)) {
        throw "Python was not found on PATH. Install Python from https://www.python.org/downloads/ and re-run."
    }

    $pyVersion = python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
    $major, $minor = $pyVersion.Split(".")
    if ([int]$major -lt 3 -or ([int]$major -eq 3 -and [int]$minor -lt 11)) {
        throw "Python 3.11+ is required; found $pyVersion."
    }

    Write-Host "=== 2. Creating Python virtual environment (.venv) ==="
    if (-not (Test-Path $VenvDir)) {
        python -m venv $VenvDir
    }
    $Python = Join-Path $VenvDir "Scripts\python.exe"
    $env:Path = (Join-Path $VenvDir "Scripts") + ";" + $env:Path
    & $Python -m pip install --upgrade pip

    Write-Host "=== 3. Installing 1D-fluid model dependencies (NumPy, SciPy, Numba, ...) ==="
    if (-not (Test-Path $Requirements)) {
        throw "requirements.txt not found at $Requirements"
    }
    & $Python -m pip install -r $Requirements
    & $Python -m pip install --no-deps -e $InstallRoot

    Write-Host "=== 4. Verifying install ==="
    & $Python -c "import numpy, scipy, numba; print(f'NumPy {numpy.__version__}, SciPy {scipy.__version__}, Numba {numba.__version__}')"

    Write-Host "===================================================="
    Write-Host "Setup complete. Activate with: .\$VenvDir\Scripts\Activate.ps1"
    Write-Host "Thread count is auto-selected at runtime (src/backend/);"
    Write-Host "the compute backend (src/backend/) auto-selects the Numba threading layer."
    Write-Host ""
    Write-Host "===================================================="
}
catch {
    Write-Error $_
    exit 1
}
