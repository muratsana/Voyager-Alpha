$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$python = Join-Path $root ".venv\Scripts\python.exe"
$pyinstaller = Join-Path $root ".venv\Scripts\pyinstaller.exe"
$dist = Join-Path $root "dist"
$build = Join-Path $root "build"
$spec = Join-Path $root "packaging\VoyagerAlpha-onefile.spec"

if (-not (Test-Path -LiteralPath $python)) {
    throw "Virtual environment not found. Create it first: python -m venv .venv"
}

if (-not (Test-Path -LiteralPath $pyinstaller)) {
    & $python -m pip install pyinstaller
}

& $python (Join-Path $root "voyager_alpha\diagnostics.py")
if ($LASTEXITCODE -ne 0) { throw "Diagnostics failed with exit code $LASTEXITCODE" }
& $python -m unittest discover -s (Join-Path $root "voyager_alpha\tests")
if ($LASTEXITCODE -ne 0) { throw "Tests failed with exit code $LASTEXITCODE" }
& $python -m pip check
if ($LASTEXITCODE -ne 0) { throw "pip check failed with exit code $LASTEXITCODE" }

Push-Location (Join-Path $root "packaging")
try {
    & $pyinstaller --noconfirm --clean --distpath $dist --workpath $build $spec
    if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}

$exe = Join-Path $dist "Voyager-Alpha.exe"
if (-not (Test-Path -LiteralPath $exe)) {
    throw "Build finished but executable was not found: $exe"
}

$hash = Get-FileHash -Algorithm SHA256 -LiteralPath $exe
Write-Host "Built: $exe"
Write-Host "SHA256: $($hash.Hash)"
