import importlib.util
import shutil
import sys
from pathlib import Path


REQUIRED_MODULES = ["numpy", "scipy", "exoplanet_core", "astropy", "astroquery", "sep", "PyQt6", "pyqtgraph"]


def module_available(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def run_diagnostics() -> int:
    print("Voyager Alpha diagnostics")
    print(f"Python: {sys.version.split()[0]}")
    missing = []
    for name in REQUIRED_MODULES:
        ok = module_available(name)
        print(f"{name}: {'OK' if ok else 'MISSING'}")
        if not ok:
            missing.append(name)

    astap = (
        shutil.which("astap.exe")
        or shutil.which("astap")
        or (r"C:\Program Files\astap\astap.exe" if Path(r"C:\Program Files\astap\astap.exe").is_file() else None)
    )
    print(f"ASTAP: {astap if astap else 'not found on PATH'}")

    if missing:
        print("")
        print("Install missing packages with:")
        print("python -m pip install -r requirements.txt")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run_diagnostics())
