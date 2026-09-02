import os
import sys


def _candidate_paths():
    base = getattr(sys, "_MEIPASS", None)
    if not base:
        return []
    return [
        os.path.join(base, "PyQt6", "Qt6", "bin"),
        os.path.join(base, "PyQt6", "Qt6", "plugins"),
        base,
    ]


for path in _candidate_paths():
    if os.path.isdir(path):
        try:
            os.add_dll_directory(path)
        except (AttributeError, OSError):
            pass
        os.environ["PATH"] = path + os.pathsep + os.environ.get("PATH", "")

