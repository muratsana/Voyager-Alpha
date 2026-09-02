import sys
import ctypes
from pathlib import Path

from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QApplication

if __package__:
    from .gui.main_window import MainWindow
else:
    try:
        from voyager_alpha.gui.main_window import MainWindow
    except ImportError:
        from gui.main_window import MainWindow

APP_USER_MODEL_ID = "tr.com.astrohub.voyageralpha"


def configure_windows_app_identity() -> None:
    """Keep the taskbar group and icon attached to Voyager Alpha."""

    if sys.platform != "win32":
        return
    try:
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(APP_USER_MODEL_ID)
    except (AttributeError, OSError):
        pass


def application_icon_path() -> Path:
    if hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / "voyager_alpha" / "assets" / "voyager-alpha.ico"
    return Path(__file__).resolve().parent / "assets" / "voyager-alpha.ico"


if __name__ == '__main__':
    configure_windows_app_identity()
    app = QApplication(sys.argv)
    app.setApplicationName("Voyager Alpha")
    app.setApplicationDisplayName("Voyager Alpha")
    app.setOrganizationName("Astrohub")
    app.setDesktopFileName(APP_USER_MODEL_ID)
    app.setWindowIcon(QIcon(str(application_icon_path())))
    window = MainWindow()
    window.showMaximized()
    sys.exit(app.exec())
