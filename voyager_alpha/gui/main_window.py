from __future__ import annotations

import sys
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QIcon, QPixmap
from PyQt6.QtWidgets import (
    QButtonGroup,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from .asteroid_workspace import AsteroidWorkspace
from .exoplanet_workspace import ExoplanetWorkspace
from .theme import APP_STYLE


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Voyager Alpha")
        self.setWindowIcon(QIcon(resource_path("assets/voyager-alpha.ico")))
        self.setMinimumSize(1280, 760)
        self.resize(1720, 980)
        self.setStyleSheet(APP_STYLE)
        self._build_ui()

    def _build_ui(self):
        root = QWidget()
        root.setObjectName("appRoot")
        self.setCentralWidget(root)
        layout = QVBoxLayout(root)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._build_header())

        self.stack = QStackedWidget()
        self.asteroid_workspace = AsteroidWorkspace(self)
        self.exoplanet_workspace = ExoplanetWorkspace(self)
        self.stack.addWidget(self.asteroid_workspace)
        self.stack.addWidget(self.exoplanet_workspace)
        layout.addWidget(self.stack, 1)

    def _build_header(self):
        frame = QFrame()
        frame.setObjectName("titleBar")
        frame.setFixedHeight(58)
        layout = QHBoxLayout(frame)
        layout.setContentsMargins(12, 5, 12, 5)
        layout.setSpacing(8)

        logo_plate = QFrame()
        logo_plate.setObjectName("logoPlate")
        logo_plate.setFixedSize(46, 46)
        logo_layout = QVBoxLayout(logo_plate)
        logo_layout.setContentsMargins(3, 3, 3, 3)
        logo = QLabel()
        logo.setAlignment(Qt.AlignmentFlag.AlignCenter)
        logo.setPixmap(
            QPixmap(resource_path("assets/voyager-alpha-logo.png")).scaled(
                40,
                40,
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
        )
        logo.setToolTip("Voyager Alpha")
        logo_layout.addWidget(logo)
        layout.addWidget(logo_plate)

        brand = QVBoxLayout()
        brand.setSpacing(0)
        self.lbl_brand_name = QLabel("Voyager Alpha")
        self.lbl_brand_name.setObjectName("brandName")
        subtitle = QLabel("Asteroid discovery | Exoplanet transit photometry")
        subtitle.setObjectName("brandSub")
        brand.addWidget(self.lbl_brand_name)
        brand.addWidget(subtitle)
        layout.addLayout(brand)
        layout.addSpacing(26)

        self.module_group = QButtonGroup(self)
        self.module_group.setExclusive(True)
        self.btn_asteroid = self._module_button("Asteroid Hunter", 0)
        self.btn_exoplanet = self._module_button("Exoplanet Inspection", 1)
        self.btn_asteroid.setChecked(True)
        layout.addWidget(self.btn_asteroid)
        layout.addWidget(self.btn_exoplanet)
        layout.addStretch(1)

        astrohub = QLabel("ASTROHUB.COM.TR")
        astrohub.setObjectName("caption")
        astrohub.setToolTip("www.astrohub.com.tr")
        layout.addWidget(astrohub)
        help_button = QPushButton("?")
        help_button.setObjectName("windowButton")
        help_button.setToolTip("Programın iki bilimsel modülü ve ölçüm sınırları hakkında bilgi")
        help_button.clicked.connect(self.show_help)
        layout.addWidget(help_button)
        settings_button = QPushButton("⚙")
        settings_button.setObjectName("windowButton")
        settings_button.setToolTip("ASTAP ve çalışma ortamı bilgileri")
        settings_button.clicked.connect(self.show_settings)
        layout.addWidget(settings_button)
        return frame

    def _module_button(self, text: str, index: int):
        button = QPushButton(text)
        button.setObjectName("moduleTab")
        button.setCheckable(True)
        button.clicked.connect(lambda _checked=False, value=index: self.stack.setCurrentIndex(value))
        self.module_group.addButton(button, index)
        return button

    def show_help(self):
        QMessageBox.information(
            self,
            "Voyager Alpha",
            "Asteroid Hunter, zaman sıralı FITS karelerinde kalibrasyon, WCS çözümü, alt piksel hizalama, "
            "statik gökyüzü çıkarımı ve zamana bağlı tracklet bağlantısı uygular. Tracklet, aynı hareketli cismin "
            "birden fazla karedeki ölçümlerinin zaman sıralı zinciridir.\n\n"
            "Exoplanet Inspection, seçilen hedef ve karşılaştırma yıldızlarından diferansiyel ışık eğrisi üretir. "
            "Her iki modülün sonuçları aday niteliğindedir; gözlem koşulları ve bağımsız doğrulama gerekir.",
        )

    def show_settings(self):
        astap = Path(r"C:\Program Files\astap\astap.exe")
        state = "Bulundu" if astap.exists() else "Bulunamadı"
        QMessageBox.information(
            self,
            "Çalışma ortamı",
            f"ASTAP: {state}\n{astap}\n\n"
            "Asteroid Guided Flow, WCS eksikse ASTAP çözümünü otomatik dener. "
            "MPC gözlemevi kodu yalnızca ADES/MPC raporu oluştururken gereklidir.",
        )

    def closeEvent(self, event):
        self.asteroid_workspace.shutdown_workers()
        self.exoplanet_workspace.shutdown_workers()
        event.accept()


def resource_path(relative: str) -> str:
    if hasattr(sys, "_MEIPASS"):
        root = Path(sys._MEIPASS) / "voyager_alpha"
    else:
        root = Path(__file__).resolve().parents[1]
    return str(root / relative)
