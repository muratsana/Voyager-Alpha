from __future__ import annotations

import csv
import os
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from astropy.wcs import WCS
from PyQt6.QtCore import QTimer, Qt
from PyQt6.QtGui import QColor, QPixmap
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QSplitter,
    QStyle,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..core.exoplanet_catalog import ExoplanetCatalog, predict_nearest_transit
from ..core.exoplanet_worker import (
    CatalogUpdateWorker,
    ExoplanetWorker,
    PhotometryPreflightWorker,
    SequencePlateSolveWorker,
)
from ..core.fits_io import read_fits_image
from ..core.metadata import inspect_sequence
from ..core.wcs_cache import WcsSolutionCache, pixel_to_sky
from .viewer import FitsViewer
from .widgets import AnalysisLog, Filmstrip, KeyValueTable, PipelineProgressStrip, camera_metadata_rows, section
from .workers import PreviewCacheWorker


class ExoplanetWorkspace(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.fits_files: list[str] = []
        self.frames = []
        self.preview_cache = {}
        self.current_index = 0
        self.current_folder = ""
        self.target_xy = None
        self.target_radec = None
        self.target_selection_frame = 0
        self.comparison_xy = []
        self.comparison_radec = []
        self.comparison_selection_frames = []
        self.selection_mode = None
        self.master_paths = {"bias": "", "dark": "", "flat": ""}
        self.preview_worker = None
        self.worker = None
        self.plate_worker = None
        self.preflight_worker = None
        self.catalog_worker = None
        self.wcs_cache = WcsSolutionCache()
        self.wcs_headers = {}
        self.catalog = ExoplanetCatalog()
        self.catalog_matches = []
        self._catalog_auto_update_checked = False
        self._catalog_update_manual = False
        self._start_after_plate_solve = False
        self._start_after_preflight = False
        self._apply_preflight_references = False
        self.preflight_result = None
        self.result = None
        self._build_ui()
        self._install_tooltips()
        self.update_actions()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(self._build_toolbar())

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(1)
        splitter.addWidget(self._build_left())
        splitter.addWidget(self._build_center())
        splitter.addWidget(self._build_right())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 0)
        splitter.setSizes([285, 1040, 330])
        root.addWidget(splitter, 1)
        root.addWidget(self._build_status())

    def _build_toolbar(self):
        frame = QFrame()
        frame.setObjectName("toolbar")
        layout = QHBoxLayout(frame)
        layout.setContentsMargins(9, 7, 9, 7)
        self.btn_open = QPushButton("Kaynak Klasör")
        self.btn_open.setObjectName("blueButton")
        self.btn_open.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon))
        self.btn_open.clicked.connect(self.load_folder)
        self.btn_open.setMinimumWidth(116)
        layout.addWidget(self.btn_open)
        self.btn_plate_solve = QPushButton("Plate Solve")
        self.btn_plate_solve.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_DialogApplyButton))
        self.btn_plate_solve.clicked.connect(self.start_plate_solve)
        layout.addWidget(self.btn_plate_solve)
        self.btn_catalog = QPushButton("Transit Katalogları")
        self.btn_catalog.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_BrowserReload))
        self.btn_catalog.clicked.connect(self.open_catalog_tab)
        layout.addWidget(self.btn_catalog)
        self.btn_start = QPushButton("Fotometriyi Başlat")
        self.btn_start.setObjectName("primaryButton")
        self.btn_start.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        self.btn_start.clicked.connect(self.start_analysis)
        self.btn_start.setMinimumWidth(154)
        layout.addWidget(self.btn_start)
        self.btn_stop = QPushButton("Durdur")
        self.btn_stop.setObjectName("dangerButton")
        self.btn_stop.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaStop))
        self.btn_stop.clicked.connect(self.stop_analysis)
        self.btn_stop.setMinimumWidth(78)
        layout.addWidget(self.btn_stop)
        layout.addStretch(1)
        self.lbl_session = QLabel("Aktif transit sekansı yok")
        self.lbl_session.setObjectName("muted")
        layout.addWidget(self.lbl_session)
        return frame

    def _build_left(self):
        panel = QFrame()
        panel.setObjectName("leftPanel")
        panel.setMinimumWidth(260)
        panel.setMaximumWidth(320)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(5)

        validation, body = section("Sequence Validation")
        self.validation = KeyValueTable()
        self.validation.set_auto_height(True)
        body.addWidget(self.validation)
        layout.addWidget(validation)

        calibration, body = section("Calibration Masters")
        self.master_labels = {}
        for kind, label in (("bias", "Master Bias"), ("dark", "Master Dark"), ("flat", "Master Flat")):
            row = QHBoxLayout()
            button = QPushButton(label)
            button.setToolTip(f"İsteğe bağlı {label.lower()} FITS kalibrasyon karesini seçer.")
            button.clicked.connect(lambda _checked=False, key=kind: self.choose_master(key))
            value = QLabel("Seçilmedi")
            value.setObjectName("muted")
            row.addWidget(button)
            row.addWidget(value, 1)
            body.addLayout(row)
            self.master_labels[kind] = value
        layout.addWidget(calibration)

        selection, body = section("Photometry Stars")
        star_buttons = QHBoxLayout()
        star_buttons.setSpacing(4)
        self.btn_target = QPushButton("Hedef Seç")
        self.btn_target.setObjectName("blueButton")
        self.btn_target.clicked.connect(lambda: self.set_selection_mode("target"))
        self.btn_auto_references = QPushButton("Otomatik Referans")
        self.btn_auto_references.setObjectName("primaryButton")
        self.btn_auto_references.clicked.connect(lambda: self.start_photometry_preflight(apply_references=True))
        star_buttons.addWidget(self.btn_target)
        star_buttons.addWidget(self.btn_auto_references)
        body.addLayout(star_buttons)
        manual_buttons = QHBoxLayout()
        manual_buttons.setSpacing(4)
        self.btn_comparison = QPushButton("Manuel Referans")
        self.btn_comparison.clicked.connect(lambda: self.set_selection_mode("comparison"))
        self.btn_clear_stars = QPushButton("Temizle")
        self.btn_clear_stars.clicked.connect(self.clear_stars)
        manual_buttons.addWidget(self.btn_comparison)
        manual_buttons.addWidget(self.btn_clear_stars)
        body.addLayout(manual_buttons)
        self.lbl_target = QLabel("Hedef: seçilmedi")
        self.lbl_comparisons = QLabel("Karşılaştırma: 0")
        body.addWidget(self.lbl_target)
        body.addWidget(self.lbl_comparisons)
        self.lbl_preflight = QLabel("Uygunluk: hedef seçilmedi")
        self.lbl_preflight.setObjectName("muted")
        self.lbl_preflight.setWordWrap(True)
        body.addWidget(self.lbl_preflight)
        layout.addWidget(selection)

        settings, body = section("Photometry Settings")
        aperture_row = QHBoxLayout()
        aperture_row.addWidget(QLabel("Aperture radius"))
        self.combo_aperture = QComboBox()
        self.combo_aperture.addItems(["4 px", "5 px", "6 px", "7 px", "8 px", "10 px"])
        self.combo_aperture.setCurrentText("6 px")
        self.combo_aperture.currentTextChanged.connect(self.invalidate_preflight)
        aperture_row.addWidget(self.combo_aperture, 1)
        body.addLayout(aperture_row)
        detrend_row = QHBoxLayout()
        detrend_row.addWidget(QLabel("Detrending"))
        self.combo_detrend = QComboBox()
        self.combo_detrend.addItems(["None", "Linear", "Quadratic"])
        self.combo_detrend.setCurrentText("Quadratic")
        detrend_row.addWidget(self.combo_detrend, 1)
        body.addLayout(detrend_row)
        linearity_row = QHBoxLayout()
        linearity_row.addWidget(QLabel("Linearity limit"))
        self.spin_linearity = QDoubleSpinBox()
        self.spin_linearity.setDecimals(0)
        self.spin_linearity.setRange(0.0, 1000000000.0)
        self.spin_linearity.setSingleStep(1000.0)
        self.spin_linearity.setSpecialValueText("Auto / FITS")
        self.spin_linearity.valueChanged.connect(self.invalidate_preflight)
        linearity_row.addWidget(self.spin_linearity, 1)
        body.addLayout(linearity_row)
        self.chk_recenter = QCheckBox("Recenter stars per frame")
        self.chk_recenter.setChecked(True)
        body.addWidget(self.chk_recenter)
        layout.addWidget(settings)

        guidance, body = section("Measurement Guardrails")
        guidance_text = (
            "Hedef ve karşılaştırma yıldızları doygun olmamalı, benzer parlaklıkta ve tüm sekans boyunca kadrajda "
            "kalmalıdır. Sonuç bir aday göstergesidir; bilimsel doğrulama değildir."
        )
        text = QLabel("Lineerlik · PSF · crowding · ensemble denetimi")
        text.setObjectName("muted")
        text.setToolTip(guidance_text)
        body.addWidget(text)
        layout.addWidget(guidance)
        layout.addStretch(1)
        return panel

    def _build_center(self):
        panel = QFrame()
        panel.setObjectName("centerPanel")
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)

        viewer_box, viewer_layout = section("FITS Photometry Viewer")
        self.viewer = FitsViewer()
        self.viewer.image_clicked.connect(self.image_clicked)
        viewer_layout.addWidget(self.viewer, 1)
        layout.addWidget(viewer_box, 3)

        self.filmstrip = Filmstrip(page_size=6)
        self.filmstrip.frame_selected.connect(self.change_frame)
        layout.addWidget(self.filmstrip)

        curve_box, curve_layout = section("Differential Light Curve")
        curve_toolbar = QHBoxLayout()
        curve_toolbar.addWidget(QLabel("Display"))
        self.combo_curve = QComboBox()
        self.combo_curve.addItems(["Detrended", "Raw relative"])
        self.combo_curve.currentTextChanged.connect(self.plot_result)
        curve_toolbar.addWidget(self.combo_curve)
        curve_toolbar.addStretch(1)
        self.lbl_curve_model = QLabel("Model waiting")
        self.lbl_curve_model.setObjectName("muted")
        curve_toolbar.addWidget(self.lbl_curve_model)
        curve_layout.addLayout(curve_toolbar)
        self.curve_plot = pg.PlotWidget(background="#090d10")
        self.curve_plot.setMinimumHeight(190)
        self.curve_plot.showGrid(x=True, y=True, alpha=0.15)
        self.curve_plot.setLabel("left", "Relative flux")
        self.curve_plot.setLabel("bottom", "Elapsed time", units="min")
        self.curve_plot.getAxis("left").setTextPen(QColor("#9aa9b0"))
        self.curve_plot.getAxis("bottom").setTextPen(QColor("#9aa9b0"))
        curve_layout.addWidget(self.curve_plot)
        layout.addWidget(curve_box, 2)
        return panel

    def _build_right(self):
        panel = QFrame()
        panel.setObjectName("rightPanel")
        panel.setMinimumWidth(300)
        panel.setMaximumWidth(390)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(9, 9, 9, 9)
        layout.setSpacing(6)
        self.right_tabs = QTabWidget()
        self.right_tabs.tabBar().setExpanding(True)
        self.right_tabs.tabBar().setUsesScrollButtons(False)
        self.right_tabs.setStyleSheet("QTabBar::tab { min-width: 48px; padding: 2px 4px; }")

        assessment_tab = QWidget()
        assessment_layout = QVBoxLayout(assessment_tab)
        assessment_layout.setContentsMargins(5, 6, 5, 6)
        assessment_layout.setSpacing(6)
        quality_box, body = section("Photometry Suitability")
        self.lbl_quality_verdict = QLabel("ÖN KONTROL BEKLİYOR")
        self.lbl_quality_verdict.setObjectName("badgeNeutral")
        body.addWidget(self.lbl_quality_verdict)
        self.preflight_table = KeyValueTable()
        self.preflight_table.setFixedHeight(145)
        body.addWidget(self.preflight_table)
        assessment_layout.addWidget(quality_box)
        result_box, body = section("Transit Assessment")
        self.lbl_verdict = QLabel("ANALİZ BEKLİYOR")
        self.lbl_verdict.setObjectName("badgeNeutral")
        body.addWidget(self.lbl_verdict)
        self.result_table = KeyValueTable()
        self.result_table.setMinimumHeight(250)
        body.addWidget(self.result_table)
        self.lbl_quality = QLabel("Quality flags: waiting")
        self.lbl_quality.setObjectName("muted")
        self.lbl_quality.setWordWrap(True)
        body.addWidget(self.lbl_quality)
        assessment_layout.addWidget(result_box, 1)
        self.right_tabs.addTab(assessment_tab, "Result")

        catalog_tab = QWidget()
        catalog_layout = QVBoxLayout(catalog_tab)
        catalog_layout.setContentsMargins(5, 6, 5, 6)
        catalog_layout.setSpacing(6)
        source_box, body = section("Catalog Sources")
        self.catalog_sources_table = QTableWidget(0, 3)
        self.catalog_sources_table.setHorizontalHeaderLabels(["Source", "Records", "State"])
        self.catalog_sources_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.catalog_sources_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.catalog_sources_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.catalog_sources_table.verticalHeader().setVisible(False)
        self.catalog_sources_table.setFixedHeight(150)
        body.addWidget(self.catalog_sources_table)
        self.btn_update_catalog = QPushButton("Tüm Katalogları Güncelle")
        self.btn_update_catalog.setObjectName("primaryButton")
        self.btn_update_catalog.clicked.connect(lambda: self.start_catalog_update(manual=True))
        body.addWidget(self.btn_update_catalog)
        catalog_layout.addWidget(source_box)

        match_box, body = section("Selected Target Match")
        self.lbl_catalog_target = QLabel("WCS çözülüp hedef seçildiğinde katalog eşleşmesi burada görünür.")
        self.lbl_catalog_target.setObjectName("muted")
        self.lbl_catalog_target.setWordWrap(True)
        body.addWidget(self.lbl_catalog_target)
        self.catalog_match_table = QTableWidget(0, 4)
        self.catalog_match_table.setHorizontalHeaderLabels(["Object", "Status", "Sep", "Transit"])
        self.catalog_match_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for column in (1, 2, 3):
            self.catalog_match_table.horizontalHeader().setSectionResizeMode(column, QHeaderView.ResizeMode.ResizeToContents)
        self.catalog_match_table.verticalHeader().setVisible(False)
        body.addWidget(self.catalog_match_table)
        note = QLabel(
            "WCS ve gözlemevi konumu varsa zaman BJD_TDB'ye çevrilir; "
            "eksikse JD_UTC açık uyarıyla korunur."
        )
        note.setObjectName("muted")
        note.setWordWrap(True)
        body.addWidget(note)
        catalog_layout.addWidget(match_box, 1)
        self.right_tabs.addTab(catalog_tab, "Catalog")

        measurement_tab = QWidget()
        measurement_layout = QVBoxLayout(measurement_tab)
        measurement_layout.setContentsMargins(5, 6, 5, 6)
        measurement_box, body = section("Frame Measurements")
        self.measurements = QTableWidget(0, 3)
        self.measurements.setHorizontalHeaderLabels(["Frame", "Target", "Comparisons"])
        self.measurements.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.measurements.verticalHeader().setVisible(False)
        body.addWidget(self.measurements)
        measurement_layout.addWidget(measurement_box, 1)
        self.right_tabs.addTab(measurement_tab, "Frames")

        log_tab = QWidget()
        log_layout = QVBoxLayout(log_tab)
        log_layout.setContentsMargins(5, 6, 5, 6)
        log_box, body = section("Analysis Log")
        self.analysis_log = AnalysisLog()
        body.addWidget(self.analysis_log)
        log_layout.addWidget(log_box, 1)
        self.right_tabs.addTab(log_tab, "Log")

        layout.addWidget(self.right_tabs, 1)

        self.btn_export = QPushButton("Işık Eğrisini CSV Aktar")
        self.btn_export.clicked.connect(self.export_csv)
        layout.addWidget(self.btn_export)
        self.refresh_catalog_status()
        return panel

    def _build_status(self):
        frame = QFrame()
        frame.setObjectName("statusBar")
        layout = QHBoxLayout(frame)
        layout.setContentsMargins(8, 4, 8, 4)
        layout.setSpacing(8)
        self.lbl_status = QLabel("Hazır")
        self.lbl_status.setObjectName("muted")
        self.lbl_status.setMinimumWidth(140)
        self.lbl_status.setMaximumWidth(220)
        self.pipeline_progress = PipelineProgressStrip(
            [
                ("source", "Kaynak"),
                ("astrometry", "Astrometri"),
                ("calibration", "Kalibrasyon"),
                ("align", "Hizalama"),
                ("photometry", "Fotometri"),
                ("detrend", "Detrend"),
                ("model", "Model"),
                ("review", "İnceleme"),
            ]
        )
        self.pipeline_progress.setMaximumWidth(720)
        self.progress = self.pipeline_progress
        self.lbl_file_status = QLabel("")
        self.lbl_file_status.setObjectName("muted")
        self.lbl_file_status.setMinimumWidth(300)
        self.lbl_file_status.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.lbl_progress_percent = QLabel("0%")
        self.lbl_progress_percent.setObjectName("statusPercent")
        self.lbl_progress_percent.setFixedWidth(42)
        self.lbl_progress_percent.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.lbl_status)
        layout.addWidget(self.pipeline_progress)
        layout.addWidget(self.lbl_file_status, 1)
        layout.addWidget(self.lbl_progress_percent)
        return frame

    def _install_tooltips(self):
        self.btn_open.setToolTip("Zaman sıralı transit FITS sekansının bulunduğu klasörü açar.")
        self.btn_plate_solve.setToolTip(
            "WCS'siz referansı ASTAP ile çözer; güvenilir yıldız kaydı bulunan diğer karelere çözümü taşır, "
            "zayıf kayıtta ASTAP'a geri döner. Orijinal FITS değiştirilmez."
        )
        self.btn_catalog.setToolTip("Doğrulanmış transitler, TOI, KOI ve K2 katalog durumunu açar.")
        self.btn_update_catalog.setToolTip("Dört resmi NASA/IPAC veri kümesini transaction güvenliğiyle yerel kataloğa günceller.")
        self.btn_start.setToolTip(
            "Eksikse önce WCS preflight çalıştırır; ardından kalibrasyon, hizalama, aperture fotometri ve "
            "diferansiyel ışık eğrisini üretir."
        )
        self.btn_stop.setToolTip("Çalışan fotometri işini kontrollü biçimde durdurur.")
        self.btn_target.setToolTip("Görüntü üzerinde transit beklenen hedef yıldızı seçmek için etkinleştirir.")
        self.btn_auto_references.setToolTip(
            "Tüm sekansı tarar; doygunluk, kamera lineerliği, SNR, PSF, crowding, parlaklık benzerliği, "
            "sekans içi kararlılık ve mümkünse Gaia/VSX kayıtlarıyla 3-8 referans seçer."
        )
        self.btn_comparison.setToolTip("Manuel referans ekler; analizden önce bilimsel ön kontrolden geçirilir.")
        self.btn_clear_stars.setToolTip("Hedef ve karşılaştırma yıldızı seçimlerini temizler.")
        self.btn_export.setToolTip("Zaman, bağıl akı, hedef akısı ve karşılaştırma akısını CSV dosyasına aktarır.")
        self.viewer.setToolTip("Tek tık seçili yıldızı işaretler; çift tık görüntüyü alana sığdırır.")
        self.combo_aperture.setToolTip("Yıldız akısının toplandığı dairesel aperture yarıçapıdır; yıldız FWHM değerinin yaklaşık 1.5-2 katı seçilir.")
        self.spin_linearity.setToolTip(
            "Kameranızın lineerlik testinde ölçülen en yüksek güvenilir ADU değeridir. 0 bırakılırsa FITS başlığı "
            "aranır; yalnız BITPIX aralığı bulunursa sonuç CAUTION olarak işaretlenir."
        )
        self.combo_detrend.setToolTip("Yavaş atmosfer ve sistem eğilimlerini sabit, doğrusal veya ikinci derece taban çizgisiyle ayırır.")
        self.chk_recenter.setToolTip("Global hizalamadan sonra her yıldızın yerel merkezini kare bazında yeniden ölçer.")
        self.combo_curve.setToolTip("Ham diferansiyel akı ile detrend edilmiş analiz eğrisi arasında geçiş yapar.")

    def showEvent(self, event):
        super().showEvent(event)
        if self._catalog_auto_update_checked:
            return
        self._catalog_auto_update_checked = True
        if os.environ.get("QT_QPA_PLATFORM", "").lower() != "offscreen":
            QTimer.singleShot(600, self.maybe_auto_update_catalog)

    def load_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Transit FITS klasörünü seç", self.current_folder)
        if folder:
            self.scan_folder(folder)

    def scan_folder(self, folder: str):
        self.shutdown_workers()
        self.current_folder = folder
        self.fits_files = find_fits_files(folder)
        self.frames = []
        self.result = None
        self.preview_cache.clear()
        self.wcs_headers.clear()
        self.curve_plot.clear()
        self.measurements.setRowCount(0)
        self.pipeline_progress.reset()
        self.lbl_progress_percent.setText("0%")
        self.lbl_file_status.clear()
        self.clear_stars()
        if not self.fits_files:
            self.log("ERROR|Source|Klasörde FITS dosyası bulunamadı")
            self.update_actions()
            return
        try:
            self.frames = inspect_sequence(self.fits_files)
            self.frames.sort(key=lambda item: (item.midpoint_jd is None, item.midpoint_jd or 0.0, item.file_path))
            self.fits_files = [item.file_path for item in self.frames]
            for file_path in self.fits_files:
                header, _method = self.wcs_cache.resolve(file_path)
                if header is not None:
                    self.wcs_headers[file_path] = header
            self.refresh_validation()
            self.filmstrip.set_frames(self.fits_files)
            self.lbl_session.setText(f"{len(self.fits_files)} frame · differential photometry")
            self.preview_worker = PreviewCacheWorker(self.fits_files, parent=self)
            self.preview_worker.preview_ready.connect(self.preview_ready)
            self.preview_worker.failed.connect(lambda message: self.log(f"ERROR|Preview|{message}"))
            self.preview_worker.start()
            self.log("INFO|Source|Transit FITS sequence loaded")
        except Exception as exc:
            self.log(f"ERROR|Metadata|{exc}")
        self.update_actions()

    def refresh_validation(self):
        timed = sum(item.midpoint_jd is not None for item in self.frames)
        shape = self.frames[0].shape if self.frames else None
        exposure = self.frames[0].exposure_seconds if self.frames else None
        wcs_count = sum(file_path in self.wcs_headers for file_path in self.fits_files)
        self.validation.set_rows([
            ("FITS Frames", str(len(self.frames)), "ok" if len(self.frames) >= 5 else "warn"),
            ("Time / WCS", f"{timed}/{len(self.frames)} · {wcs_count}/{len(self.frames)}", "ok" if timed == len(self.frames) and wcs_count == len(self.frames) else "warn"),
            ("Frame / Exposure", f"{shape[1]} x {shape[0]} · {exposure:.3f} s" if shape and exposure else "-", None),
        ] + (camera_metadata_rows(self.frames[0].camera) if self.frames else []))

    def preview_ready(self, index, display, original_shape, thumbnail):
        self.preview_cache[int(index)] = (display, tuple(original_shape))
        self.filmstrip.set_thumbnail(int(index), QPixmap.fromImage(thumbnail))
        if int(index) == self.current_index:
            self.refresh_frame()

    def change_frame(self, index: int):
        if not self.fits_files:
            return
        self.current_index = max(0, min(int(index), len(self.fits_files) - 1))
        self.filmstrip.select(self.current_index)
        self.refresh_frame()

    def refresh_frame(self):
        cached = self.preview_cache.get(self.current_index)
        if cached:
            self.viewer.set_display_image(cached[0], cached[1])
            known = []
            for index, value in enumerate(self.comparison_xy):
                radec = self.comparison_radec[index] if index < len(self.comparison_radec) else None
                display_xy = self._radec_to_current_pixel(radec) if radec else value
                if display_xy is not None:
                    known.append({"x": display_xy[0], "y": display_xy[1]})
            target_display = self._radec_to_current_pixel(self.target_radec) if self.target_radec else self.target_xy
            candidates = [] if target_display is None else [{"x": target_display[0], "y": target_display[1]}]
            self.viewer.set_overlays(known=known, candidates=candidates)
        elif self.fits_files:
            try:
                data = read_fits_image(self.fits_files[self.current_index])
                self.viewer.set_image(data, stretch_mode="Auto STF")
            except Exception as exc:
                self.log(f"ERROR|Viewer|{Path(self.fits_files[self.current_index]).name}: {exc}")

    def set_selection_mode(self, mode: str):
        if not self.fits_files:
            return
        self.selection_mode = mode
        self.lbl_status.setText("Görüntüde hedef yıldızı seçin" if mode == "target" else "Görüntüde karşılaştırma yıldızı seçin")

    def image_clicked(self, x: float, y: float):
        radec = self._current_pixel_to_radec(x, y)
        reference_xy = (self._radec_to_reference_pixel(radec) if radec else None) or (x, y)
        if self.selection_mode == "target":
            self.target_xy = reference_xy
            self.target_radec = radec
            self.target_selection_frame = self.current_index
            if radec:
                self.lbl_target.setText(f"Hedef: {radec[0]:.5f}°, {radec[1]:+.5f}°")
            else:
                self.lbl_target.setText(f"Hedef: X {x:.1f}, Y {y:.1f} · WCS bekliyor")
            self.selection_mode = None
            self.match_selected_target()
            self.invalidate_preflight()
            QTimer.singleShot(0, lambda: self.start_photometry_preflight(apply_references=True))
        elif self.selection_mode == "comparison":
            if all(np.hypot(reference_xy[0] - px, reference_xy[1] - py) > 8.0 for px, py in self.comparison_xy):
                self.comparison_xy.append(reference_xy)
                self.comparison_radec.append(radec)
                self.comparison_selection_frames.append(self.current_index)
            self.lbl_comparisons.setText(f"Karşılaştırma: {len(self.comparison_xy)}")
            self.invalidate_preflight()
        self.refresh_frame()
        self.update_actions()

    def clear_stars(self):
        self.target_xy = None
        self.target_radec = None
        self.target_selection_frame = 0
        self.comparison_xy = []
        self.comparison_radec = []
        self.comparison_selection_frames = []
        self.selection_mode = None
        self.preflight_result = None
        if hasattr(self, "lbl_target"):
            self.lbl_target.setText("Hedef: seçilmedi")
            self.lbl_comparisons.setText("Karşılaştırma: 0")
            self.lbl_preflight.setText("Uygunluk: hedef seçilmedi")
            self.lbl_quality_verdict.setText("ÖN KONTROL BEKLİYOR")
            self.preflight_table.set_rows([])
            self.viewer.clear_overlays()
        if hasattr(self, "catalog_match_table"):
            self.catalog_match_table.setRowCount(0)
            self.lbl_catalog_target.setText("WCS çözülüp hedef seçildiğinde katalog eşleşmesi burada görünür.")
        self.update_actions()

    def _current_pixel_to_radec(self, x: float, y: float):
        if not self.fits_files:
            return None
        header = self.wcs_headers.get(self.fits_files[self.current_index])
        if header is None:
            return None
        try:
            return pixel_to_sky(header, x, y)
        except Exception:
            return None

    def _radec_to_reference_pixel(self, radec):
        if not radec or not self.fits_files:
            return None
        header = self.wcs_headers.get(self.fits_files[0])
        if header is None:
            return None
        try:
            x, y = WCS(header).celestial.all_world2pix(radec[0], radec[1], 0)
            return float(x), float(y)
        except Exception:
            return None

    def _radec_to_current_pixel(self, radec):
        if not radec or not self.fits_files:
            return None
        header = self.wcs_headers.get(self.fits_files[self.current_index])
        if header is None:
            return None
        try:
            x, y = WCS(header).celestial.all_world2pix(radec[0], radec[1], 0)
            return float(x), float(y)
        except Exception:
            return None

    def choose_master(self, kind: str):
        path, _ = QFileDialog.getOpenFileName(self, f"Master {kind.title()} seç", self.current_folder, "FITS (*.fit *.fits *.fts *.fits.gz)")
        if path:
            self.master_paths[kind] = path
            self.master_labels[kind].setText(Path(path).name)
            self.master_labels[kind].setToolTip(path)
            self.invalidate_preflight()

    def invalidate_preflight(self, *_args):
        self.preflight_result = None
        if hasattr(self, "lbl_preflight") and self.target_xy is not None:
            self.lbl_preflight.setText("Uygunluk: yeniden kontrol gerekli")
        if hasattr(self, "lbl_quality_verdict"):
            self.lbl_quality_verdict.setText("ÖN KONTROL GÜNCEL DEĞİL")
        self.update_actions()

    def start_photometry_preflight(self, _checked=False, *, apply_references=False, start_after=False):
        if self.target_xy is None or len(self.fits_files) < 5:
            if apply_references:
                QMessageBox.information(self, "Fotometri Ön Kontrolü", "Önce transit hedef yıldızını seçin.")
            return
        if self.preflight_worker and self.preflight_worker.isRunning():
            self._apply_preflight_references = self._apply_preflight_references or bool(apply_references)
            self._start_after_preflight = self._start_after_preflight or bool(start_after)
            return
        self._apply_preflight_references = bool(apply_references)
        self._start_after_preflight = bool(start_after)
        reference_header = self.wcs_headers.get(self.fits_files[0]) if self.fits_files else None
        self.preflight_worker = PhotometryPreflightWorker(
            self.fits_files,
            self.target_xy,
            aperture_radius=float(self.combo_aperture.currentText().split()[0]),
            linearity_limit_adu=self.spin_linearity.value(),
            master_paths=self.master_paths,
            reference_header=reference_header,
            catalog_lookup=reference_header is not None,
            parent=self,
        )
        self.preflight_worker.progress.connect(self.preflight_progress)
        self.preflight_worker.log.connect(self.log)
        self.preflight_worker.result_ready.connect(self.preflight_ready)
        self.preflight_worker.failed.connect(self.preflight_failed)
        self.preflight_worker.finished_scan.connect(self.preflight_finished)
        self.preflight_worker.start()
        self.lbl_preflight.setText("Uygunluk: sekans taranıyor")
        self.lbl_quality_verdict.setText("ÖN KONTROL ÇALIŞIYOR")
        self.pipeline_progress.reset()
        self.pipeline_progress.set_active_stage("photometry")
        self.update_actions()

    def preflight_progress(self, value, message):
        self.pipeline_progress.setValue(int(value))
        self.lbl_progress_percent.setText(f"{int(value)}%")
        self.lbl_file_status.setText(message)
        self.lbl_file_status.setToolTip(message)

    def preflight_ready(self, result):
        self.preflight_result = result
        target = result.target
        if self._apply_preflight_references:
            self.comparison_xy = list(result.recommended_xy)
            self.comparison_selection_frames = [0] * len(self.comparison_xy)
            header = self.wcs_headers.get(self.fits_files[0]) if self.fits_files else None
            self.comparison_radec = [pixel_to_sky(header, *xy) if header is not None else None for xy in self.comparison_xy]
            self.lbl_comparisons.setText(f"Karşılaştırma: {len(self.comparison_xy)} otomatik")
            self.refresh_frame()
        verdict_name = {
            "PASS": "HEDEF UYGUN",
            "CAUTION": "HEDEF KOŞULLU",
            "FAIL": "HEDEF UYGUN DEĞİL",
        }.get(target.status, "HEDEF BELİRSİZ")
        self.lbl_quality_verdict.setText(verdict_name)
        self.lbl_quality_verdict.setObjectName(
            "badgeKnown" if target.status == "PASS" else "badgeCandidate" if target.status == "CAUTION" else "badgeRejected"
        )
        self.lbl_quality_verdict.style().unpolish(self.lbl_quality_verdict)
        self.lbl_quality_verdict.style().polish(self.lbl_quality_verdict)
        peak = f"{target.peak_adu:.0f} / {result.linearity.limit_adu:.0f} ADU"
        self.preflight_table.set_rows([
            ("Target", target.status, "ok" if target.status == "PASS" else "warn"),
            ("Peak / limit", peak, "warn" if (target.peak_fraction or 0) >= 0.9 else None),
            ("PSF / SNR", f"{target.median_fwhm_px:.2f} px · {target.median_snr:.0f}", None),
            ("Valid frames", f"{target.valid_fraction * 100:.0f}%", "ok" if target.valid_fraction >= 0.95 else "warn"),
            ("References", f"{len(result.recommended_xy)} selected", "ok" if len(result.recommended_xy) >= 3 else "warn"),
            ("Linearity", result.linearity.source, "ok" if result.linearity.verified else "warn"),
        ])
        detail = "; ".join([*target.messages, *result.warnings]) or "Tüm kalite kontrolleri geçti"
        self.lbl_preflight.setText(f"Uygunluk: {target.status} · {detail}")
        self.lbl_preflight.setToolTip(detail)
        self.pipeline_progress.complete()
        self.lbl_progress_percent.setText("100%")
        self.log(
            f"{'INFO' if result.analysis_allowed else 'WARN'}|Photometry QC|"
            f"Target {target.status}; peak {peak}; FWHM {target.median_fwhm_px:.2f}px; "
            f"{len(result.recommended_xy)} references"
        )

    def preflight_failed(self, message):
        self.preflight_result = None
        self.lbl_preflight.setText(f"Uygunluk: kontrol başarısız · {message}")
        self.lbl_quality_verdict.setText("ÖN KONTROL BAŞARISIZ")
        if self.pipeline_progress.active_key:
            self.pipeline_progress.set_stage_state(self.pipeline_progress.active_key, "warning")

    def preflight_finished(self):
        should_start = self._start_after_preflight and self.preflight_result is not None
        self._start_after_preflight = False
        self._apply_preflight_references = False
        self.update_actions()
        if should_start:
            if not self.preflight_result.analysis_allowed:
                QMessageBox.warning(
                    self,
                    "Fotometri uygun değil",
                    "Hedef veya referans ensemble bilimsel kalite kapısını geçemedi. Uygunluk ayrıntılarını kontrol edin.",
                )
            else:
                QTimer.singleShot(0, self.start_analysis)

    def start_plate_solve(self, _checked=False, *, start_after=False):
        if not self.fits_files:
            QMessageBox.information(self, "Plate Solve", "Önce transit FITS sekansını yükleyin.")
            return
        if self.plate_worker and self.plate_worker.isRunning():
            self._start_after_plate_solve = self._start_after_plate_solve or start_after
            return
        self._start_after_plate_solve = bool(start_after)
        self.plate_worker = SequencePlateSolveWorker(self.fits_files, cache=self.wcs_cache, parent=self)
        self.plate_worker.progress.connect(self.plate_solve_progress)
        self.plate_worker.log.connect(self.log)
        self.plate_worker.result_ready.connect(self.plate_solve_ready)
        self.plate_worker.finished_scan.connect(self.plate_solve_finished)
        self.plate_worker.start()
        self.pipeline_progress.reset()
        self.pipeline_progress.set_active_stage("astrometry")
        self.lbl_status.setText("Astrometri preflight çalışıyor")
        self.log("INFO|Astrometry|Sequence WCS preflight started")
        self.update_actions()

    def plate_solve_progress(self, value: int, message: str):
        self.pipeline_progress.setValue(int(value))
        self.lbl_progress_percent.setText(f"{int(value)}%")
        self.lbl_file_status.setText(message)
        self.lbl_file_status.setToolTip(message)

    def plate_solve_ready(self, result):
        self.wcs_headers = dict(result.headers)
        self.refresh_validation()
        self._normalize_star_selections_to_reference()
        self.invalidate_preflight()
        self.match_selected_target()
        if result.failures:
            self.lbl_status.setText(f"WCS eksik: {len(result.failures)} kare")
            self.log(f"ERROR|Astrometry|{len(result.failures)} frame could not be solved")
        else:
            propagated = sum(method.startswith("propagated") for method in result.methods.values())
            direct = sum(method == "astap" for method in result.methods.values())
            self.lbl_status.setText("Tüm karelerin WCS çözümü hazır")
            self.log(f"INFO|Astrometry|WCS complete: {direct} ASTAP, {propagated} propagated")
            self.pipeline_progress.complete()
            self.lbl_progress_percent.setText("100%")

    def plate_solve_finished(self):
        should_start = self._start_after_plate_solve and len(self.wcs_headers) == len(self.fits_files)
        self._start_after_plate_solve = False
        self.update_actions()
        if should_start:
            QTimer.singleShot(0, self.start_analysis)
        elif self.fits_files and len(self.wcs_headers) != len(self.fits_files):
            QMessageBox.warning(
                self,
                "Astrometri tamamlanmadı",
                "Bazı kareler çözülemedi. Fotometri başlamadı; Analysis Log içindeki ASTAP hatalarını kontrol edin.",
            )

    def _normalize_star_selections_to_reference(self):
        if self.target_xy is not None and self.target_radec is None and self.fits_files:
            index = min(self.target_selection_frame, len(self.fits_files) - 1)
            header = self.wcs_headers.get(self.fits_files[index])
            if header is not None:
                self.target_radec = pixel_to_sky(header, *self.target_xy)
                self.target_xy = self._radec_to_reference_pixel(self.target_radec) or self.target_xy
        normalized_xy = []
        normalized_radec = []
        for index, xy in enumerate(self.comparison_xy):
            radec = self.comparison_radec[index] if index < len(self.comparison_radec) else None
            if radec is None and self.fits_files:
                frame_index = self.comparison_selection_frames[index] if index < len(self.comparison_selection_frames) else 0
                frame_index = min(frame_index, len(self.fits_files) - 1)
                header = self.wcs_headers.get(self.fits_files[frame_index])
                if header is not None:
                    radec = pixel_to_sky(header, *xy)
            normalized_radec.append(radec)
            normalized_xy.append(self._radec_to_reference_pixel(radec) or xy)
        self.comparison_radec = normalized_radec
        self.comparison_xy = normalized_xy
        self.refresh_frame()

    def open_catalog_tab(self):
        self.right_tabs.setCurrentIndex(1)
        self.refresh_catalog_status()

    def maybe_auto_update_catalog(self):
        if self.catalog.is_stale(max_age_days=7):
            self.start_catalog_update(manual=False)

    def start_catalog_update(self, _checked=False, *, manual=True):
        if self.catalog_worker and self.catalog_worker.isRunning():
            return
        self._catalog_update_manual = bool(manual)
        self.catalog_worker = CatalogUpdateWorker(self.catalog, self)
        self.catalog_worker.progress.connect(self.catalog_update_progress)
        self.catalog_worker.log.connect(self.log)
        self.catalog_worker.result_ready.connect(self.catalog_update_ready)
        self.catalog_worker.failed.connect(self.catalog_update_failed)
        self.catalog_worker.finished_scan.connect(self.catalog_update_finished)
        self.pipeline_progress.reset()
        self.pipeline_progress.set_active_stage("source")
        self.catalog_worker.start()
        self.right_tabs.setCurrentIndex(1)
        self.lbl_catalog_target.setText("Resmi transit katalogları güncelleniyor...")
        self.update_actions()

    def catalog_update_progress(self, value: int, source_key: str, message: str):
        self.lbl_status.setText("Transit katalogları güncelleniyor")
        self.lbl_file_status.setText(message)
        self.lbl_file_status.setToolTip(f"{source_key}: {message}")
        self.lbl_progress_percent.setText(f"{int(value)}%")
        self.pipeline_progress.setValue(int(value))

    def catalog_update_ready(self, counts):
        self.refresh_catalog_status()
        self.match_selected_target()
        self.lbl_status.setText(f"Transit kataloğu hazır: {sum(counts.values()):,} kaynak kaydı")
        self.lbl_progress_percent.setText("100%")
        self.pipeline_progress.complete()
        if self._catalog_update_manual:
            QMessageBox.information(
                self,
                "Transit Katalogları",
                "Güncelleme tamamlandı.\n\n"
                + "\n".join(f"{key}: {value:,}" for key, value in counts.items()),
            )

    def catalog_update_failed(self, message: str):
        self.refresh_catalog_status()
        self.lbl_status.setText("Katalog güncellemesi tamamlanamadı")
        if self.pipeline_progress.active_key:
            self.pipeline_progress.set_stage_state(self.pipeline_progress.active_key, "warning")
        if self._catalog_update_manual:
            QMessageBox.warning(
                self,
                "Transit Katalogları",
                f"Güncelleme tamamlanamadı. Daha önceki başarılı veriler korundu.\n\n{message}",
            )

    def catalog_update_finished(self):
        self.update_actions()

    def refresh_catalog_status(self):
        if not hasattr(self, "catalog_sources_table"):
            return
        rows = self.catalog.source_status()
        self.catalog_sources_table.setRowCount(len(rows))
        for row_index, row in enumerate(rows):
            updated = str(row["updated_utc"])[:10] if row["updated_utc"] else "Never"
            state = f"{row['status']} · {updated}"
            values = (str(row["display_name"]), f"{int(row['record_count']):,}", state)
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                if row["error"]:
                    item.setToolTip(str(row["error"]))
                self.catalog_sources_table.setItem(row_index, column, item)

    def match_selected_target(self):
        if not hasattr(self, "catalog_match_table") or self.target_xy is None:
            return
        if self.target_radec is None:
            reference = self.wcs_headers.get(self.fits_files[0]) if self.fits_files else None
            if reference is None:
                self.lbl_catalog_target.setText("Hedef seçildi; katalog eşleşmesi için plate solve gerekiyor.")
                return
            self.target_radec = pixel_to_sky(reference, *self.target_xy)
        ra, dec = self.target_radec
        self.lbl_catalog_target.setText(f"RA {ra:.6f}° · Dec {dec:+.6f}° · 5 arcmin search")
        self.catalog_matches = self.catalog.cone_search(ra, dec, radius_arcmin=5.0, limit=25)
        self.catalog_match_table.setRowCount(len(self.catalog_matches))
        start_jd = min((frame.midpoint_jd for frame in self.frames if frame.midpoint_jd is not None), default=None)
        end_jd = max((frame.midpoint_jd for frame in self.frames if frame.midpoint_jd is not None), default=None)
        labels = {
            "confirmed": "VERIFIED",
            "candidate": "CANDIDATE",
            "unverified": "UNVERIFIED",
            "false_positive": "FALSE POSITIVE",
        }
        for row, match in enumerate(self.catalog_matches):
            prediction = predict_nearest_transit(match, start_jd, end_jd) if start_jd is not None and end_jd is not None else None
            if prediction is None:
                transit_text = "No ephemeris"
            elif prediction.overlaps_observation:
                transit_text = "IN WINDOW"
            else:
                transit_text = f"Δ {prediction.offset_minutes:+.0f} min"
            values = (
                f"{match.name} · {match.source_key}",
                labels.get(match.disposition, match.disposition.upper()),
                f"{match.separation_arcsec:.1f}″",
                transit_text,
            )
            for column, value in enumerate(values):
                self.catalog_match_table.setItem(row, column, QTableWidgetItem(value))
        if not self.catalog_matches:
            self.catalog_match_table.setRowCount(0)
            self.lbl_catalog_target.setText(
                f"RA {ra:.6f}° · Dec {dec:+.6f}° · desteklenen kataloglarda 5 arcmin içinde kayıt yok"
            )

    def start_analysis(self):
        if self.worker and self.worker.isRunning():
            return
        if self.target_xy is None:
            QMessageBox.information(self, "Yıldız seçimi", "Önce transit hedef yıldızını seçin.")
            return
        if len(self.wcs_headers) != len(self.fits_files):
            self.start_plate_solve(start_after=True)
            return
        if self.preflight_result is None:
            self.start_photometry_preflight(apply_references=True, start_after=True)
            return
        if not self.preflight_result.analysis_allowed:
            QMessageBox.warning(
                self,
                "Fotometri uygun değil",
                "Hedef yıldız doygun, yetersiz örneklenmiş veya yeterli kararlı referans bulunamadı. "
                "Photometry Suitability bölümünü kontrol edin.",
            )
            return
        self.result = None
        self.measurements.setRowCount(len(self.fits_files))
        self.worker = ExoplanetWorker(
            self.fits_files,
            self.target_xy,
            self.comparison_xy,
            aperture_radius=float(self.combo_aperture.currentText().split()[0]),
            detrend_order={"None": 0, "Linear": 1, "Quadratic": 2}[self.combo_detrend.currentText()],
            recenter=self.chk_recenter.isChecked(),
            linearity_limit_adu=self.preflight_result.linearity.limit_adu,
            master_paths=self.master_paths,
            parent=self,
        )
        self.worker.progress.connect(self.update_progress)
        self.worker.file_progress.connect(self.file_progress)
        self.worker.log.connect(self.log)
        self.worker.measurement_ready.connect(self.measurement_ready)
        self.worker.result_ready.connect(self.result_ready)
        self.worker.finished_scan.connect(self.analysis_finished)
        self.pipeline_progress.reset()
        self.pipeline_progress.set_active_stage("calibration")
        self.lbl_progress_percent.setText("0%")
        self.worker.start()
        self.lbl_verdict.setText("ANALİZ ÇALIŞIYOR")
        self.update_actions()

    def stop_analysis(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.lbl_status.setText("Durdurma bekleniyor")
        if self.plate_worker and self.plate_worker.isRunning():
            self.plate_worker.stop()
            self._start_after_plate_solve = False
            self.lbl_status.setText("Astrometri durduruluyor")
        if self.catalog_worker and self.catalog_worker.isRunning():
            self.catalog_worker.stop()
            self.lbl_status.setText("Katalog güncellemesi durduruluyor")
        if self.preflight_worker and self.preflight_worker.isRunning():
            self.preflight_worker.stop()
            self._start_after_preflight = False
            self.lbl_status.setText("Fotometri ön kontrolü durduruluyor")

    def update_progress(self, value: int, message: str):
        self.progress.setValue(int(value))
        self.lbl_status.setText(message)
        self.lbl_status.setToolTip(message)
        self.lbl_progress_percent.setText(f"{int(value)}%")
        if value < 8:
            stage = "source"
        elif value < 18:
            stage = "calibration"
        elif value < 42:
            stage = "align"
        elif value < 92:
            stage = "photometry"
        elif value < 96:
            stage = "detrend"
        elif value < 100:
            stage = "model"
        else:
            stage = "review"
        if self.pipeline_progress.active_key != stage:
            self.pipeline_progress.set_active_stage(stage)
        self.progress.setValue(int(value))
        if value >= 100:
            self.pipeline_progress.complete()

    def file_progress(self, current, total, name, stage):
        value = int(current / max(total, 1) * 100)
        self.pipeline_progress.set_file_value(value)
        text = f"{stage}: {current}/{total} · {name}"
        self.lbl_file_status.setText(text)
        self.lbl_file_status.setToolTip(text)

    def measurement_ready(self, index, target, comparison):
        values = (str(index + 1), f"{target:.1f}", f"{comparison:.1f}")
        for column, value in enumerate(values):
            self.measurements.setItem(index, column, QTableWidgetItem(value))

    def result_ready(self, result):
        self.result = result
        verdict = "TRANSİT ADAYI" if result.transit_candidate else "BELİRGİN TRANSİT YOK"
        self.lbl_verdict.setText(verdict)
        self.lbl_verdict.setObjectName("badgeCandidate" if result.transit_candidate else "badgeKnown")
        self.lbl_verdict.style().unpolish(self.lbl_verdict)
        self.lbl_verdict.style().polish(self.lbl_verdict)
        self.result_table.set_rows([
            ("Depth", f"{result.depth * 100:.4f}%", "warn" if result.transit_candidate else None),
            ("Depth uncertainty", f"{result.depth_uncertainty * 100:.4f}%", None),
            ("Detection SNR", f"{result.detection_snr:.2f}", "ok" if result.detection_snr >= 5 else "warn"),
            ("Duration", f"{result.duration_minutes:.2f} min" if result.duration_minutes else "-", None),
            ("Mid transit", f"{result.mid_transit_jd:.7f}" if result.mid_transit_jd else "-", None),
            ("Scatter", f"{result.scatter * 100:.4f}%", None),
            ("Time system", result.time_system, "warn" if "UTC" in result.time_system else "ok"),
            ("Rp / R*", f"{result.model_fit.radius_ratio:.4f}" if result.model_fit and result.model_fit.success else "-", None),
            ("Impact b", f"{result.model_fit.impact_parameter:.3f}" if result.model_fit and result.model_fit.success else "-", None),
            ("Delta BIC", f"{result.model_fit.delta_bic:.2f}" if result.model_fit and result.model_fit.success else "-", "ok" if result.model_fit and (result.model_fit.delta_bic or 0) >= 6 else "warn"),
        ])
        weights = ", ".join(f"C{index + 1} {weight * 100:.0f}%" for index, weight in enumerate(result.comparison_weights))
        flags = ", ".join(result.quality_flags) if result.quality_flags else "none"
        self.lbl_quality.setText(f"Comparison weights: {weights or '-'}\nQuality flags: {flags}")
        self.plot_result()

    def plot_result(self):
        if self.result is None:
            return
        result = self.result
        times = np.asarray(result.times_jd, dtype=float)
        minutes = (times - times[0]) * 1440.0
        detrended = self.combo_curve.currentText() == "Detrended"
        flux = np.asarray(result.relative_flux if detrended else result.raw_relative_flux, dtype=float)
        uncertainty = np.asarray(result.flux_uncertainty, dtype=float)
        valid = np.asarray(result.valid_mask, dtype=bool)
        self.curve_plot.clear()
        self.curve_plot.addItem(
            pg.ErrorBarItem(
                x=minutes[valid],
                y=flux[valid],
                height=2.0 * uncertainty[valid],
                beam=max(float(np.ptp(minutes)) * 0.002, 0.02),
                pen=pg.mkPen("#526b73", width=0.8),
            )
        )
        self.curve_plot.plot(
            minutes[valid],
            flux[valid],
            pen=pg.mkPen("#2dd7df", width=1.2),
            symbol="o",
            symbolSize=5,
            symbolBrush="#2dd7df",
        )
        if np.any(~valid):
            self.curve_plot.plot(minutes[~valid], flux[~valid], pen=None, symbol="x", symbolSize=7, symbolPen="#ef6257")
        if detrended and result.model_fit and result.model_fit.success:
            model = np.asarray(result.model_fit.model_flux, dtype=float)
            order = np.argsort(minutes)
            self.curve_plot.plot(minutes[order], model[order], pen=pg.mkPen("#efbc3f", width=2.0))
            self.lbl_curve_model.setText(f"Limb-darkened fit · ΔBIC {result.model_fit.delta_bic:.1f}")
        else:
            self.lbl_curve_model.setText("Detrended model" if not detrended else "Physical model unavailable")
        self.curve_plot.addLine(y=1.0, pen=pg.mkPen("#65757d", width=1, style=Qt.PenStyle.DashLine))

    def analysis_finished(self):
        self.lbl_status.setText("Transit analizi tamamlandı" if self.result else "Analiz sonuç üretmeden durdu")
        if self.result:
            self.pipeline_progress.complete()
            self.lbl_progress_percent.setText("100%")
        elif self.pipeline_progress.active_key:
            self.pipeline_progress.set_stage_state(self.pipeline_progress.active_key, "warning")
        self.update_actions()

    def export_csv(self):
        if self.result is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Işık eğrisini aktar", self.current_folder, "CSV (*.csv)")
        if not path:
            return
        if not path.lower().endswith(".csv"):
            path += ".csv"
        with open(path, "w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    self.result.time_system,
                    "raw_relative_flux",
                    "detrended_flux",
                    "flux_uncertainty",
                    "target_flux",
                    "comparison_ensemble",
                    "valid",
                    "model_flux",
                ]
            )
            model = self.result.model_fit.model_flux if self.result.model_fit else [float("nan")] * len(self.result.times_jd)
            writer.writerows(
                zip(
                    self.result.times_jd,
                    self.result.raw_relative_flux,
                    self.result.relative_flux,
                    self.result.flux_uncertainty,
                    self.result.target_flux,
                    self.result.comparison_flux,
                    self.result.valid_mask,
                    model,
                )
            )
        self.log(f"INFO|Export|Light curve exported: {Path(path).name}")

    def update_actions(self):
        running = bool(self.worker and self.worker.isRunning())
        solving = bool(self.plate_worker and self.plate_worker.isRunning())
        preflighting = bool(self.preflight_worker and self.preflight_worker.isRunning())
        updating_catalog = bool(self.catalog_worker and self.catalog_worker.isRunning())
        ready = len(self.fits_files) >= 5 and self.target_xy is not None and bool(self.comparison_xy)
        if hasattr(self, "btn_start"):
            self.btn_open.setEnabled(not running and not solving and not preflighting)
            self.btn_start.setEnabled(ready and not running and not solving and not preflighting)
            self.btn_plate_solve.setEnabled(bool(self.fits_files) and not running and not solving and not preflighting)
            self.btn_auto_references.setEnabled(self.target_xy is not None and not running and not solving and not preflighting)
            self.btn_catalog.setEnabled(not updating_catalog)
            self.btn_update_catalog.setEnabled(not updating_catalog)
            self.btn_stop.setEnabled(running or solving or updating_catalog or preflighting)
            self.btn_export.setEnabled(self.result is not None)

    def log(self, payload: str):
        self.analysis_log.append_message(payload)

    def shutdown_workers(self):
        for worker in (self.preview_worker, self.worker, self.plate_worker, self.preflight_worker, self.catalog_worker):
            if worker and worker.isRunning():
                if hasattr(worker, "stop"):
                    worker.stop()
                worker.wait(3000)


def find_fits_files(folder: str) -> list[str]:
    suffixes = (".fit", ".fits", ".fts", ".fits.gz")
    root = Path(folder)
    direct = [path for path in root.iterdir() if path.is_file() and path.name.lower().endswith(suffixes)]
    files = direct if direct else [path for path in root.rglob("*") if path.is_file() and path.name.lower().endswith(suffixes)]
    return sorted((str(path) for path in files), key=lambda value: Path(value).name.lower())
