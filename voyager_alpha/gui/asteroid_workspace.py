from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from PyQt6.QtCore import QTimer, Qt
from PyQt6.QtGui import QDoubleValidator, QIntValidator, QPixmap
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QStyle,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ..core.ades import render_ades_psv_draft
from ..core.fits_io import read_fits_image
from ..core.metadata import inspect_sequence
from ..core.pipeline import AsteroidWorker, materialize_sequence_frame
from ..core.reporting import render_tracklet_html_report
from .viewer import FitsViewer, sequence_stf_levels, stretch_image
from .widgets import (
    AnalysisLog,
    Filmstrip,
    KeyValueTable,
    PipelineProgressStrip,
    SegmentedControl,
    WorkflowRail,
    camera_metadata_rows,
    metric_row,
    section,
)
from .workers import DetectionOptimizationWorker, PlateSolveWorker, PreviewCacheWorker, SyntheticTrackWorker


class AsteroidWorkspace(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("workspaceRoot")
        self.fits_files: list[str] = []
        self.sequence_frames = []
        self.sequence_result = None
        self.frame_candidates: dict[int, list] = {}
        self.playback_cache: dict[int, tuple[np.ndarray, tuple[int, int]]] = {}
        self.difference_cache: dict[int, tuple[np.ndarray, tuple[int, int]]] = {}
        self.preview_worker = None
        self.preview_cache_complete = False
        self.preview_alignment_stats: dict[int, tuple[str, float, int]] = {}
        self.preview_generation = 0
        self.analysis_worker = None
        self.plate_worker = None
        self.synthetic_worker = None
        self.optimization_worker = None
        self.current_index = 0
        self.current_folder = ""
        self.current_tracklet = None
        self.current_known = None
        self.sequence_levels = None
        self.residual_levels = None
        self.manual_levels = (0.1, 0.35, 99.9)
        self.synthetic_result = None
        self.review_history = []
        self.analysis_mode = "guided"
        self.master_paths = {"bias": "", "dark": "", "flat": ""}
        self.custom_detection = {
            "sigma": 5.0,
            "min_pixels": 5,
            "min_frames": 3,
            "max_sources": 24,
            "edge_margin": 6,
            "expected_fwhm_px": 3.0,
            "min_motion_px_per_frame": 1.5,
            "max_step_px": 35.0,
            "min_median_snr": 5.1,
            "max_fit_rms_px": 1.8,
            "strong_fit_rms_px": 0.9,
            "match_tolerance_px": 2.8,
            "min_track_occupancy": 0.5,
            "max_missing_gap_frames": 3,
            "max_artifact_fraction": 0.34,
            "persistence_fraction": 0.18,
        }
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.next_frame)
        self._build_ui()
        self._install_tooltips()
        self.workflow.reset()
        self.update_action_state()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(self._build_toolbar())

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(1)
        self.left_panel = self._build_left_panel()
        self.center_panel = self._build_center_panel()
        self.right_panel = self._build_right_panel()
        splitter.addWidget(self.left_panel)
        splitter.addWidget(self.center_panel)
        splitter.addWidget(self.right_panel)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 0)
        splitter.setSizes([285, 1040, 330])
        root.addWidget(splitter, 1)
        root.addWidget(self._build_status_bar())

    def _build_toolbar(self):
        frame = QFrame()
        frame.setObjectName("toolbar")
        layout = QHBoxLayout(frame)
        layout.setContentsMargins(9, 7, 9, 7)
        layout.setSpacing(7)
        self.btn_open = QPushButton("Kaynak Klasör")
        self.btn_open.setObjectName("blueButton")
        self.btn_open.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon))
        self.btn_open.clicked.connect(self.load_folder)
        self.btn_open.setMinimumWidth(116)
        layout.addWidget(self.btn_open)
        self.btn_rescan = QPushButton("Yeniden Tara")
        self.btn_rescan.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_BrowserReload))
        self.btn_rescan.clicked.connect(self.rescan)
        self.btn_rescan.setMinimumWidth(108)
        layout.addWidget(self.btn_rescan)
        self.btn_start = QPushButton("Guided Analizi Başlat")
        self.btn_start.setObjectName("primaryButton")
        self.btn_start.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        self.btn_start.clicked.connect(self.start_analysis)
        self.btn_start.setMinimumWidth(166)
        layout.addWidget(self.btn_start)
        self.btn_stop = QPushButton("Durdur")
        self.btn_stop.setObjectName("dangerButton")
        self.btn_stop.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaStop))
        self.btn_stop.clicked.connect(self.stop_analysis)
        self.btn_stop.setMinimumWidth(78)
        layout.addWidget(self.btn_stop)
        layout.addSpacing(8)
        self.flow_control = SegmentedControl([("guided", "Guided Flow"), ("manual", "Manual Flow")])
        self.flow_control.changed.connect(self.set_flow_mode)
        layout.addWidget(self.flow_control)
        layout.addStretch(1)
        self.lbl_session = QLabel("Aktif sekans yok")
        self.lbl_session.setObjectName("muted")
        layout.addWidget(self.lbl_session)
        self.lbl_session_state = QLabel("READY")
        self.lbl_session_state.setObjectName("badgeNeutral")
        layout.addWidget(self.lbl_session_state)
        return frame

    def _build_left_panel(self):
        panel = QFrame()
        panel.setObjectName("leftPanel")
        panel.setMinimumWidth(260)
        panel.setMaximumWidth(320)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(5)

        head = QHBoxLayout()
        title = QLabel("Sequence Validation")
        title.setObjectName("panelTitle")
        self.lbl_sequence_state = QLabel("● Waiting")
        self.lbl_sequence_state.setObjectName("muted")
        head.addWidget(title)
        head.addStretch(1)
        head.addWidget(self.lbl_sequence_state)
        layout.addLayout(head)
        self.sequence_table = KeyValueTable()
        self.sequence_table.set_auto_height(True)
        layout.addWidget(self.sequence_table)
        self.workflow = WorkflowRail()
        self.workflow.step_clicked.connect(self.workflow_step_clicked)
        self.workflow.setVisible(False)
        layout.addWidget(self.workflow)

        calibration, calibration_layout = section("Calibration Masters")
        self.master_labels = {}
        for row, (key, label) in enumerate((("bias", "Master Bias"), ("dark", "Master Dark"), ("flat", "Master Flat"))):
            line = QHBoxLayout()
            button = QPushButton(label)
            button.setFixedWidth(102)
            button.setToolTip(f"İsteğe bağlı {label.lower()} FITS kalibrasyon karesini seçer.")
            button.clicked.connect(lambda _checked=False, value=key: self.choose_master(value))
            value_label = QLabel("Seçilmedi")
            value_label.setObjectName("muted")
            line.addWidget(button)
            line.addWidget(value_label, 1)
            calibration_layout.addLayout(line)
            self.master_labels[key] = value_label
        layout.addWidget(calibration)

        detection, detection_layout = section("Detection Profile")
        self.combo_profile = QComboBox()
        self.combo_profile.addItems(["Documented Workflow", "Conservative", "Balanced", "Auto Optimized", "Custom"])
        detection_layout.addWidget(self.combo_profile)
        profile_actions = QHBoxLayout()
        profile_actions.setSpacing(5)
        self.btn_optimize_detection = QPushButton("Otomatik Ayar")
        self.btn_optimize_detection.setObjectName("blueButton")
        self.btn_optimize_detection.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_DialogApplyButton))
        self.btn_optimize_detection.clicked.connect(self.start_detection_optimization)
        self.btn_detection_settings = QPushButton("Gelişmiş")
        self.btn_detection_settings.clicked.connect(self.edit_detection_settings)
        profile_actions.addWidget(self.btn_optimize_detection)
        profile_actions.addWidget(self.btn_detection_settings)
        detection_layout.addLayout(profile_actions)
        self.lbl_optimization = QLabel("Otomatik profil henüz ölçülmedi")
        self.lbl_optimization.setObjectName("muted")
        detection_layout.addWidget(self.lbl_optimization)
        self.chk_known = QCheckBox("Bilinen cisimler")
        self.chk_known.setChecked(True)
        self.chk_candidates = QCheckBox("Residual adaylar")
        self.chk_candidates.setChecked(True)
        self.chk_track = QCheckBox("Seçili tracklet izi")
        self.chk_track.setChecked(True)
        self.chk_grid = QCheckBox("Gökyüzü ızgarası")
        for checkbox in (self.chk_known, self.chk_candidates, self.chk_track):
            checkbox.stateChanged.connect(self.refresh_overlays)
        self.chk_grid.stateChanged.connect(lambda state: self.viewer.set_grid_visible(bool(state)))
        detection_layout.addWidget(self.chk_known)
        detection_layout.addWidget(self.chk_candidates)
        detection_layout.addWidget(self.chk_track)
        detection_layout.addWidget(self.chk_grid)
        layout.addWidget(detection)

        self.manual_panel, manual_layout = section("Manual Commands")
        self.btn_validate = QPushButton("Sekansı Doğrula")
        self.btn_validate.clicked.connect(self.validate_loaded_sequence)
        self.btn_solve = QPushButton("Seçili Kareyi Plate Solve")
        self.btn_solve.clicked.connect(self.solve_selected_frame)
        self.btn_known_run = QPushButton("Generate / Recover Known")
        self.btn_known_run.clicked.connect(self.start_known_recovery)
        self.btn_discover_run = QPushButton("Discover Unknown Movers")
        self.btn_discover_run.clicked.connect(self.start_unknown_discovery)
        for button in (self.btn_validate, self.btn_solve, self.btn_known_run, self.btn_discover_run):
            manual_layout.addWidget(button)
        self.manual_panel.setVisible(False)
        layout.addWidget(self.manual_panel)
        layout.addStretch(1)
        return panel

    def _build_center_panel(self):
        panel = QFrame()
        panel.setObjectName("centerPanel")
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(7, 7, 7, 6)
        layout.setSpacing(6)

        toolbar = QFrame()
        toolbar.setObjectName("viewerToolbar")
        tools = QHBoxLayout(toolbar)
        tools.setContentsMargins(6, 4, 6, 4)
        tools.setSpacing(5)
        tools.addWidget(QLabel("View"))
        self.combo_view = QComboBox()
        self.combo_view.addItems(["Blink Review", "Single Frame", "Tracklet Follow"])
        self.combo_view.setMinimumWidth(112)
        tools.addWidget(self.combo_view)
        tools.addWidget(QLabel("Mode"))
        self.combo_mode = QComboBox()
        self.combo_mode.addItems(["Original", "Difference", "Synthetic Track"])
        self.combo_mode.setMinimumWidth(108)
        self.combo_mode.currentTextChanged.connect(self.refresh_current_frame)
        tools.addWidget(self.combo_mode)
        tools.addWidget(QLabel("Stretch"))
        self.combo_stretch = QComboBox()
        self.combo_stretch.addItems(["Auto STF", "Asinh", "Linear", "Manual STF"])
        self.combo_stretch.setMinimumWidth(102)
        self.combo_stretch.currentTextChanged.connect(self.stretch_changed)
        tools.addWidget(self.combo_stretch)
        self.chk_invert = QCheckBox("Invert")
        self.chk_invert.stateChanged.connect(self.refresh_current_frame)
        tools.addWidget(self.chk_invert)
        tools.addStretch(1)
        self.lbl_fov = QLabel("FOV: waiting")
        self.lbl_fov.setObjectName("muted")
        tools.addWidget(self.lbl_fov)
        for icon, callback, tip in (
            (QStyle.StandardPixmap.SP_ArrowUp, lambda: self.viewer.zoom_by(0.75), "Yakınlaştır"),
            (QStyle.StandardPixmap.SP_ArrowDown, lambda: self.viewer.zoom_by(1.25), "Uzaklaştır"),
            (QStyle.StandardPixmap.SP_FileDialogDetailedView, self.viewer_auto_fit, "Görüntüyü alana sığdır"),
        ):
            button = QToolButton()
            button.setIcon(self.style().standardIcon(icon))
            button.setToolTip(tip)
            button.clicked.connect(callback)
            tools.addWidget(button)
        layout.addWidget(toolbar)

        self.viewer = FitsViewer()
        layout.addWidget(self.viewer, 1)

        film_row = QHBoxLayout()
        self.filmstrip = Filmstrip(page_size=6)
        self.filmstrip.frame_selected.connect(self.change_frame)
        film_row.addWidget(self.filmstrip, 1)
        blink_controls = QFrame()
        blink_controls.setObjectName("filmstripBar")
        blink_layout = QGridLayout(blink_controls)
        blink_layout.setContentsMargins(8, 6, 8, 6)
        blink_layout.addWidget(QLabel("Blink Speed"), 0, 0, 1, 3)
        self.combo_blink_speed = QComboBox()
        self.combo_blink_speed.addItems(["1 fps", "2 fps", "4 fps", "6 fps", "10 fps", "15 fps", "20 fps"])
        self.combo_blink_speed.setCurrentText("6 fps")
        self.combo_blink_speed.currentTextChanged.connect(self.update_blink_speed)
        blink_layout.addWidget(self.combo_blink_speed, 1, 0, 1, 3)
        self.btn_prev = QToolButton()
        self.btn_prev.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaSkipBackward))
        self.btn_prev.clicked.connect(self.previous_frame)
        self.btn_play = QToolButton()
        self.btn_play.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        self.btn_play.clicked.connect(self.toggle_play)
        self.btn_next = QToolButton()
        self.btn_next.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaSkipForward))
        self.btn_next.clicked.connect(self.next_frame)
        blink_layout.addWidget(self.btn_prev, 2, 0)
        blink_layout.addWidget(self.btn_play, 2, 1)
        blink_layout.addWidget(self.btn_next, 2, 2)
        self.lbl_frame_counter = QLabel("0 / 0")
        self.lbl_frame_counter.setAlignment(Qt.AlignmentFlag.AlignCenter)
        blink_layout.addWidget(self.lbl_frame_counter, 3, 0, 1, 3)
        film_row.addWidget(blink_controls)
        layout.addLayout(film_row)

        results_split = QSplitter(Qt.Orientation.Horizontal)
        log_frame, log_layout = section("Analysis Log")
        self.log_table = AnalysisLog()
        log_layout.addWidget(self.log_table)
        results_split.addWidget(log_frame)
        result_frame, result_layout = section("Session Results")
        self.result_tabs = QTabWidget()
        self.tracklet_table = self._make_tracklet_table()
        self.known_table = self._make_known_table()
        self.result_tabs.addTab(self.tracklet_table, "Discover Results")
        self.result_tabs.addTab(self.known_table, "Generate / Recovery")
        result_layout.addWidget(self.result_tabs)
        results_split.addWidget(result_frame)
        results_split.setSizes([500, 700])
        results_split.setFixedHeight(230)
        layout.addWidget(results_split)
        return panel

    def _build_right_panel(self):
        panel = QFrame()
        panel.setObjectName("rightPanel")
        panel.setMinimumWidth(300)
        panel.setMaximumWidth(390)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 7, 8, 7)
        layout.setSpacing(6)
        title = QLabel("Tracklet Evidence")
        title.setObjectName("panelTitle")
        layout.addWidget(title)
        self.lbl_evidence_id = QLabel("Tracklet seçilmedi")
        self.lbl_evidence_id.setStyleSheet("font-size:18px;font-weight:600;color:#eef3f5;")
        layout.addWidget(self.lbl_evidence_id)
        self.lbl_verdict = QLabel("HAREKETLİ CİSİM KANITI BEKLENİYOR")
        self.lbl_verdict.setObjectName("badgeNeutral")
        self.lbl_verdict.setWordWrap(True)
        layout.addWidget(self.lbl_verdict)

        evidence, evidence_layout = section("Measurements")
        evidence_layout.setSpacing(2)
        self.evidence_values = {}
        for key, label in (
            ("frames", "Frames"),
            ("arc", "Time arc"),
            ("motion", "Motion"),
            ("pa", "Position angle"),
            ("rms", "Fit RMS"),
            ("snr", "Median SNR"),
            ("confidence", "Triage score"),
        ):
            row, value = metric_row(label)
            evidence_layout.addWidget(row)
            self.evidence_values[key] = value
        layout.addWidget(evidence)

        crossmatch, cross_layout = section("Cross-match & Quality")
        cross_layout.setSpacing(3)
        row, self.lbl_skybot = metric_row("SkyBoT", "Waiting")
        cross_layout.addWidget(row)
        row, self.lbl_artifacts = metric_row("Artifact flags", "None")
        cross_layout.addWidget(row)
        row, self.lbl_review_state = metric_row("Review", "Unreviewed")
        cross_layout.addWidget(row)
        layout.addWidget(crossmatch)

        review, review_layout = section("Review")
        buttons = QHBoxLayout()
        self.btn_accept = QPushButton("✓  Accept")
        self.btn_accept.setObjectName("reviewAccept")
        self.btn_accept.clicked.connect(lambda: self.set_review_status("accepted"))
        self.btn_reject = QPushButton("✕  Reject")
        self.btn_reject.setObjectName("reviewReject")
        self.btn_reject.clicked.connect(lambda: self.set_review_status("rejected"))
        buttons.addWidget(self.btn_accept)
        buttons.addWidget(self.btn_reject)
        review_layout.addLayout(buttons)
        self.btn_follow = QPushButton("?  Needs follow-up")
        self.btn_follow.setObjectName("reviewFollow")
        self.btn_follow.clicked.connect(lambda: self.set_review_status("follow_up"))
        review_layout.addWidget(self.btn_follow)
        self.btn_undo_review = QPushButton("Son kararı geri al")
        self.btn_undo_review.clicked.connect(self.undo_review)
        review_layout.addWidget(self.btn_undo_review)
        self.notes = QPlainTextEdit()
        self.notes.setPlaceholderText("İnceleme notu...")
        self.notes.setFixedHeight(42)
        self.notes.textChanged.connect(self.notes_changed)
        review_layout.addWidget(self.notes)
        layout.addWidget(review)

        actions, action_layout = section("Confirmation & Export")
        confirm_buttons = QHBoxLayout()
        confirm_buttons.setSpacing(5)
        self.btn_center = QPushButton("Seçili Objeyi Ortala")
        self.btn_center.clicked.connect(self.center_selected_tracklet)
        self.btn_synthetic = QPushButton("Synthetic Track")
        self.btn_synthetic.setObjectName("primaryButton")
        self.btn_synthetic.clicked.connect(self.run_synthetic_track)
        confirm_buttons.addWidget(self.btn_center)
        confirm_buttons.addWidget(self.btn_synthetic)
        action_layout.addLayout(confirm_buttons)
        self.lbl_synthetic = QLabel("Sönük aday doğrulaması çalıştırılmadı")
        self.lbl_synthetic.setObjectName("muted")
        self.lbl_synthetic.setToolTip("Synthetic Track sonucu burada gösterilir.")
        action_layout.addWidget(self.lbl_synthetic)

        export_buttons = QHBoxLayout()
        export_buttons.setSpacing(5)
        self.btn_export_csv = QPushButton("CSV")
        self.btn_export_csv.clicked.connect(self.export_csv)
        self.btn_export_report = QPushButton("HTML")
        self.btn_export_report.clicked.connect(self.export_report)
        self.btn_export_ades = QPushButton("ADES")
        self.btn_export_ades.clicked.connect(self.export_ades)
        export_buttons.addWidget(self.btn_export_csv)
        export_buttons.addWidget(self.btn_export_report)
        export_buttons.addWidget(self.btn_export_ades)
        action_layout.addLayout(export_buttons)
        self.btn_export_ades.setToolTip(
            "Yalnızca Accept edilmiş ve WCS doğrulanmış ölçümleri inceleme gerektiren ADES taslağına aktarır."
        )
        layout.addWidget(actions)
        layout.addStretch(1)
        return panel

    def _build_status_bar(self):
        frame = QFrame()
        frame.setObjectName("statusBar")
        layout = QHBoxLayout(frame)
        layout.setContentsMargins(8, 4, 8, 4)
        layout.setSpacing(8)
        self.lbl_status = QLabel("Hazır")
        self.lbl_status.setObjectName("muted")
        self.lbl_status.setMinimumWidth(140)
        self.lbl_status.setMaximumWidth(220)
        self.lbl_status.setToolTip("Aktif işlem aşaması")
        self.pipeline_progress = PipelineProgressStrip(
            [
                ("source", "Kaynak"),
                ("calibration", "Kalibrasyon"),
                ("solve", "WCS"),
                ("reference", "Model"),
                ("detect", "Tespit"),
                ("link", "Tracklet"),
                ("known", "Katalog"),
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

    def _make_tracklet_table(self):
        columns = ["ID", "Class", "Frames", "Arc", "Motion", "PA", "Fit RMS", "SNR", "Score", "Review"]
        table = QTableWidget(0, len(columns))
        table.setHorizontalHeaderLabels(columns)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        table.setAlternatingRowColors(True)
        table.setWordWrap(False)
        table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        table.verticalHeader().setDefaultSectionSize(24)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        table.currentCellChanged.connect(self.tracklet_selected)
        return table

    def _make_known_table(self):
        columns = ["Object", "Type", "Pred Mag", "Local SNR", "Offset", "Visible", "Confidence"]
        table = QTableWidget(0, len(columns))
        table.setHorizontalHeaderLabels(columns)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        table.setAlternatingRowColors(True)
        table.setWordWrap(False)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        table.currentCellChanged.connect(self.known_selected)
        return table

    def _install_tooltips(self):
        tips = {
            self.btn_open: "Zaman sıralı science FITS karelerinin bulunduğu klasörü açar.",
            self.btn_rescan: "Aynı klasörü yeniden tarar ve önceki oturum sonucunu temizler.",
            self.btn_start: "Kalibrasyon, astrometri, hizalama, bilinen cisim arama ve bilinmeyen hareket taramasını sırayla çalıştırır.",
            self.btn_stop: "Çalışan analizi kontrollü biçimde durdurur; kaynak dosyaları değiştirmez.",
            self.flow_control: "Guided Flow güvenli bilimsel sırayı uygular; Manual Flow araçları bağımsız çalıştırmanızı sağlar.",
            self.combo_profile: "Conservative profil yanlış pozitifleri azaltır; Custom değerler Advanced Settings içindedir.",
            self.btn_optimize_detection: "Örnek karelerde gürültü, FWHM, yıldız yoğunluğu, WCS görüntü ölçeği, kadans ve residual yoğunluğunu ölçerek tespit profilini otomatik ayarlar.",
            self.btn_detection_settings: "Kaynak eşiği, minimum alan, tracklet kare sayısı ve kare başına aday sınırını düzenler.",
            self.combo_blink_speed: "Önbelleğe alınmış sekansı kalite kaybetmeden seçilen hızda oynatır.",
            self.combo_stretch: "Yalnızca ekran görünümünü değiştirir; bilimsel piksel değerlerini değiştirmez.",
            self.chk_grid: "İnceleme için düşük yoğunluklu gökyüzü ızgarasını açar.",
            self.btn_validate: "DATE-OBS, pozlama, WCS ve kare boyutu tutarlılığını yeniden denetler.",
            self.btn_solve: "Yalnızca seçili FITS karesini ASTAP ile geçici bir kopya üzerinde plate solve eder.",
            self.btn_known_run: "WCS alanındaki katalog tahminlerini ve yerel görüntü karşılıklarını öncelikli olarak inceler.",
            self.btn_discover_run: "Katalog sorgusunu atlayıp bilinmeyen hareket zinciri taramasını çalıştırır.",
            self.btn_prev: "Önceki FITS karesine gider.",
            self.btn_play: "Blink oynatımını başlatır veya duraklatır.",
            self.btn_next: "Sonraki FITS karesine gider.",
            self.btn_accept: "Seçili tracklet'i insan incelemesinden geçmiş aday olarak işaretler.",
            self.btn_reject: "Seçili tracklet'i reddeder; ölçümler oturumda korunur.",
            self.btn_follow: "Seçili tracklet için ek gözlem veya yeniden inceleme gerektiğini işaretler.",
            self.btn_undo_review: "Son Accept, Reject veya Needs follow-up kararını geri alır.",
            self.btn_center: "Görüntüyü seçili tracklet'in mevcut konumuna yakınlaştırır.",
            self.btn_synthetic: "Seçili tracklet hızına göre kareleri kaydırıp istifleyerek sönük hareketli sinyali doğrular.",
            self.btn_export_csv: "Tüm tracklet ölçümlerini ve inceleme durumlarını CSV olarak kaydeder.",
            self.btn_export_report: "Sekans kalitesi ve tracklet kanıtlarını HTML rapora aktarır.",
            self.btn_export_ades: "MPC gönderimine hazır olmayan, manuel doğrulama gerektiren ADES taslağı üretir.",
        }
        for widget, text in tips.items():
            widget.setToolTip(text)

    def set_flow_mode(self, mode: str):
        manual = mode == "manual"
        self.manual_panel.setVisible(manual)
        self.workflow.setVisible(False)
        self.btn_start.setText("Tam Pipeline" if manual else "Guided Analizi Başlat")
        self.log(f"INFO|Workflow|{'Manual Flow' if manual else 'Guided Flow'} etkin")

    def load_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "FITS sekans klasörü seç", self.current_folder or "")
        if folder:
            self.scan_folder(folder)

    def rescan(self):
        if self.current_folder:
            self.scan_folder(self.current_folder)

    def scan_folder(self, folder: str):
        self.shutdown_workers()
        self.timer.stop()
        self.current_folder = folder
        self.fits_files = find_fits_files(folder)
        self.sequence_result = None
        self.sequence_frames = []
        self.frame_candidates.clear()
        self.playback_cache.clear()
        self.preview_cache_complete = False
        self.preview_alignment_stats.clear()
        self.difference_cache.clear()
        self.synthetic_result = None
        self.current_tracklet = None
        self.current_known = None
        self.current_index = 0
        self.tracklet_table.setRowCount(0)
        self.known_table.setRowCount(0)
        self.filmstrip.set_frames(self.fits_files)
        self.lbl_session.setText(Path(folder).name if self.fits_files else "FITS bulunamadı")
        if not self.fits_files:
            self.log(f"WARN|Source|FITS bulunamadı: {folder}")
            self.viewer.clear_image()
            self.update_validation([])
            self.update_action_state()
            return
        self.log(f"INFO|Source|{len(self.fits_files)} FITS karesi yüklendi")
        self.validate_loaded_sequence()
        try:
            first = np.asarray(read_fits_image(self.fits_files[0]), dtype=np.float32)
            self.sequence_levels = sequence_stf_levels(first)
            self.viewer.set_image(first, stretch_mode=self.combo_stretch.currentText(), sequence_levels=None, keep_view=False)
        except Exception as exc:
            self.log(f"ERROR|Viewer|İlk FITS görüntülenemedi: {exc}")
        self._start_preview_cache()
        self.update_action_state()

    def _start_preview_cache(self, sequence_result=None):
        if self.preview_worker and self.preview_worker.isRunning():
            self.preview_worker.stop()
            self.preview_worker.wait(5000)
        self.timer.stop()
        self.playback_cache.clear()
        self.preview_cache_complete = False
        self.preview_alignment_stats.clear()
        self.preview_generation += 1
        generation = self.preview_generation
        self.preview_worker = PreviewCacheWorker(
            self.fits_files,
            sequence_result=sequence_result,
            parent=self,
        )
        self.preview_worker.preview_ready.connect(
            lambda index, display, shape, thumbnail, token=generation: self.preview_ready(
                index, display, shape, thumbnail, token
            )
        )
        self.preview_worker.alignment_ready.connect(
            lambda index, method, rms, count, token=generation: self.preview_alignment_ready(
                index, method, rms, count, token
            )
        )
        self.preview_worker.failed.connect(
            lambda message, token=generation: self.log(f"ERROR|Preview|{message}")
            if token == self.preview_generation
            else None
        )
        self.preview_worker.complete.connect(
            lambda token=generation: self.preview_cache_finished(token)
        )
        self.preview_worker.start()
        self.lbl_status.setText("Yıldız hizalı Auto STF blink önbelleği hazırlanıyor")
        self.update_action_state()

    def preview_alignment_ready(self, index, method, rms, matched_stars, generation=None):
        if generation is not None and generation != self.preview_generation:
            return
        self.preview_alignment_stats[int(index)] = (str(method), float(rms), int(matched_stars))
        if "rot180" in method:
            self.log(
                f"INFO|Blink Align|Kare {int(index) + 1}: meridian flip algılandı ve 180 derece düzeltildi "
                f"({int(matched_stars)} yıldız, RMS {float(rms):.2f} px)"
            )
        if not method.startswith("star") and method != "reference":
            self.log(
                f"WARN|Blink Align|Kare {int(index) + 1}: yeterli yıldız eşleşmedi; "
                "bu kare yalnızca faz korelasyonuyla hizalandı"
            )
        elif method.startswith("star") and float(rms) > 1.5:
            self.log(f"WARN|Blink Align|Kare {int(index) + 1}: yıldız hizalama RMS {float(rms):.2f} px")

    def preview_cache_finished(self, generation=None):
        if generation is not None and generation != self.preview_generation:
            return
        self.preview_cache_complete = len(self.playback_cache) == len(self.fits_files) and bool(self.fits_files)
        star_aligned = sum(method.startswith("star") for method, _rms, _count in self.preview_alignment_stats.values())
        if self.preview_cache_complete:
            self.log(
                f"INFO|Blink|{len(self.playback_cache)} kare hizalandı; her kareye bağımsız Auto STF uygulandı "
                f"({star_aligned} yıldız çözümü)"
            )
            self.lbl_status.setText("Yıldız hizalı Auto STF blink hazır")
        else:
            self.log(
                f"ERROR|Blink|Önbellek eksik: {len(self.playback_cache)}/{len(self.fits_files)} kare; blink başlatılmadı"
            )
            self.lbl_status.setText("Blink önbelleği tamamlanamadı")
        self.update_action_state()

    def validate_loaded_sequence(self):
        if not self.fits_files:
            return
        try:
            self.sequence_frames = inspect_sequence(self.fits_files)
            self.fits_files = [frame.file_path for frame in self.sequence_frames]
            self.filmstrip.set_frames(self.fits_files)
            self.update_validation(self.sequence_frames)
            self.workflow.set_state("source", "done")
            self.workflow.set_state("calibration", "active")
            self.log(f"INFO|Metadata|Sekans doğrulandı: {len(self.sequence_frames)} kare")
        except Exception as exc:
            self.log(f"ERROR|Metadata|{exc}")

    def update_validation(self, frames):
        if not frames:
            self.sequence_table.set_rows([("FITS Frames", "0", "warn"), ("Sequence", "Waiting", None)])
            self.lbl_sequence_state.setText("● Waiting")
            return
        total = len(frames)
        valid_time = sum(frame.midpoint_jd is not None for frame in frames)
        wcs_count = sum(frame.has_wcs for frame in frames)
        exposure_values = [frame.exposure_seconds for frame in frames if frame.exposure_seconds is not None]
        exposure = f"{np.median(exposure_values):.1f} s" if exposure_values else "Missing"
        cadence = "--"
        times = [frame.midpoint_jd for frame in frames if frame.midpoint_jd is not None]
        if len(times) > 1:
            cadence = f"{np.median(np.diff(times)) * 86400.0:.1f} s"
        flags = sum(len(frame.quality_flags) for frame in frames)
        self.sequence_table.set_rows(
            [
                ("FITS Frames", str(total), "cyan"),
                ("Time / WCS", f"{valid_time}/{total} · {wcs_count}/{total}", "ok" if valid_time == total and wcs_count == total else "warn"),
                ("Exposure / Cadence", f"{exposure} · {cadence}", None),
                ("Frame / Flags", f"{frames[0].shape[1]} x {frames[0].shape[0]} · {flags}", "ok" if flags == 0 else "warn"),
            ]
            + camera_metadata_rows(frames[0].camera)
        )
        blocking = any("invalid_time" in frame.quality_flags or "shape_mismatch" in frame.quality_flags for frame in frames)
        self.lbl_sequence_state.setText("● Review" if blocking else "● Sequence OK")
        self.lbl_sequence_state.setObjectName("error" if blocking else "ok")

    def preview_ready(self, index, display, original_shape, thumbnail, generation=None):
        if generation is not None and generation != self.preview_generation:
            return
        self.playback_cache[int(index)] = (display, tuple(original_shape))
        self.filmstrip.set_thumbnail(int(index), QPixmap.fromImage(thumbnail))
        if index == self.current_index and self.timer.isActive():
            self.viewer.set_display_image(display, original_shape)

    def start_analysis(self, *, analysis_mode="guided", match_known_objects=True):
        if len(self.fits_files) < 3:
            QMessageBox.information(self, "Sekans gerekli", "En az üç zaman damgalı FITS karesi seçin.")
            return
        if self.analysis_worker and self.analysis_worker.isRunning():
            return
        settings = self.detection_settings()
        self.analysis_mode = analysis_mode
        self.timer.stop()
        self.result_tabs.setCurrentIndex(0)
        self.tracklet_table.setRowCount(0)
        self.known_table.setRowCount(0)
        self.difference_cache.clear()
        self.residual_levels = None
        self.pipeline_progress.reset()
        self.lbl_progress_percent.setText("0%")
        self.lbl_session_state.setText("ANALYZING")
        self.lbl_session_state.setObjectName("badgeCandidate")
        self.analysis_worker = AsteroidWorker(
            self.fits_files,
            sigma=settings["sigma"],
            min_pix=settings["min_pixels"],
            min_tracklet_frames=settings["min_frames"],
            max_sources_per_frame=settings["max_sources"],
            edge_margin=settings.get("edge_margin", 10),
            expected_fwhm_px=settings.get("expected_fwhm_px"),
            min_motion_px_per_frame=settings.get("min_motion_px_per_frame", 0.7),
            max_step_px=settings.get("max_step_px", 35.0),
            min_median_snr=settings.get("min_median_snr", 4.5),
            max_fit_rms_px=settings.get("max_fit_rms_px", 1.8),
            strong_fit_rms_px=settings.get("strong_fit_rms_px", 0.9),
            match_tolerance_px=settings.get("match_tolerance_px", 2.8),
            min_track_occupancy=settings.get("min_track_occupancy", 0.5),
            max_missing_gap_frames=settings.get("max_missing_gap_frames", 3),
            max_artifact_fraction=settings.get("max_artifact_fraction", 0.5),
            persistence_fraction=settings.get("persistence_fraction", 0.12),
            match_known_objects=match_known_objects,
            master_bias_path=self.master_paths["bias"],
            master_dark_path=self.master_paths["dark"],
            master_flat_path=self.master_paths["flat"],
            auto_plate_solve=True,
        )
        self.analysis_worker.progress.connect(self.update_progress)
        self.analysis_worker.file_progress.connect(self.update_file_progress)
        self.analysis_worker.stage_changed.connect(self.stage_changed)
        self.analysis_worker.frame_done.connect(self.analysis_frame_done)
        self.analysis_worker.sequence_ready.connect(self.analysis_sequence_ready)
        self.analysis_worker.known_objects_done.connect(self.populate_known_table)
        self.analysis_worker.tracklets_done.connect(self.populate_tracklet_table)
        self.analysis_worker.result_ready.connect(self.analysis_result_ready)
        self.analysis_worker.log.connect(self.log)
        self.analysis_worker.finished_scan.connect(self.analysis_finished)
        self.analysis_worker.start()
        self.log(f"INFO|Workflow|Analysis mode: {analysis_mode}")
        self.update_action_state()

    def start_known_recovery(self):
        self.start_analysis(analysis_mode="known-recovery", match_known_objects=True)

    def start_unknown_discovery(self):
        self.start_analysis(analysis_mode="unknown-discovery", match_known_objects=False)

    def stop_analysis(self):
        if self.analysis_worker and self.analysis_worker.isRunning():
            self.analysis_worker.stop()
            self.lbl_status.setText("Durdurma isteği gönderildi")
        if self.optimization_worker and self.optimization_worker.isRunning():
            self.optimization_worker.stop()
            self.lbl_status.setText("Optimizasyon durduruluyor")

    def stage_changed(self, key, label):
        mapping = {
            "metadata": "source",
            "solve": "solve",
            "reference": "reference",
            "detect": "detect",
            "link": "link",
            "known": "known",
            "complete": "review",
        }
        active = mapping.get(key)
        ordered = [item[0] for item in WorkflowRail.STEPS]
        if active in ordered:
            active_index = ordered.index(active)
            for index, step in enumerate(ordered):
                self.workflow.set_state(step, "done" if index < active_index else "active" if index == active_index else "pending")
        if active:
            self.pipeline_progress.set_active_stage(active)
        if key == "complete":
            self.pipeline_progress.complete()
            self.lbl_progress_percent.setText("100%")
        elif key == "error":
            self.pipeline_progress.set_stage_state(self.pipeline_progress.active_key, "error")
        self.lbl_status.setText(label)
        self.lbl_status.setToolTip(label)

    def analysis_frame_done(self, index, payload):
        self.frame_candidates[index] = payload.get("detections", [])
        preview = payload.get("residual_preview")
        if preview is not None:
            display, _ = stretch_image(preview, mode="Auto STF", sequence_levels=None)
            original_shape = self.sequence_frames[index].shape if index < len(self.sequence_frames) else preview.shape
            self.difference_cache[index] = (np.asarray(display * 255.0, dtype=np.uint8), original_shape)

    def analysis_sequence_ready(self, frames):
        self.sequence_frames = frames
        self.update_validation(frames)

    def analysis_result_ready(self, result):
        self.sequence_result = result
        self.sequence_frames = result.frames
        self.frame_candidates = {
            int(frame_index): list(detections)
            for frame_index, detections in result.detections_by_frame.items()
        }
        self.populate_tracklet_table(result.tracklets)
        self.populate_known_table(result.known_objects)
        self.result_tabs.setCurrentIndex(1 if self.analysis_mode == "known-recovery" else 0)
        self.update_fov()
        self.lbl_session_state.setText("REVIEW READY")
        self.lbl_session_state.setObjectName("badgeKnown" if result.tracklets else "badgeNeutral")
        self.workflow.set_state("review", "active")
        self._start_preview_cache(sequence_result=result)
        self.refresh_current_frame()

    def analysis_finished(self):
        self.update_action_state()
        if self.progress.value() < 100 and not self.sequence_result:
            self.lbl_session_state.setText("STOPPED")
            self.lbl_session_state.setObjectName("badgeRejected")

    def update_progress(self, value, message):
        self.progress.setValue(int(value))
        text = str(message).split("|", 2)[-1]
        self.lbl_status.setText(text)
        self.lbl_status.setToolTip(text)
        self.lbl_progress_percent.setText(f"{int(value)}%")

    def update_file_progress(self, current, total, name, stage):
        file_value = int(current / max(total, 1) * 100)
        self.pipeline_progress.set_file_value(file_value)
        text = f"{stage}: {current}/{total} · {name}"
        self.lbl_file_status.setText(text)
        self.lbl_file_status.setToolTip(text)

    def populate_tracklet_table(self, tracklets):
        self.tracklet_table.setRowCount(len(tracklets))
        for row, tracklet in enumerate(tracklets):
            summary = tracklet.to_summary()
            motion = summary["motion_arcsec_per_min"]
            if summary["known_match"]:
                class_label = "Known object"
            elif summary["classification"] == "review_candidate":
                class_label = "Borderline review"
            else:
                class_label = "Potential discovery"
            values = [
                summary["id"],
                class_label,
                f"{summary['frames']} / {len(self.fits_files)}",
                f"{summary['arc_minutes']:.2f} min",
                f"{motion:.2f}″/min" if motion is not None else f"{summary['motion_px_per_frame']:.2f} px/frame",
                f"{summary['position_angle_deg']:.1f}°",
                f"{summary['fit_rms_px']:.2f} px",
                f"{summary['median_snr']:.1f}",
                f"{summary['confidence']:.2f}",
                summary["review_status"],
            ]
            for column, value in enumerate(values):
                self.tracklet_table.setItem(row, column, QTableWidgetItem(value))
        potential = sum(item.classification == "unknown_candidate" and not item.known_match for item in tracklets)
        borderline = sum(item.classification == "review_candidate" and not item.known_match for item in tracklets)
        self.result_tabs.setTabText(0, f"Discover · {potential} potential / {borderline} borderline")

    def populate_known_table(self, predictions):
        self.known_table.setRowCount(len(predictions))
        for row, item in enumerate(predictions):
            values = [
                item.name,
                item.object_type,
                f"{item.magnitude:.2f}" if item.magnitude is not None else "--",
                f"{item.local_snr:.1f}" if item.local_snr is not None else "--",
                f"{item.local_offset_px:.1f} px" if item.local_offset_px is not None else "--",
                {
                    "high_confidence_match": "High-confidence match",
                    "plausible_visible_match": "Plausible visible match",
                    "predicted_visual_confirmation_weak": "Visual confirmation weak",
                    "predicted_in_field": "Predicted in field",
                }.get(item.status, "Predicted in field"),
                f"{item.confidence:.2f}",
            ]
            for column, value in enumerate(values):
                self.known_table.setItem(row, column, QTableWidgetItem(value))
        recovered = sum(item.visible for item in predictions)
        self.result_tabs.setTabText(1, f"Recovery · {recovered} recovered / {len(predictions) - recovered} missed")

    def tracklet_selected(self, row, _column, _previous_row, _previous_column):
        if not self.sequence_result or row < 0 or row >= len(self.sequence_result.tracklets):
            return
        self.current_tracklet = self.sequence_result.tracklets[row]
        self.current_known = None
        self.show_tracklet_evidence(self.current_tracklet)
        representative = self.current_tracklet.detections[len(self.current_tracklet.detections) // 2]
        self.change_frame(representative.frame_index)

    def known_selected(self, row, _column, _previous_row, _previous_column):
        if not self.sequence_result or row < 0 or row >= len(self.sequence_result.known_objects):
            return
        item = self.sequence_result.known_objects[row]
        self.current_known = item
        self.current_tracklet = None
        self.lbl_evidence_id.setText(item.name)
        self.lbl_verdict.setText("BİLİNEN CİSİM · GÖRÜNTÜDE DOĞRULANDI" if item.visible else "BİLİNEN CİSİM · TAHMİN, ZAYIF GÖRÜNTÜ KANITI")
        self.lbl_verdict.setObjectName("badgeKnown" if item.visible else "badgeNeutral")
        self._repolish(self.lbl_verdict)
        self.lbl_skybot.setText("Catalog prediction")
        self.lbl_artifacts.setText("None")
        self.lbl_review_state.setText("Not required")
        if item.x is not None and item.y is not None:
            self.viewer.center_on(item.x, item.y, 90)
        self.refresh_overlays()
        self.update_action_state()

    def show_tracklet_evidence(self, tracklet):
        summary = tracklet.to_summary()
        self.lbl_evidence_id.setText(f"Tracklet {summary['id']}")
        if summary["known_match"]:
            verdict, style = f"BİLİNEN CİSİM · {summary['known_match']['name']}", "badgeKnown"
        elif summary["classification"] == "review_candidate":
            verdict, style = "SINIRDA ADAY · EK BLINK / SYNTHETIC TRACK GEREKİR", "badgeNeutral"
        elif summary["artifact_flags"]:
            verdict, style = "BİLİNMEYEN ADAY · ARTEFAKT KONTROLÜ GEREKİR", "badgeCandidate"
        else:
            verdict, style = "BİLİNMEYEN HAREKETLİ CİSİM ADAYI", "badgeCandidate"
        self.lbl_verdict.setText(verdict)
        self.lbl_verdict.setObjectName(style)
        self._repolish(self.lbl_verdict)
        values = {
            "frames": f"{summary['frames']} / {len(self.fits_files)}",
            "arc": f"{summary['arc_minutes']:.2f} min",
            "motion": f"{summary['motion_arcsec_per_min']:.2f}″/min" if summary["motion_arcsec_per_min"] is not None else f"{summary['motion_px_per_frame']:.2f} px/frame",
            "pa": f"{summary['position_angle_deg']:.1f}°",
            "rms": f"{summary['fit_rms_px']:.2f} px",
            "snr": f"{summary['median_snr']:.1f}",
            "confidence": f"{summary['confidence']:.2f}",
        }
        for key, value in values.items():
            self.evidence_values[key].setText(value)
        self.lbl_skybot.setText(summary["known_match"]["name"] if summary["known_match"] else "No catalog match")
        self.lbl_artifacts.setText(", ".join(summary["artifact_flags"]) if summary["artifact_flags"] else "None")
        self.lbl_review_state.setText(summary["review_status"])
        self.notes.blockSignals(True)
        self.notes.setPlainText(tracklet.reviewer_notes)
        self.notes.blockSignals(False)
        self.refresh_overlays()
        self.update_action_state()

    def change_frame(self, index: int):
        if not self.fits_files:
            return
        self.current_index = max(0, min(int(index), len(self.fits_files) - 1))
        self.filmstrip.select(self.current_index)
        self.lbl_frame_counter.setText(f"{self.current_index + 1} / {len(self.fits_files)}")
        self.refresh_current_frame()

    def refresh_current_frame(self):
        if not self.fits_files:
            return
        mode = self.combo_mode.currentText()
        try:
            if self.timer.isActive():
                cache = self.difference_cache if mode == "Difference" else self.playback_cache
                if self.current_index in cache:
                    display, shape = cache[self.current_index]
                    self.viewer.set_display_image(display, shape)
                    self.refresh_overlays()
                    return
            if mode == "Original" and self.sequence_result is None and self.current_index in self.playback_cache:
                display, shape = self.playback_cache[self.current_index]
                self.viewer.set_display_image(display, shape)
                self.refresh_overlays()
                return
            if mode == "Synthetic Track" and self.synthetic_result is not None:
                self.viewer.set_image(self.synthetic_result.image, stretch_mode="Auto STF", keep_view=False)
            elif self.sequence_result is not None:
                data = materialize_sequence_frame(
                    self.sequence_result,
                    self.current_index,
                    residual=mode == "Difference",
                )
                self.viewer.set_image(
                    data,
                    stretch_mode=self.combo_stretch.currentText(),
                    sequence_levels=None,
                    manual_levels=self.manual_levels,
                    invert=self.chk_invert.isChecked(),
                )
            else:
                data = read_fits_image(self.fits_files[self.current_index])
                self.viewer.set_image(
                    data,
                    stretch_mode=self.combo_stretch.currentText(),
                    sequence_levels=None,
                    manual_levels=self.manual_levels,
                    invert=self.chk_invert.isChecked(),
                )
            self.refresh_overlays()
        except Exception as exc:
            self.log(f"ERROR|Viewer|{Path(self.fits_files[self.current_index]).name}: {exc}")

    def refresh_overlays(self):
        known = []
        candidates = self.frame_candidates.get(self.current_index, [])
        if self.sequence_result and self.current_index == self.sequence_result.reference_index:
            known = self.sequence_result.known_objects
        self.viewer.set_overlays(
            known=known,
            candidates=candidates,
            selected_tracklet=self.current_tracklet,
            show_known=self.chk_known.isChecked(),
            show_candidates=self.chk_candidates.isChecked(),
            show_track=self.chk_track.isChecked(),
        )
        self.viewer.set_grid_visible(self.chk_grid.isChecked())

    def toggle_play(self):
        if not self.fits_files:
            return
        if self.combo_mode.currentText() == "Synthetic Track":
            self.combo_mode.setCurrentText("Original")
        if self.timer.isActive():
            self.timer.stop()
            self.btn_play.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
            self.refresh_current_frame()
        else:
            mode = self.combo_mode.currentText()
            cache = self.difference_cache if mode == "Difference" else self.playback_cache
            if len(cache) != len(self.fits_files):
                self.lbl_status.setText("Blink için hizalı Auto STF önbelleği bekleniyor")
                self.log(f"WARN|Blink|Önbellek hazır değil: {len(cache)}/{len(self.fits_files)} kare")
                return
            self.timer.start(self.current_blink_interval())
            self.btn_play.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPause))

    def current_blink_interval(self):
        try:
            fps = int(self.combo_blink_speed.currentText().split()[0])
        except Exception:
            fps = 6
        return max(40, int(1000 / max(fps, 1)))

    def update_blink_speed(self):
        if self.timer.isActive():
            self.timer.start(self.current_blink_interval())

    def next_frame(self):
        if self.fits_files:
            self.change_frame((self.current_index + 1) % len(self.fits_files))

    def previous_frame(self):
        if self.fits_files:
            self.change_frame((self.current_index - 1) % len(self.fits_files))

    def stretch_changed(self, mode):
        if mode == "Manual STF":
            dialog = ManualStfDialog(self)
            if dialog.exec() == QDialog.DialogCode.Accepted:
                self.manual_levels = dialog.values()
            else:
                self.combo_stretch.setCurrentText("Auto STF")
                return
        self.refresh_current_frame()

    def choose_master(self, kind: str):
        path, _ = QFileDialog.getOpenFileName(
            self,
            f"Master {kind} FITS seç",
            self.current_folder,
            "FITS (*.fit *.fits *.fts *.fits.gz)",
        )
        if not path:
            return
        try:
            data = read_fits_image(path)
            if self.sequence_frames and data.shape != self.sequence_frames[0].shape:
                raise ValueError(f"Boyut uyuşmuyor: {data.shape} != {self.sequence_frames[0].shape}")
            self.master_paths[kind] = path
            self.master_labels[kind].setText(Path(path).name)
            self.master_labels[kind].setToolTip(path)
            self.log(f"INFO|Calibration|Master {kind} seçildi: {Path(path).name}")
            self.workflow.set_state("calibration", "done")
            self.workflow.set_state("solve", "active")
        except Exception as exc:
            QMessageBox.warning(self, "Kalibrasyon", str(exc))

    def edit_detection_settings(self):
        dialog = DetectionSettingsDialog(self.custom_detection, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.custom_detection = dialog.values()
            self.combo_profile.setCurrentText("Custom")

    def detection_settings(self):
        profile = self.combo_profile.currentText()
        if profile == "Documented Workflow":
            return dict(self.custom_detection)
        if profile == "Conservative":
            return {
                **self.custom_detection,
                "sigma": 5.5,
                "min_pixels": 5,
                "min_frames": 5 if len(self.fits_files) > 15 else 4,
                "max_sources": 20,
                "min_median_snr": 5.6,
                "max_fit_rms_px": 1.2,
                "strong_fit_rms_px": 0.75,
                "min_track_occupancy": 0.65,
            }
        if profile == "Balanced":
            return {
                **self.custom_detection,
                "sigma": 5.0,
                "min_pixels": 5,
                "min_frames": 4 if len(self.fits_files) > 8 else 3,
                "max_sources": 28,
                "min_median_snr": 5.0,
                "max_fit_rms_px": 1.6,
                "strong_fit_rms_px": 0.95,
                "min_track_occupancy": 0.55,
            }
        return dict(self.custom_detection)

    def start_detection_optimization(self):
        if len(self.fits_files) < 3:
            QMessageBox.information(self, "Otomatik Ayar", "Önce en az üç zaman sıralı FITS karesi yükleyin.")
            return
        if self.optimization_worker and self.optimization_worker.isRunning():
            return
        self.optimization_worker = DetectionOptimizationWorker(self.fits_files, self.master_paths, self)
        self.optimization_worker.progress.connect(self.optimization_progress)
        self.optimization_worker.result_ready.connect(self.optimization_ready)
        self.optimization_worker.failed.connect(self.optimization_failed)
        self.optimization_worker.finished.connect(self.update_action_state)
        self.optimization_worker.start()
        self.pipeline_progress.reset()
        self.pipeline_progress.set_active_stage("detect")
        self.lbl_session_state.setText("OPTIMIZING")
        self.lbl_session_state.setObjectName("badgeCandidate")
        self._repolish(self.lbl_session_state)
        self.log("INFO|Optimize|Görüntüye özel tespit profili ölçülüyor")
        self.update_action_state()

    def optimization_progress(self, value: int, message: str):
        self.pipeline_progress.setValue(value)
        self.lbl_progress_percent.setText(f"{value}%")
        self.lbl_status.setText(message)
        self.lbl_status.setToolTip(message)

    def optimization_ready(self, result):
        self.custom_detection = dict(result.settings)
        self.combo_profile.setCurrentText("Auto Optimized")
        summary = result.summary()
        scale = result.metrics["pixel_scale_arcsec"]
        scale_text = f"{float(scale):.2f}″/px" if isinstance(scale, (int, float)) else "WCS yok"
        self.lbl_optimization.setText(
            f"σ{result.settings['sigma']:.1f} · FWHM {result.metrics['fwhm_px']:.2f}px · {scale_text}"
        )
        details = " · ".join(f"{key}: {value}" for key, value in result.metrics.items())
        self.lbl_optimization.setToolTip(details)
        self.pipeline_progress.complete()
        self.lbl_progress_percent.setText("100%")
        self.lbl_status.setText("Otomatik profil uygulandı")
        self.lbl_session_state.setText("PROFILE READY")
        self.lbl_session_state.setObjectName("badgeKnown")
        self._repolish(self.lbl_session_state)
        self.log(f"INFO|Optimize|{summary}")
        for warning in result.warnings:
            self.log(f"WARN|Optimize|{warning}")
        QMessageBox.information(
            self,
            "Otomatik Tespit Profili",
            f"Profil uygulandı.\n\n{summary}\n"
            f"Residual P90: {result.metrics['residuals_p90']} aday/kare\n"
            f"Hizalama RMS: {result.metrics['registration_rms_px']} px\n"
            f"Kare başına üst sınır: {result.settings['max_sources']}",
        )
        self.update_action_state()

    def optimization_failed(self, message: str):
        self.lbl_status.setText("Otomatik profil üretilemedi")
        self.lbl_optimization.setText("Ölçüm başarısız")
        self.lbl_optimization.setToolTip(message)
        self.pipeline_progress.set_stage_state("detect", "error")
        self.lbl_session_state.setText("OPTIMIZE ERROR")
        self.lbl_session_state.setObjectName("badgeRejected")
        self._repolish(self.lbl_session_state)
        self.log(f"ERROR|Optimize|{message}")
        QMessageBox.warning(self, "Otomatik Ayar", message)
        self.update_action_state()

    def solve_selected_frame(self):
        if not self.fits_files or (self.plate_worker and self.plate_worker.isRunning()):
            return
        self.plate_worker = PlateSolveWorker(self.fits_files[self.current_index], self)
        self.plate_worker.result_ready.connect(self.plate_solve_finished)
        self.lbl_status.setText("ASTAP plate solve çalışıyor")
        self.plate_worker.start()
        self.update_action_state()

    def plate_solve_finished(self, result):
        if result.success:
            self.log("INFO|Astrometry|Seçili kare geçici kopya üzerinde çözüldü; orijinal FITS değiştirilmedi")
            self.lbl_status.setText("Plate solve başarılı")
        else:
            self.log(f"ERROR|Astrometry|Plate solve başarısız: {result.message}")
            QMessageBox.warning(self, "Plate Solve", result.message or "ASTAP çözüm üretemedi")
        self.update_action_state()

    def run_synthetic_track(self):
        if not self.sequence_result or self.current_tracklet is None:
            QMessageBox.information(self, "Tracklet seçin", "Synthetic Track için önce bir tracklet seçin.")
            return
        if self.synthetic_worker and self.synthetic_worker.isRunning():
            return
        self.lbl_synthetic.setText("Kareler aday hızına göre merkezleniyor...")
        self.synthetic_worker = SyntheticTrackWorker(self.sequence_result, self.current_tracklet, parent=self)
        self.synthetic_worker.progress.connect(
            lambda current, total: self.pipeline_progress.set_file_value(int(current / max(total, 1) * 100))
        )
        self.synthetic_worker.result_ready.connect(self.synthetic_track_ready)
        self.synthetic_worker.failed.connect(self.synthetic_track_failed)
        self.synthetic_worker.start()
        self.update_action_state()

    def synthetic_track_ready(self, result):
        self.synthetic_result = result
        verdict = "güçlü merkez sinyali" if result.snr >= 6.0 and result.peak_offset_px <= 3.0 else "zayıf veya merkez dışı sinyal"
        self.lbl_synthetic.setText(
            f"{result.used_frames} kare · SNR {result.snr:.1f} · offset {result.peak_offset_px:.2f} px · {verdict}"
        )
        self.combo_mode.setCurrentText("Synthetic Track")
        self.log(f"INFO|Synthetic Track|SNR {result.snr:.2f}, offset {result.peak_offset_px:.2f} px")
        self.update_action_state()

    def synthetic_track_failed(self, message):
        self.lbl_synthetic.setText(f"Synthetic Track başarısız: {message}")
        self.log(f"ERROR|Synthetic Track|{message}")
        self.update_action_state()

    def set_review_status(self, status: str):
        if self.current_tracklet is None:
            return
        previous = self.current_tracklet.review_status
        self.review_history.append((self.current_tracklet, previous))
        self.current_tracklet.review_status = status
        self.lbl_review_state.setText(status)
        self.populate_tracklet_table(self.sequence_result.tracklets)
        row = self.sequence_result.tracklets.index(self.current_tracklet)
        self.tracklet_table.selectRow(row)
        self.log(f"INFO|Review|{self.current_tracklet.tracklet_id}: {status}")
        self.update_action_state()

    def undo_review(self):
        if not self.review_history:
            return
        tracklet, status = self.review_history.pop()
        tracklet.review_status = status
        if tracklet is self.current_tracklet:
            self.lbl_review_state.setText(status)
        self.populate_tracklet_table(self.sequence_result.tracklets)
        self.log(f"INFO|Review|{tracklet.tracklet_id}: karar geri alındı")
        self.update_action_state()

    def notes_changed(self):
        if self.current_tracklet is not None:
            self.current_tracklet.reviewer_notes = self.notes.toPlainText()

    def center_selected_tracklet(self):
        if self.current_tracklet is None:
            return
        detection = min(
            self.current_tracklet.detections,
            key=lambda item: abs(item.frame_index - self.current_index),
        )
        self.viewer.center_on(detection.x, detection.y, 75)

    def export_csv(self):
        if not self.sequence_result or not self.sequence_result.tracklets:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Tracklet CSV", "voyager-alpha-tracklets.csv", "CSV (*.csv)")
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8-sig") as handle:
            fieldnames = ["tracklet_id", "classification", "review", "frame", "x", "y", "ra", "dec", "snr", "fit_rms_px", "motion_arcsec_min", "known_match", "notes"]
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for tracklet in self.sequence_result.tracklets:
                for detection in tracklet.detections:
                    writer.writerow(
                        {
                            "tracklet_id": tracklet.tracklet_id,
                            "classification": tracklet.classification,
                            "review": tracklet.review_status,
                            "frame": detection.frame_index + 1,
                            "x": f"{detection.x:.3f}",
                            "y": f"{detection.y:.3f}",
                            "ra": "" if detection.ra is None else f"{detection.ra:.8f}",
                            "dec": "" if detection.dec is None else f"{detection.dec:.8f}",
                            "snr": f"{detection.snr:.3f}",
                            "fit_rms_px": f"{tracklet.fit_rms_px:.3f}",
                            "motion_arcsec_min": "" if tracklet.motion_arcsec_per_min is None else f"{tracklet.motion_arcsec_per_min:.4f}",
                            "known_match": tracklet.known_match["name"] if tracklet.known_match else "",
                            "notes": tracklet.reviewer_notes,
                        }
                    )
        self.log(f"INFO|Export|CSV kaydedildi: {path}")

    def export_report(self):
        if not self.sequence_result:
            return
        path, _ = QFileDialog.getSaveFileName(self, "HTML Report", "voyager-alpha-report.html", "HTML (*.html)")
        if path:
            Path(path).write_text(render_tracklet_html_report(self.sequence_result.frames, self.sequence_result.tracklets), encoding="utf-8")
            self.log(f"INFO|Export|HTML rapor kaydedildi: {path}")

    def export_ades(self):
        if not self.sequence_result:
            return
        accepted = [item for item in self.sequence_result.tracklets if item.review_status == "accepted"]
        if not accepted:
            QMessageBox.information(self, "ADES", "Önce en az bir tracklet'i Accept olarak işaretleyin.")
            return
        code, ok = QInputDialog.getText(
            self,
            "MPC Observatory Code",
            "MPC tarafından atanmış 3 karakterli gözlemevi kodu. Resmi kodunuz yoksa boş bırakın:",
        )
        if not ok:
            return
        code = code.strip().upper()
        if code and len(code) != 3:
            QMessageBox.warning(self, "MPC Code", "Observatory Code tam 3 karakter olmalıdır.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "ADES Draft", "voyager-alpha-ades-draft.psv", "PSV (*.psv)")
        if path:
            Path(path).write_text(render_ades_psv_draft(accepted, self.sequence_result.frames, code), encoding="utf-8")
            self.log(f"INFO|Export|ADES taslağı kaydedildi: {path}")

    def workflow_step_clicked(self, key: str):
        messages = {
            "source": "FITS klasörü seçin; DATE-OBS, pozlama ve boyut tutarlılığı doğrulanır.",
            "calibration": "Bias, dark ve flat isteğe bağlıdır; dark poz süresine göre ölçeklenir.",
            "solve": "WCS yoksa referans kare ASTAP ile geçici kopya üzerinde çözülür.",
            "reference": "Yıldız eşleşmeleriyle affine hizalama ve statik gökyüzü modeli kurulur.",
            "known": "SkyBoT tahminleri yerel piksel sinyaliyle doğrulanır.",
            "detect": "Residual kaynaklar gerçek zaman ekseninde hareket zincirlerine bağlanır.",
            "review": "Bilinmeyen adaylar blink ve Synthetic Track ile insan incelemesine sunulur.",
        }
        self.lbl_status.setText(messages.get(key, ""))

    def update_fov(self):
        if not self.sequence_result:
            return
        try:
            from astropy.wcs import WCS
            from astropy.wcs.utils import proj_plane_pixel_scales

            wcs = WCS(self.sequence_result.reference_header).celestial
            scales = np.abs(proj_plane_pixel_scales(wcs))
            height, width = self.sequence_result.frames[0].shape
            self.lbl_fov.setText(f"FOV: {width * scales[0]:.2f}° × {height * scales[1]:.2f}°")
        except Exception:
            self.lbl_fov.setText("FOV: pixel-space")

    def viewer_auto_fit(self):
        self.viewer.auto_range()

    def update_action_state(self):
        running = bool(self.analysis_worker and self.analysis_worker.isRunning())
        optimizing = bool(self.optimization_worker and self.optimization_worker.isRunning())
        synthetic_running = bool(self.synthetic_worker and self.synthetic_worker.isRunning())
        plate_running = bool(self.plate_worker and self.plate_worker.isRunning())
        has_files = bool(self.fits_files)
        selected = self.current_tracklet is not None
        self.btn_open.setEnabled(not running and not optimizing)
        self.btn_start.setEnabled(has_files and not running and not optimizing)
        self.btn_stop.setEnabled(running or optimizing)
        self.btn_rescan.setEnabled(bool(self.current_folder) and not running and not optimizing)
        self.btn_optimize_detection.setEnabled(len(self.fits_files) >= 3 and not running and not optimizing)
        self.btn_solve.setEnabled(has_files and not running and not optimizing and not plate_running)
        mode = self.combo_mode.currentText() if has_files else "Original"
        active_cache = self.difference_cache if mode == "Difference" else self.playback_cache
        self.btn_play.setEnabled(has_files and len(active_cache) == len(self.fits_files))
        self.btn_prev.setEnabled(has_files)
        self.btn_next.setEnabled(has_files)
        for button in (self.btn_accept, self.btn_reject, self.btn_follow, self.btn_center):
            button.setEnabled(selected)
        self.btn_synthetic.setEnabled(selected and self.sequence_result is not None and not synthetic_running)
        self.btn_undo_review.setEnabled(bool(self.review_history))
        self.btn_export_csv.setEnabled(bool(self.sequence_result and self.sequence_result.tracklets))
        self.btn_export_report.setEnabled(self.sequence_result is not None)
        self.btn_export_ades.setEnabled(bool(self.sequence_result and any(item.review_status == "accepted" for item in self.sequence_result.tracklets)))

    def log(self, payload: str):
        self.log_table.append_message(payload)

    def _repolish(self, widget):
        widget.style().unpolish(widget)
        widget.style().polish(widget)

    def shutdown_workers(self):
        self.timer.stop()
        for worker in (self.preview_worker, self.analysis_worker, self.synthetic_worker, self.plate_worker, self.optimization_worker):
            if worker and worker.isRunning():
                if hasattr(worker, "stop"):
                    worker.stop()
                worker.wait(3000)


class DetectionSettingsDialog(QDialog):
    def __init__(self, values, parent=None):
        super().__init__(parent)
        self._base_values = dict(values)
        self.setWindowTitle("Discover Pipeline Settings")
        self.setMinimumWidth(430)
        form = QFormLayout(self)
        self.sigma = QLineEdit(str(values["sigma"]))
        self.sigma.setValidator(QDoubleValidator(3.0, 12.0, 2, self))
        self.min_pixels = QLineEdit(str(values["min_pixels"]))
        self.min_pixels.setValidator(QIntValidator(5, 100, self))
        self.min_frames = QLineEdit(str(values["min_frames"]))
        self.min_frames.setValidator(QIntValidator(3, 20, self))
        self.max_sources = QLineEdit(str(values["max_sources"]))
        self.max_sources.setValidator(QIntValidator(8, 100, self))
        self.min_motion = QLineEdit(str(values.get("min_motion_px_per_frame", 1.5)))
        self.min_motion.setValidator(QDoubleValidator(0.05, 20.0, 3, self))
        self.max_step = QLineEdit(str(values.get("max_step_px", 35.0)))
        self.max_step.setValidator(QDoubleValidator(1.0, 200.0, 2, self))
        self.max_fit_rms = QLineEdit(str(values.get("max_fit_rms_px", 1.8)))
        self.max_fit_rms.setValidator(QDoubleValidator(0.3, 5.0, 2, self))
        self.track_occupancy = QLineEdit(str(values.get("min_track_occupancy", 0.58)))
        self.track_occupancy.setValidator(QDoubleValidator(0.3, 1.0, 2, self))
        form.addRow("Residual threshold (sigma)", self.sigma)
        form.addRow("Minimum source area (px)", self.min_pixels)
        form.addRow("Minimum linked frames", self.min_frames)
        form.addRow("Maximum residuals / frame", self.max_sources)
        form.addRow("Minimum motion (px/frame)", self.min_motion)
        form.addRow("Maximum step (px/frame)", self.max_step)
        form.addRow("Maximum linear-fit RMS (px)", self.max_fit_rms)
        form.addRow("Minimum track occupancy (0-1)", self.track_occupancy)
        note = QLabel("Dokümante başlangıç: hybrid detector, 5 sigma, FWHM 3 px, en fazla 24 residual, en az 3 kare, 1.5 px seed hareketi ve 2.8 px yeniden eşleme.")
        note.setWordWrap(True)
        note.setObjectName("muted")
        form.addRow(note)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    def values(self):
        values = dict(self._base_values)
        values.update({
            "sigma": float(self.sigma.text()),
            "min_pixels": int(self.min_pixels.text()),
            "min_frames": int(self.min_frames.text()),
            "max_sources": int(self.max_sources.text()),
            "min_motion_px_per_frame": float(self.min_motion.text()),
            "max_step_px": float(self.max_step.text()),
            "max_fit_rms_px": float(self.max_fit_rms.text()),
            "strong_fit_rms_px": min(float(self.max_fit_rms.text()) * 0.62, float(self.max_fit_rms.text())),
            "min_track_occupancy": float(self.track_occupancy.text()),
        })
        return values


class ManualStfDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Manual STF")
        form = QFormLayout(self)
        self.black = QLineEdit("0.1")
        self.midtone = QLineEdit("0.35")
        self.white = QLineEdit("99.9")
        for editor in (self.black, self.midtone, self.white):
            editor.setValidator(QDoubleValidator(0.0, 100.0, 4, self))
        form.addRow("Black percentile", self.black)
        form.addRow("Midtone (0-1)", self.midtone)
        form.addRow("White percentile", self.white)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Apply | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    def values(self):
        return float(self.black.text()), float(self.midtone.text()), float(self.white.text())


def find_fits_files(folder: str) -> list[str]:
    root = Path(folder)
    extensions = (".fit", ".fits", ".fts", ".fits.gz")
    direct = [path for path in root.iterdir() if path.is_file() and path.name.lower().endswith(extensions)]
    files = direct if direct else [path for path in root.rglob("*") if path.is_file() and path.name.lower().endswith(extensions)]
    science = []
    for path in files:
        lowered = path.name.lower()
        if any(token in lowered for token in ("masterdark", "master_dark", "masterflat", "master_flat", "masterbias", "master_bias")):
            continue
        science.append(str(path))
    return sorted(science, key=lambda value: Path(value).name.lower())
