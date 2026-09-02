from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from PyQt6.QtCore import QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QPainter, QPen, QPixmap
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)


class PipelineProgressStrip(QWidget):
    """Compact, segmented pipeline status with a secondary per-file track."""

    COLORS = (
        "#3b82c4",
        "#2f9db5",
        "#25aa9a",
        "#45ad68",
        "#b89d38",
        "#c07a36",
        "#9b6ac6",
        "#42aab1",
    )

    def __init__(self, stages: list[tuple[str, str]], parent=None):
        super().__init__(parent)
        self.stages = list(stages)
        self.states = {key: "pending" for key, _label in self.stages}
        self.active_key = self.stages[0][0] if self.stages else ""
        self._value = 0
        self._file_value = 0
        self.setMinimumWidth(520)
        self.setFixedHeight(34)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setToolTip(
            "Renkli sektörler analiz aşamalarını gösterir. Alt ince çizgi, aktif aşamadaki dosya ilerlemesidir."
        )
        if self.active_key:
            self.states[self.active_key] = "active"

    def reset(self):
        self._value = 0
        self._file_value = 0
        self.active_key = self.stages[0][0] if self.stages else ""
        self.states = {key: "pending" for key, _label in self.stages}
        if self.active_key:
            self.states[self.active_key] = "active"
        self.update()

    def setValue(self, value: int):
        self._value = max(0, min(100, int(value)))
        self.update()

    def value(self) -> int:
        return self._value

    def set_file_value(self, value: int):
        self._file_value = max(0, min(100, int(value)))
        self.update()

    def set_active_stage(self, key: str):
        keys = [stage_key for stage_key, _label in self.stages]
        if key not in keys:
            return
        active_index = keys.index(key)
        self.active_key = key
        for index, stage_key in enumerate(keys):
            self.states[stage_key] = "done" if index < active_index else "active" if index == active_index else "pending"
        self._file_value = 0
        self.update()

    def set_stage_state(self, key: str, state: str):
        if key in self.states and state in {"pending", "active", "done", "warning", "error"}:
            self.states[key] = state
            if state == "active":
                self.active_key = key
            self.update()

    def complete(self):
        self._value = 100
        self._file_value = 100
        for key in self.states:
            self.states[key] = "done"
        self.active_key = ""
        self.update()

    def paintEvent(self, event):
        del event
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        outer = QRectF(0.5, 2.5, max(1.0, self.width() - 1.0), 28.0)
        painter.setPen(QPen(QColor("#2b3940"), 1.0))
        painter.setBrush(QColor("#10171b"))
        painter.drawRoundedRect(outer, 4.0, 4.0)
        if not self.stages:
            return

        margin = 4.0
        gap = 3.0
        count = len(self.stages)
        available = max(1.0, outer.width() - margin * 2.0 - gap * (count - 1))
        segment_width = available / count
        font = painter.font()
        font.setPointSizeF(8.2)
        painter.setFont(font)
        metrics = painter.fontMetrics()

        for index, (key, label) in enumerate(self.stages):
            x = outer.left() + margin + index * (segment_width + gap)
            segment = QRectF(x, outer.top() + 4.0, segment_width, 18.0)
            state = self.states.get(key, "pending")
            color = QColor(self.COLORS[index % len(self.COLORS)])
            if state == "error":
                color = QColor("#d6534a")
            elif state == "warning":
                color = QColor("#d2a52f")

            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor("#1b262b"))
            painter.drawRoundedRect(segment, 2.5, 2.5)
            if state == "done":
                fill_ratio = 1.0
            elif state in {"active", "warning", "error"}:
                fill_ratio = max(0.13, self._active_fill(index))
            else:
                fill_ratio = 0.0
            if fill_ratio > 0:
                fill = QRectF(segment.left(), segment.top(), max(3.0, segment.width() * fill_ratio), segment.height())
                painter.setBrush(color)
                painter.drawRoundedRect(fill, 2.5, 2.5)

            if segment.width() >= 42:
                text = metrics.elidedText(label, Qt.TextElideMode.ElideRight, max(10, int(segment.width() - 10)))
                painter.setPen(QColor("#f4f8f9") if state != "pending" else QColor("#81939b"))
                painter.drawText(segment.adjusted(5.0, 0.0, -5.0, 0.0), Qt.AlignmentFlag.AlignCenter, text)

        track = QRectF(outer.left() + margin, outer.bottom() - 3.0, outer.width() - margin * 2.0, 2.0)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor("#26343a"))
        painter.drawRoundedRect(track, 1.0, 1.0)
        if self._file_value > 0:
            file_fill = QRectF(track.left(), track.top(), track.width() * self._file_value / 100.0, track.height())
            painter.setBrush(QColor("#63d3db"))
            painter.drawRoundedRect(file_fill, 1.0, 1.0)

    def _active_fill(self, index: int) -> float:
        count = max(1, len(self.stages))
        normalized = self._value / 100.0 * count
        return max(0.0, min(1.0, normalized - index))


class SegmentedControl(QFrame):
    changed = pyqtSignal(str)

    def __init__(self, options: list[tuple[str, str]], parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)
        self.group = QButtonGroup(self)
        self.group.setExclusive(True)
        self.buttons = {}
        for index, (key, label) in enumerate(options):
            button = QPushButton(label)
            button.setCheckable(True)
            button.setObjectName("segment")
            button.clicked.connect(lambda _checked=False, value=key: self.changed.emit(value))
            self.group.addButton(button)
            self.buttons[key] = button
            layout.addWidget(button)
            if index == 0:
                button.setChecked(True)

    def value(self) -> str:
        for key, button in self.buttons.items():
            if button.isChecked():
                return key
        return next(iter(self.buttons))

    def set_value(self, key: str):
        if key in self.buttons:
            self.buttons[key].setChecked(True)
            self.changed.emit(key)


class WorkflowRail(QFrame):
    step_clicked = pyqtSignal(str)

    STEPS = [
        ("source", "1", "Kaynak ve sekans"),
        ("calibration", "2", "Kalibrasyon"),
        ("solve", "3", "Astrometri"),
        ("reference", "4", "Hizalama ve model"),
        ("known", "5", "Bilinen cisimler"),
        ("detect", "6", "Bilinmeyenleri tara"),
        ("review", "7", "Kanıt inceleme"),
    ]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("section")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(7, 7, 7, 7)
        layout.setSpacing(1)
        title = QLabel("GUIDED WORKFLOW")
        title.setObjectName("caption")
        layout.addWidget(title)
        self.buttons = {}
        for key, number, label in self.STEPS:
            button = QPushButton(f"{number}   {label}")
            button.setObjectName("workflowStep")
            button.setProperty("stepState", "pending")
            button.clicked.connect(lambda _checked=False, value=key: self.step_clicked.emit(value))
            self.buttons[key] = button
            layout.addWidget(button)

    def set_state(self, key: str, state: str):
        button = self.buttons.get(key)
        if button is None:
            return
        button.setProperty("stepState", state)
        button.style().unpolish(button)
        button.style().polish(button)

    def reset(self):
        for key in self.buttons:
            self.set_state(key, "pending")
        self.set_state("source", "active")


class KeyValueTable(QTableWidget):
    def __init__(self, parent=None):
        super().__init__(0, 2, parent)
        self.setHorizontalHeaderLabels(["Metric", "Value"])
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setShowGrid(False)
        self.verticalHeader().setDefaultSectionSize(23)
        self._auto_height = False

    def set_auto_height(self, enabled: bool = True):
        self._auto_height = bool(enabled)
        policy = Qt.ScrollBarPolicy.ScrollBarAlwaysOff if self._auto_height else Qt.ScrollBarPolicy.ScrollBarAsNeeded
        self.setVerticalScrollBarPolicy(policy)
        self.setHorizontalScrollBarPolicy(policy)
        self._sync_height()

    def _sync_height(self):
        if not self._auto_height:
            return
        row_height = self.verticalHeader().defaultSectionSize()
        content_height = self.rowCount() * row_height
        self.setFixedHeight(max(48, content_height + self.frameWidth() * 2 + 2))

    def set_rows(self, rows: list[tuple[str, str, str | None]]):
        self.setRowCount(len(rows))
        for row, (label, value, tone) in enumerate(rows):
            left = QTableWidgetItem(label)
            right = QTableWidgetItem(value)
            left.setToolTip(label)
            right.setToolTip(value)
            right.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            if tone == "ok":
                right.setForeground(QColor("#4dc56a"))
            elif tone == "warn":
                right.setForeground(QColor("#efbc3f"))
            elif tone == "error":
                right.setForeground(QColor("#ef6257"))
            elif tone == "cyan":
                right.setForeground(QColor("#2dd7df"))
            self.setItem(row, 0, left)
            self.setItem(row, 1, right)
        self._sync_height()


def camera_metadata_rows(camera) -> list[tuple[str, str, str | None]]:
    camera_name = camera.instrument or camera.detector or "Header missing"
    sensor_parts = []
    if camera.detector and camera.detector != camera_name:
        sensor_parts.append(camera.detector)
    if camera.pixel_size_x_um is not None:
        if camera.pixel_size_y_um is not None and abs(camera.pixel_size_y_um - camera.pixel_size_x_um) > 1e-3:
            sensor_parts.append(f"{camera.pixel_size_x_um:g}x{camera.pixel_size_y_um:g} um")
        else:
            sensor_parts.append(f"{camera.pixel_size_x_um:g} um")
    if camera.binning_x and camera.binning_y:
        sensor_parts.append(f"{camera.binning_x}x{camera.binning_y} bin")

    capture_parts = []
    if camera.gain is not None:
        capture_parts.append(f"G {camera.gain:g}")
    if camera.offset is not None:
        capture_parts.append(f"O {camera.offset:g}")
    if camera.sensor_temperature_c is not None:
        capture_parts.append(f"{camera.sensor_temperature_c:g} C")
    if camera.filter_name:
        capture_parts.append(camera.filter_name)

    optics_parts = []
    if camera.image_scale_arcsec_px is not None:
        optics_parts.append(f"{camera.image_scale_arcsec_px:.3f} arcsec/px")
    if camera.focal_length_mm is not None:
        optics_parts.append(f"{camera.focal_length_mm:g} mm FL")
    if camera.aperture_mm is not None:
        optics_parts.append(f"{camera.aperture_mm:g} mm D")
    if camera.readout_mode:
        optics_parts.append(camera.readout_mode)
    elif camera.bayer_pattern:
        optics_parts.append(f"Bayer {camera.bayer_pattern}")

    return [
        ("Camera", camera_name, None if camera_name != "Header missing" else "warn"),
        ("Sensor / Bin", " · ".join(sensor_parts) or "Header missing", None if sensor_parts else "warn"),
        ("Capture", " · ".join(capture_parts) or "Header missing", None if capture_parts else "warn"),
        ("Scale / Optics", " · ".join(optics_parts) or "Header missing", None if optics_parts else "warn"),
    ]


class AnalysisLog(QTableWidget):
    def __init__(self, parent=None):
        super().__init__(0, 4, parent)
        self.setHorizontalHeaderLabels(["UTC", "Level", "Source", "Message"])
        self.verticalHeader().setVisible(False)
        self.setWordWrap(False)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.setAlternatingRowColors(True)
        self.verticalHeader().setDefaultSectionSize(23)
        header = self.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        self.setColumnWidth(0, 70)
        self.setColumnWidth(1, 62)
        self.setColumnWidth(2, 95)

    def append_message(self, payload: str):
        level, source, message = parse_log_payload(payload)
        row = self.rowCount()
        self.insertRow(row)
        values = [datetime.now(UTC).strftime("%H:%M:%S"), level, source, message.replace("\n", " ")]
        color = {
            "INFO": QColor("#4dc56a"),
            "WARN": QColor("#efbc3f"),
            "ERROR": QColor("#ef6257"),
        }.get(level, QColor("#9fb1b9"))
        for column, value in enumerate(values):
            item = QTableWidgetItem(value)
            item.setToolTip(value)
            if column == 1:
                item.setForeground(color)
            self.setItem(row, column, item)
        if self.rowCount() > 1000:
            self.removeRow(0)
        self.scrollToBottom()


class ThumbnailCard(QToolButton):
    """Responsive filmstrip card with an edge-to-edge cropped preview."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._thumbnail = QPixmap()
        self._label = ""
        self.setMinimumSize(68, 94)
        self.setMaximumHeight(100)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    def set_thumbnail(self, pixmap: QPixmap | None):
        self._thumbnail = QPixmap(pixmap) if pixmap is not None else QPixmap()
        self.update()

    def set_label(self, label: str):
        self._label = str(label)
        self.update()

    def enterEvent(self, event):
        self.update()
        super().enterEvent(event)

    def leaveEvent(self, event):
        self.update()
        super().leaveEvent(event)

    def paintEvent(self, event):
        del event
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        selected = bool(self.property("selected"))
        hovered = self.underMouse()
        border = QColor("#31c8d0") if selected else QColor("#58717a") if hovered else QColor("#314047")
        background = QColor("#14292e") if selected else QColor("#10171b")
        outer = QRectF(0.5, 0.5, max(1.0, self.width() - 1.0), max(1.0, self.height() - 1.0))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(background)
        painter.drawRoundedRect(outer, 3.0, 3.0)

        label_height = 21.0
        image_rect = outer.adjusted(1.0, 1.0, -1.0, -label_height)
        if not self._thumbnail.isNull() and image_rect.width() > 1 and image_rect.height() > 1:
            scaled = self._thumbnail.scaled(
                max(1, int(image_rect.width())),
                max(1, int(image_rect.height())),
                Qt.AspectRatioMode.KeepAspectRatioByExpanding,
                Qt.TransformationMode.SmoothTransformation,
            )
            painter.save()
            painter.setClipRect(image_rect)
            x = int(image_rect.center().x() - scaled.width() / 2)
            y = int(image_rect.center().y() - scaled.height() / 2)
            painter.drawPixmap(x, y, scaled)
            painter.restore()
        else:
            painter.fillRect(image_rect, QColor("#090d0f"))

        label_rect = QRectF(outer.left() + 2.0, outer.bottom() - label_height, outer.width() - 4.0, label_height)
        painter.fillRect(label_rect, QColor("#14292e") if selected else QColor("#10171b"))
        font = painter.font()
        font.setPointSizeF(8.2)
        font.setWeight(600 if selected else 500)
        painter.setFont(font)
        painter.setPen(QColor("#f2ffff") if selected else QColor("#c2cdd1"))
        painter.drawText(label_rect, Qt.AlignmentFlag.AlignCenter, self._label)

        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.setPen(QPen(border, 2.0 if selected else 1.0))
        painter.drawRoundedRect(outer, 3.0, 3.0)


class Filmstrip(QFrame):
    frame_selected = pyqtSignal(int)

    def __init__(self, page_size=6, parent=None):
        super().__init__(parent)
        self.setObjectName("filmstripBar")
        self.setFixedHeight(108)
        self.max_page_size = int(page_size)
        self.page_size = self.max_page_size
        self.page = 0
        self.paths: list[str] = []
        self.thumbnails: dict[int, QPixmap] = {}
        self.selected_index = 0
        root = QHBoxLayout(self)
        root.setContentsMargins(4, 4, 4, 4)
        root.setSpacing(4)
        self.btn_previous = QToolButton()
        self.btn_previous.setText("‹")
        self.btn_previous.setFixedWidth(26)
        self.btn_previous.setToolTip("Önceki kare sayfası")
        self.btn_previous.clicked.connect(lambda: self.change_page(-1))
        root.addWidget(self.btn_previous)
        self.buttons = []
        for slot in range(self.page_size):
            button = ThumbnailCard()
            button.setObjectName("thumbnailButton")
            button.clicked.connect(lambda _checked=False, position=slot: self._slot_clicked(position))
            self.buttons.append(button)
            root.addWidget(button, 1)
        self.btn_next = QToolButton()
        self.btn_next.setText("›")
        self.btn_next.setFixedWidth(26)
        self.btn_next.setToolTip("Sonraki kare sayfası")
        self.btn_next.clicked.connect(lambda: self.change_page(1))
        root.addWidget(self.btn_next)
        self.refresh()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        available = max(1, self.width() - 72)
        responsive_size = max(1, min(self.max_page_size, available // 145))
        if responsive_size != self.page_size:
            self.page_size = responsive_size
            self.page = self.selected_index // self.page_size
            self.refresh()

    def set_frames(self, paths: list[str]):
        self.paths = list(paths)
        self.page = 0
        self.selected_index = 0
        self.thumbnails.clear()
        self.refresh()

    def set_thumbnail(self, index: int, pixmap: QPixmap):
        self.thumbnails[int(index)] = pixmap
        if self.page * self.page_size <= index < (self.page + 1) * self.page_size:
            self.refresh()

    def select(self, index: int):
        if not self.paths:
            return
        self.selected_index = max(0, min(int(index), len(self.paths) - 1))
        self.page = self.selected_index // self.page_size
        self.refresh()

    def change_page(self, direction: int):
        max_page = max(0, (len(self.paths) - 1) // self.page_size)
        self.page = max(0, min(max_page, self.page + int(direction)))
        self.refresh()

    def refresh(self):
        start = self.page * self.page_size
        for slot, button in enumerate(self.buttons):
            if slot >= self.page_size:
                button.setVisible(False)
                continue
            index = start + slot
            if index >= len(self.paths):
                button.setVisible(False)
                continue
            button.setVisible(True)
            name = Path(self.paths[index]).name
            button.set_label(f"{index + 1} / {len(self.paths)}")
            button.setToolTip(f"{name}\n{self.paths[index]}")
            pixmap = self.thumbnails.get(index)
            button.set_thumbnail(pixmap)
            button.setProperty("selected", index == self.selected_index)
            button.update()
        max_page = max(0, (len(self.paths) - 1) // self.page_size)
        self.btn_previous.setEnabled(self.page > 0)
        self.btn_next.setEnabled(self.page < max_page)

    def _slot_clicked(self, slot: int):
        index = self.page * self.page_size + slot
        if index < len(self.paths):
            self.selected_index = index
            self.refresh()
            self.frame_selected.emit(index)


def section(title: str, parent=None) -> tuple[QFrame, QVBoxLayout]:
    frame = QFrame(parent)
    frame.setObjectName("section")
    layout = QVBoxLayout(frame)
    layout.setContentsMargins(8, 7, 8, 8)
    layout.setSpacing(6)
    label = QLabel(title)
    label.setObjectName("sectionTitle")
    layout.addWidget(label)
    return frame, layout


def metric_row(label: str, value: str = "--") -> tuple[QWidget, QLabel]:
    widget = QWidget()
    layout = QHBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(8)
    left = QLabel(label)
    left.setObjectName("muted")
    right = QLabel(value)
    right.setObjectName("value")
    right.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    layout.addWidget(left)
    layout.addStretch(1)
    layout.addWidget(right)
    return widget, right


def parse_log_payload(payload: str) -> tuple[str, str, str]:
    parts = str(payload).split("|", 2)
    if len(parts) == 3 and parts[0] in {"INFO", "WARN", "ERROR"}:
        return parts[0], parts[1], parts[2]
    lowered = str(payload).lower()
    if "hata" in lowered or "error" in lowered or "başarısız" in lowered:
        return "ERROR", "System", str(payload)
    if "uyarı" in lowered or "warn" in lowered:
        return "WARN", "System", str(payload)
    return "INFO", "System", str(payload)
