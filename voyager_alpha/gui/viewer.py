from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QFrame, QVBoxLayout


class FitsViewer(QFrame):
    image_clicked = pyqtSignal(float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("imageViewport")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.plot_widget = pg.PlotWidget(background="#050709")
        self.plot_widget.setMenuEnabled(False)
        self.plot_widget.hideAxis("left")
        self.plot_widget.hideAxis("bottom")
        self.plot_widget.setAspectLocked(True)
        self.plot_widget.setMouseEnabled(x=True, y=True)
        self.plot_widget.getViewBox().setDefaultPadding(0.01)
        self.image_item = pg.ImageItem(axisOrder="row-major")
        self.plot_widget.addItem(self.image_item)
        self.plot_widget.scene().sigMouseClicked.connect(self._scene_clicked)
        layout.addWidget(self.plot_widget)

        self.raw_data = None
        self.image_shape = None
        self.display_data = None
        self.sequence_levels: tuple[float, float] | None = None
        self._first_image = True
        self._overlay_items = []
        self._grid_items = []

        self.known_scatter = pg.ScatterPlotItem(
            size=15,
            pen=pg.mkPen(QColor("#4dc56a"), width=1.3),
            brush=pg.mkBrush(0, 0, 0, 0),
            symbol="o",
        )
        self.candidate_scatter = pg.ScatterPlotItem(
            size=14,
            pen=pg.mkPen(QColor("#ef5f51"), width=1.35),
            brush=pg.mkBrush(239, 95, 81, 22),
            symbol="s",
        )
        self.track_scatter = pg.ScatterPlotItem(
            size=11,
            pen=pg.mkPen(QColor("#2dd7df"), width=1.4),
            brush=pg.mkBrush(45, 215, 223, 40),
            symbol="o",
        )
        self.plot_widget.addItem(self.known_scatter)
        self.plot_widget.addItem(self.candidate_scatter)
        self.plot_widget.addItem(self.track_scatter)

        self.empty_label = pg.TextItem(
            html='<div style="color:#71818a;font-size:13px;">FITS sekansı bekleniyor</div>',
            anchor=(0.5, 0.5),
        )
        self.plot_widget.addItem(self.empty_label)
        self.legend_label = pg.TextItem(
            html=(
                '<div style="background:rgba(6,9,11,185);padding:5px 7px;color:#cbd4d8;font-size:10px;">'
                '<span style="color:#4dc56a;">○ Bilinen</span>&nbsp;&nbsp;'
                '<span style="color:#ef5f51;">□ Tek-kare residual</span>&nbsp;&nbsp;'
                '<span style="color:#2dd7df;">— Tracklet</span></div>'
            ),
            anchor=(0, 0),
        )
        self.plot_widget.addItem(self.legend_label)
        self.legend_label.setVisible(False)

    def set_image(
        self,
        data,
        *,
        stretch_mode="Auto STF",
        sequence_levels=None,
        manual_levels=None,
        invert=False,
        keep_view=True,
    ):
        array = np.asarray(data, dtype=np.float32)
        self.raw_data = array
        self.image_shape = array.shape
        display, levels = stretch_image(
            array,
            mode=stretch_mode,
            sequence_levels=sequence_levels,
            manual_levels=manual_levels,
            invert=invert,
        )
        self.display_data = display
        self.image_item.setImage(display, autoLevels=False, levels=(0.0, 1.0))
        self.empty_label.setVisible(False)
        self.legend_label.setVisible(True)
        self.legend_label.setPos(8, 8)
        if self._first_image or not keep_view:
            self.auto_range()
            self._first_image = False
        return levels

    def set_display_image(self, display, original_shape=None, *, keep_view=True):
        array = np.asarray(display)
        if array.dtype == np.uint8:
            array = array.astype(np.float32) / 255.0
        else:
            array = np.asarray(array, dtype=np.float32)
        self.display_data = array
        shape = tuple(original_shape) if original_shape is not None else array.shape
        self.raw_data = None
        self.image_shape = shape
        self.image_item.setImage(array, autoLevels=False, levels=(0.0, 1.0))
        if array.shape != shape:
            self.image_item.setRect(QRectF(0.0, 0.0, float(shape[1]), float(shape[0])))
        self.empty_label.setVisible(False)
        self.legend_label.setVisible(True)
        self.legend_label.setPos(8, 8)
        if self._first_image or not keep_view:
            self.auto_range()
            self._first_image = False

    def clear_image(self):
        self.raw_data = None
        self.image_shape = None
        self.display_data = None
        self.image_item.clear()
        self.clear_overlays()
        self.empty_label.setVisible(True)
        self.empty_label.setPos(0, 0)
        self.legend_label.setVisible(False)

    def set_overlays(
        self,
        *,
        known=None,
        candidates=None,
        selected_tracklet=None,
        show_known=True,
        show_candidates=True,
        show_track=True,
    ):
        known = known or []
        candidates = candidates or []
        self._clear_dynamic_items()
        known_points = []
        if show_known:
            for item in known:
                x, y = _xy(item)
                if x is not None and y is not None:
                    known_points.append({"pos": (x, y)})
        candidate_points = []
        if show_candidates:
            for item in candidates:
                x, y = _xy(item)
                if x is not None and y is not None:
                    candidate_points.append({"pos": (x, y)})
        self.known_scatter.setData(known_points)
        self.candidate_scatter.setData(candidate_points)

        track_points = []
        if show_track and selected_tracklet is not None:
            detections = sorted(selected_tracklet.detections, key=lambda item: item.frame_index)
            xs = [float(item.x) for item in detections]
            ys = [float(item.y) for item in detections]
            if len(xs) >= 2:
                line = pg.PlotDataItem(
                    xs,
                    ys,
                    pen=pg.mkPen(QColor("#2dd7df"), width=1.6),
                    connect="all",
                )
                self.plot_widget.addItem(line)
                self._overlay_items.append(line)
            track_points = [{"pos": (x, y)} for x, y in zip(xs, ys)]
            if xs:
                label = pg.TextItem(
                    html=(
                        '<div style="background:rgba(6,9,11,205);padding:2px 5px;'
                        f'color:#2dd7df;font-size:10px;">{selected_tracklet.tracklet_id}</div>'
                    ),
                    anchor=(0, 1),
                )
                label.setPos(xs[-1] + 5, ys[-1] - 5)
                self.plot_widget.addItem(label)
                self._overlay_items.append(label)
        self.track_scatter.setData(track_points)

    def set_grid_visible(self, visible: bool):
        for item in self._grid_items:
            self.plot_widget.removeItem(item)
        self._grid_items.clear()
        if not visible or self.image_shape is None:
            return
        height, width = self.image_shape
        pen = pg.mkPen((45, 162, 177, 72), width=0.8)
        for fraction in (0.2, 0.4, 0.6, 0.8):
            vertical = pg.InfiniteLine(pos=width * fraction, angle=90, pen=pen, movable=False)
            horizontal = pg.InfiniteLine(pos=height * fraction, angle=0, pen=pen, movable=False)
            self.plot_widget.addItem(vertical)
            self.plot_widget.addItem(horizontal)
            self._grid_items.extend((vertical, horizontal))

    def clear_overlays(self):
        self.known_scatter.clear()
        self.candidate_scatter.clear()
        self.track_scatter.clear()
        self._clear_dynamic_items()

    def auto_range(self):
        self.plot_widget.getViewBox().autoRange(padding=0.01)

    def zoom_by(self, factor: float):
        self.plot_widget.getViewBox().scaleBy((factor, factor))

    def one_to_one(self):
        if self.image_shape is None:
            return
        height, width = self.image_shape
        view = self.plot_widget.getViewBox()
        view.setRange(xRange=(0, width), yRange=(0, height), padding=0)

    def center_on(self, x: float, y: float, radius: float = 70.0):
        self.plot_widget.getViewBox().setRange(
            xRange=(x - radius, x + radius),
            yRange=(y - radius, y + radius),
            padding=0,
        )

    def _scene_clicked(self, event):
        if event.double():
            self.auto_range()
            return
        if event.button() != Qt.MouseButton.LeftButton:
            return
        point = self.plot_widget.getViewBox().mapSceneToView(event.scenePos())
        if self.image_shape is None:
            return
        if 0 <= point.x() < self.image_shape[1] and 0 <= point.y() < self.image_shape[0]:
            self.image_clicked.emit(float(point.x()), float(point.y()))

    def _clear_dynamic_items(self):
        for item in self._overlay_items:
            self.plot_widget.removeItem(item)
        self._overlay_items.clear()


def stretch_image(
    data: np.ndarray,
    *,
    mode: str,
    sequence_levels=None,
    manual_levels=None,
    invert=False,
) -> tuple[np.ndarray, tuple[float, float]]:
    array = np.asarray(data, dtype=np.float32)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return np.zeros_like(array, dtype=np.float32), (0.0, 1.0)
    normalized_mode = mode.lower()
    if manual_levels is not None and "manual" in normalized_mode:
        black_percentile, midtone, white_percentile = manual_levels
        black_percentile = float(np.clip(black_percentile, 0.0, 99.999))
        white_percentile = float(np.clip(white_percentile, black_percentile + 0.001, 100.0))
        lo, hi = np.nanpercentile(finite, (black_percentile, white_percentile))
    else:
        if sequence_levels is not None:
            lo, hi = sequence_levels
        elif "linear" in normalized_mode:
            lo, hi = np.nanpercentile(finite, (0.1, 99.9))
        else:
            median = float(np.nanmedian(finite))
            mad = float(np.nanmedian(np.abs(finite - median)))
            sigma = max(1.4826 * mad, 1e-6)
            lo = max(float(np.nanpercentile(finite, 0.02)), median - 2.2 * sigma)
            hi = min(float(np.nanpercentile(finite, 99.98)), median + 14.0 * sigma)
            if hi <= lo:
                lo, hi = np.nanpercentile(finite, (0.1, 99.9))
        midtone = 0.34
    hi = max(float(hi), float(lo) + 1e-9)
    scaled = np.clip((array - float(lo)) / (hi - float(lo)), 0.0, 1.0)
    if "linear" in normalized_mode:
        display = scaled
    elif "asinh" in normalized_mode:
        display = np.arcsinh(scaled * 7.0) / np.arcsinh(7.0)
    else:
        gamma = 2.0 ** ((float(midtone) - 0.5) * 4.0)
        display = np.power(scaled, gamma)
    if invert:
        display = 1.0 - display
    return np.nan_to_num(display, copy=False).astype(np.float32), (float(lo), float(hi))


def sequence_stf_levels(data: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(data, dtype=np.float32)
    finite = finite[np.isfinite(finite)]
    median = float(np.median(finite))
    mad = float(np.median(np.abs(finite - median)))
    sigma = max(1.4826 * mad, 1e-6)
    lo = max(float(np.percentile(finite, 0.02)), median - 2.2 * sigma)
    hi = min(float(np.percentile(finite, 99.98)), median + 14.0 * sigma)
    if hi <= lo:
        lo, hi = np.percentile(finite, (0.1, 99.9))
    return float(lo), float(hi)


def _xy(item):
    if isinstance(item, dict):
        return item.get("x"), item.get("y")
    return getattr(item, "x", None), getattr(item, "y", None)
