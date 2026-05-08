"""Dialogs used by the Reconstruct FID tab."""

from __future__ import annotations

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)
import matplotlib.pyplot as plt


class RecFIDManualEchoTimeDialog(QDialog):
    """Manual echo-time assignment for SE file lists."""

    def __init__(self, files, echo_times=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Manual Echo Time")
        self.resize(700, 400)
        echo_times = echo_times or []

        self.table = QTableWidget(len(files), 2, self)
        self.table.setHorizontalHeaderLabels(["File", "Echo Time, μs"])
        for row, file_path in enumerate(files):
            file_item = QTableWidgetItem(str(file_path))
            file_item.setFlags(file_item.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(row, 0, file_item)
            value = echo_times[row] if row < len(echo_times) else ""
            self.table.setItem(row, 1, QTableWidgetItem(str(value)))
        self.table.resizeColumnsToContents()

        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout = QHBoxLayout()
        button_layout.addStretch(1)
        button_layout.addWidget(ok_button)
        button_layout.addWidget(cancel_button)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Enter echo times for each loaded SE file."))
        layout.addWidget(self.table)
        layout.addLayout(button_layout)

    def echo_times(self):
        values = []
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 1)
            if item is None or item.text().strip() == "":
                raise ValueError("Every SE file needs an echo time.")
            values.append(float(item.text()))
        return values


class RecFIDTimeAnalysisDialog(QDialog):
    """Optional T2 map visualization for RecFID time-range analysis."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.canvas = FigureCanvas(plt.figure())
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(100)
        self.slider.setValue(100)
        self.slider.valueChanged.connect(self.update_plot)
        self.slider_label = QLabel("T₂ Threshold: 100%")
        self.rainbow_checkbox = QCheckBox("Rainbow")
        self.rainbow_checkbox.stateChanged.connect(self.update_plot)
        self.save_button = QPushButton("Save Data")
        self.save_button.clicked.connect(self.save_data)

        controls = QHBoxLayout()
        controls.addWidget(self.slider_label)
        controls.addWidget(self.slider)
        controls.addWidget(self.rainbow_checkbox)
        controls.addWidget(self.save_button)

        layout = QVBoxLayout(self)
        layout.addWidget(self.canvas)
        layout.addLayout(controls)

        self.start_range = None
        self.finish_range = None
        self.t2 = None
        self.m2 = None

    def plot_data(self, start_range, finish_range, t2, m2):
        self.start_range = np.asarray(start_range)
        self.finish_range = np.asarray(finish_range)
        self.t2 = np.asarray(t2)
        self.m2 = np.asarray(m2)
        self.update_plot()

    def update_plot(self):
        if self.t2 is None or self.m2 is None:
            return
        max_t2 = np.nanmax(self.t2)
        if not np.isfinite(max_t2) or max_t2 == 0:
            max_t2 = 1
        threshold = self.slider.value() / 100.0 * max_t2
        self.slider_label.setText(f"T₂ Threshold: {threshold:.2f}")
        t2_masked = np.ma.masked_where((self.t2 < 5) | (self.t2 > threshold), self.t2)
        cmap = "rainbow" if self.rainbow_checkbox.isChecked() else "RdYlGn"
        self.canvas.figure.clf()
        gs = self.canvas.figure.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.1)
        ax_main = self.canvas.figure.add_subplot(gs[0])
        ax_hist = self.canvas.figure.add_subplot(gs[1])
        c = ax_main.pcolormesh(self.finish_range, self.start_range, t2_masked, shading="auto", cmap=cmap)
        self.canvas.figure.colorbar(c, ax=ax_main, label="T₂ Value")
        ax_main.set_xlabel("Finish Range")
        ax_main.set_ylabel("Start Range")
        ax_main.set_xlim(self.finish_range[0], self.finish_range[-1])
        ax_main.set_ylim(self.start_range[0], self.start_range[-1])
        compressed = t2_masked.compressed()
        if compressed.size:
            ax_hist.hist(compressed, bins=20, color="blue")
            ax_hist.set_xlim(float(np.min(compressed)), float(np.max(compressed)))
        ax_hist.set_ylabel("Count")
        ax_hist.yaxis.set_label_position("right")
        ax_hist.yaxis.tick_right()
        ax_hist.set_xlabel("T₂ value")
        ax_main.set_title("T₂ Distribution")
        self.canvas.draw()

    def save_data(self):
        if self.start_range is None or self.finish_range is None or self.t2 is None or self.m2 is None:
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Data", "", "CSV Files (*.csv);;Text Files (*.txt)"
        )
        if not file_path:
            return
        data = np.column_stack(
            (
                self.start_range.repeat(len(self.finish_range)),
                np.tile(self.finish_range, len(self.start_range)),
                self.m2.flatten(),
                self.t2.flatten(),
            )
        )
        np.savetxt(file_path, data, delimiter=",", header="Start Range,Finish Range,M2Value,T2 Value", comments="")
