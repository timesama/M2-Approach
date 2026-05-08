"""Dialogs used by the Reconstruct FID tab."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QCheckBox, QDialog, QFileDialog, QHBoxLayout, QLabel, QPushButton, QSlider, QVBoxLayout


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
        """Render the supplied M2/T2 grids into the dialog canvas."""
        self.start_range = np.asarray(start_range)
        self.finish_range = np.asarray(finish_range)
        self.t2 = np.asarray(t2)
        self.m2 = np.asarray(m2)
        self.update_plot()

    def update_plot(self):
        """Refresh the heatmap using the current threshold and color-map settings."""
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
        grid = self.canvas.figure.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.1)
        ax_main = self.canvas.figure.add_subplot(grid[0])
        ax_hist = self.canvas.figure.add_subplot(grid[1])
        heatmap = ax_main.pcolormesh(self.finish_range, self.start_range, t2_masked, shading="auto", cmap=cmap)
        self.canvas.figure.colorbar(heatmap, ax=ax_main, label="T₂ Value")
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
        """Export the currently displayed time-analysis grid."""
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

    def closeEvent(self, event):
        plt.close(self.canvas.figure)
        self.deleteLater()
        event.accept()
