import logging
import os
import re

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem

from calculations import dq_temp_signal
from controllers.base_tab_controller import BaseTabController
from controllers.table_columns import DQTempColumns
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


class DQTempTabController(BaseTabController):
    def connect_signals(self):
        self.ui.DQTemp_Table_Results.horizontalHeader().sectionDoubleClicked.connect(
            lambda index: self.parent.renameSection(self.ui.DQTemp_Table_Results, index)
        )

    def update_DQ_comparison(self):
        if not self.parent.selected_DQfiles:
            self._warn_no_data("No DQ temperature data available. Select DQ comparison files first.")
            return

        with busy_cursor():
            table = self.ui.DQTemp_Table_Results
            valid_files = []
            failed_files = []
            logger.info("DQ Temp loading %d files", len(self.parent.selected_DQfiles))
            table.setRowCount(0)

            for file_path in self.parent.selected_DQfiles:
                filename = os.path.basename(file_path)
                logger.info(
                    "DQ Temp processing file %d/%d: %s",
                    len(valid_files) + len(failed_files) + 1,
                    len(self.parent.selected_DQfiles),
                    filename,
                )

                try:
                    dq_temp_signal.load_distribution_file(file_path)
                except Exception:
                    logger.exception("DQ Temp invalid input file: %s", file_path)
                    failed_files.append(filename)
                    continue

                row = table.rowCount()
                table.insertRow(row)
                valid_files.append(file_path)
                table.setItem(row, DQTempColumns.FOLDER, QTableWidgetItem(filename))
                table.setItem(row, DQTempColumns.NAME, self._x_axis_item(filename, row))

            self.parent.selected_DQfiles = valid_files
            table.resizeColumnsToContents()

            if failed_files:
                self._warn_failed_files(
                    "Some DQ temperature files could not be loaded and were skipped.",
                    failed_files,
                )

            if not valid_files:
                self._clear_plots()
                return

            self.launch()

    def launch(self):
        if not self.parent.selected_DQfiles:
            self._warn_no_data("No DQ temperature data available. Select DQ comparison files first.")
            return

        logger.info("DQ Temp launching comparison fit")
        self.parent.dq_t2 = {}
        failed_files = []

        for row, file_path in enumerate(self.parent.selected_DQfiles):
            try:
                self.parent.dq_t2[row] = dq_temp_signal.load_distribution_file(file_path)
            except Exception:
                logger.exception("DQ Temp launch failed for file: %s", file_path)
                self.parent.dq_t2[row] = None
                failed_files.append(os.path.basename(file_path))

        plot_failed_files = self.update_DQ_comparison_plot(show_warning=False)
        all_failed_files = list(dict.fromkeys(failed_files + plot_failed_files))
        if all_failed_files:
            self._warn_failed_files(
                "Some DQ temperature files could not be fitted. They were filled with NaN.",
                all_failed_files,
            )

    def update_DQ_comparison_plot(self, show_warning=True):
        table = self.ui.DQTemp_Table_Results
        if table.rowCount() == 0:
            self._clear_plots()
            if show_warning:
                self._warn_no_data("Cannot update DQ temperature plot because the table is empty.")
            return []

        if not self.parent.dq_t2:
            self._clear_plots()
            if show_warning:
                self._warn_no_data("No DQ temperature data available. Select DQ comparison files first.")
            return []

        cmap = pg.ColorMap([0, len(self.parent.dq_t2)], [pg.mkColor("b"), pg.mkColor("r")])
        self.parent.dq_comparison_distribution = {
            "File name": [],
            "X axis": [],
            "Center": [],
            "FWHM": [],
            "Lorentz ratio": [],
            "Fitting type": [],
            "T2 limit": [],
        }
        self._clear_plots()
        legend = self.ui.DQTemp_PlotWidget_T2Distribution.addLegend(offset=(0, 0))
        legend_center = self.ui.DQTemp_PlotWidget_CenterVsXAxis.addLegend()

        if legend is not None:
            legend.clear()
            legend.setPen((0, 0, 0))
            legend.anchor(itemPos=(1, 0), parentPos=(1, 0))

        if legend_center is not None:
            legend_center.clear()
            legend_center.setPen((0, 0, 0))

        comparison_axis = []
        center_gauss = []
        center_lorenz = []
        center_voigt = []
        center_derivative = []
        failed_files = []

        for row, (key, data) in zip(range(table.rowCount()), self.parent.dq_t2.items()):
            file_name_item = table.item(row, DQTempColumns.NAME)
            file_item = table.item(row, DQTempColumns.FOLDER)
            file_name = file_name_item.text() if file_name_item is not None else str(row + 1)
            display_name = file_item.text() if file_item is not None else file_name

            if file_name == "hide":
                continue

            comparison_value = self._comparison_value(file_name, row)
            comparison_axis.append(comparison_value)

            if data is None:
                self._set_failed_fit_values(row)
                self._append_nan_fit_values(center_gauss, center_lorenz, center_voigt, center_derivative)
                failed_files.append(display_name)
                continue

            try:
                fit_result = dq_temp_signal.fit_distribution(data)
            except Exception:
                logger.exception("DQ Temp fit failed for row %d: %s", row, display_name)
                self._set_failed_fit_values(row)
                self._append_nan_fit_values(center_gauss, center_lorenz, center_voigt, center_derivative)
                failed_files.append(display_name)
                continue

            center_gauss.append(fit_result["center_gauss"])
            center_lorenz.append(fit_result["center_lorenz"])
            center_voigt.append(fit_result["center_voigt"])
            center_derivative.append(fit_result["center_derivative"])
            self._plot_distribution(key, fit_result, file_name, cmap)
            self._write_fit_values(row, fit_result)

        table.resizeColumnsToContents()
        self._plot_centers(comparison_axis, center_gauss, center_lorenz, center_voigt, center_derivative)

        if failed_files and show_warning:
            self._warn_failed_files(
                "Some DQ temperature files could not be fitted. They were filled with NaN.",
                failed_files,
            )

        return failed_files

    def _x_axis_item(self, filename, row):
        item_name = QTableWidgetItem()
        pattern = r"Table_DQ_(-?[0-9]+).*.csv"
        match = re.search(pattern, filename)
        if match is not None:
            item_name.setData(Qt.EditRole, float(match.group(1)))
        else:
            item_name.setData(Qt.EditRole, float(row + 1))

        return item_name

    def _comparison_value(self, file_name, row):
        try:
            return float(file_name)
        except Exception:
            self.ui.DQTemp_Table_Results.setItem(row, DQTempColumns.NAME, QTableWidgetItem("NaN"))
            return np.nan

    def _plot_distribution(self, key, fit_result, file_name, cmap):
        color = tuple(cmap.map(key))
        self.ui.DQTemp_PlotWidget_T2Distribution.plot(
            fit_result["t2_linearized"],
            fit_result["dq_normalized"],
            pen=None,
            symbolPen=None,
            symbol="o",
            symbolBrush=color,
            symbolSize=10,
            name=file_name,
        )
        self.ui.DQTemp_PlotWidget_T2Distribution.plot(
            fit_result["t2_fit"],
            fit_result["voigt_fit"],
            pen=color,
        )
        self.ui.DQTemp_PlotWidget_PolyFit.plot(
            fit_result["t2_linearized"],
            fit_result["dq_normalized"],
            pen=None,
            symbolPen=None,
            symbol="o",
            symbolBrush=color,
            symbolSize=10,
        )
        self.ui.DQTemp_PlotWidget_PolyFit.plot(
            fit_result["derivative_x"],
            fit_result["derivative_y"],
            pen=color,
        )

    def _plot_centers(self, comparison_axis, center_gauss, center_lorenz, center_voigt, center_derivative):
        self._plot_center_series(comparison_axis, center_gauss, "r", "Gaus")
        self._plot_center_series(comparison_axis, center_lorenz, "b", "Lorenz")
        self._plot_center_series(comparison_axis, center_voigt, "k", "Voigt")
        self._plot_center_series(comparison_axis, center_derivative, "g", "Derivative")

    def _plot_center_series(self, comparison_axis, center_values, color, name):
        x_values, y_values = dq_temp_signal.finite_xy(comparison_axis, center_values)
        if len(x_values) == 0:
            return

        self.ui.DQTemp_PlotWidget_CenterVsXAxis.plot(
            x_values,
            y_values,
            pen=color,
            symbolPen=None,
            symbol="o",
            symbolBrush=color,
            name=name,
        )

    def _write_fit_values(self, row, fit_result):
        self._set_num(row, DQTempColumns.CENTER_GAUSS, fit_result["center_gauss"])
        self._set_num(row, DQTempColumns.CENTER_LORENZ, fit_result["center_lorenz"])
        self._set_num(row, DQTempColumns.CENTER_VOIGT, fit_result["center_voigt"])
        self._set_num(row, DQTempColumns.CENTER_Y, fit_result["center_derivative"])
        self._set_num(row, DQTempColumns.FWHM_GAUSS, fit_result["fwhm_gauss"])
        self._set_num(row, DQTempColumns.FWHM_LORENZ, fit_result["fwhm_lorenz"])
        self._set_num(row, DQTempColumns.FWHM_VOIGT, fit_result["fwhm_voigt"])

    def _set_failed_fit_values(self, row):
        for column in (
            DQTempColumns.CENTER_GAUSS,
            DQTempColumns.CENTER_LORENZ,
            DQTempColumns.CENTER_VOIGT,
            DQTempColumns.CENTER_Y,
            DQTempColumns.FWHM_GAUSS,
            DQTempColumns.FWHM_LORENZ,
            DQTempColumns.FWHM_VOIGT,
        ):
            self._set_num(row, column, np.nan)

    def _set_num(self, row, column, value):
        item = QTableWidgetItem()
        item.setData(Qt.EditRole, round(float(value), 2) if np.isfinite(value) else np.nan)
        self.ui.DQTemp_Table_Results.setItem(row, column, item)

    def _append_nan_fit_values(self, center_gauss, center_lorenz, center_voigt, center_derivative):
        center_gauss.append(np.nan)
        center_lorenz.append(np.nan)
        center_voigt.append(np.nan)
        center_derivative.append(np.nan)

    def _clear_plots(self):
        self.ui.DQTemp_PlotWidget_T2Distribution.clear()
        self.ui.DQTemp_PlotWidget_CenterVsXAxis.clear()
        self.ui.DQTemp_PlotWidget_PolyFit.clear()

    def _warn_no_data(self, message):
        QMessageBox.warning(self.parent, "No DQ temperature data", message, QMessageBox.Ok)

    def _warn_failed_files(self, message, failed_files):
        preview = "\n".join(failed_files[:5])
        if len(failed_files) > 5:
            preview += f"\n...and {len(failed_files) - 5} more"

        QMessageBox.warning(
            self.parent,
            "DQ temperature processing warning",
            f"{message}\n\n{preview}",
            QMessageBox.Ok,
        )
