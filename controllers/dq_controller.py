import re

import numpy as np
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem

import Calculator as Cal
from calculations import dq_signal
from controllers.base_tab_controller import BaseTabController


class DQTabController(BaseTabController):
    def connect_signals(self):
        self.ui.DQ_RadioButton_LogScale.clicked.connect(lambda *_args: self.plot_fit_from_user())
        self.ui.DQ_Table_Data.currentItemChanged.connect(lambda *_args: self.update_graphs())
        self.ui.DQ_Table_Data.itemSelectionChanged.connect(self.update_graphs)
        self.ui.DQ_ComboBox_FitFunction.activated.connect(lambda *_args: self.plot_fit_from_user())
        self.ui.DQ_DoubleSpinBox_FilterFrom.editingFinished.connect(self.update_graphs_from_user)
        self.ui.DQ_DoubleSpinBox_FilterTo.editingFinished.connect(self.update_graphs_from_user)

    def process_processed_file(self, i, filename, x, y, z, m2, t2, file_path):
        match = re.search(r"_(\d+\.\d+)_", filename)
        dq_time = match.group(1) if match else "0"
        amplitude = Cal._calculate_amplitude(y, z)
        dq = Cal.calculate_DQ_intensity(x, amplitude)
        self.ui.DQ_Table_Data.setRowCount(i)
        self.parent.fill_table(self.ui.DQ_Table_Data, dq_time, dq, m2, t2, i)

        if self.ui.radioButton.isChecked():
            self.parent.save_figures(file_path, dq_time)

        self.update_graphs()

    def update_graphs_from_user(self):
        if not self._has_table_data("Cannot update DQ plot because the table is empty."):
            return

        self.update_graphs(show_warning=True)

    def update_graphs(self, show_warning=False):
        table = self.ui.DQ_Table_Data
        if table.rowCount() == 0:
            self._clear_plots()
            if show_warning:
                self._warn_no_data("Cannot update DQ plot because the table is empty.")
            return

        dq_time = self._read_numeric_column(0, show_warning=show_warning)
        t2_values = self._read_numeric_column(3, show_warning=show_warning)
        if dq_time is None or t2_values is None:
            return

        if len(dq_time) > 1 and len(t2_values) > 1:
            if self.linearization(show_warning=show_warning):
                self.plot_fit(show_warning=False)
        else:
            self.dq_t2_graph(show_warning=show_warning)

    def dq_t2_graph(self, show_warning=False):
        if not self._has_table_data("Cannot update DQ plot because the table is empty.", show_warning):
            return

        x = self._read_numeric_column(0, show_warning=show_warning)
        y = self._read_numeric_column(3, show_warning=show_warning)
        if x is None or y is None:
            return

        self.ui.DQ_PlotWidget_T2.clear()
        self.ui.DQ_PlotWidget_T2.plot(
            x,
            y,
            pen=None,
            symbol="o",
            symbolPen=None,
            symbolBrush=(255, 0, 0, 255),
            symbolSize=10,
        )
        self.highlight_selected_point_widget_1()

    def linearization(self, show_warning=False):
        if not self._has_table_data("No DQ data available. Load or analyze DQ files first.", show_warning):
            return False

        time_min = self.ui.DQ_DoubleSpinBox_FilterFrom.value()
        time_max = self.ui.DQ_DoubleSpinBox_FilterTo.value()
        dq_time = self._read_numeric_column(0, show_warning=show_warning)
        t2 = self._read_numeric_column(3, show_warning=show_warning)
        dq = self._read_numeric_column(1, show_warning=show_warning)
        if dq_time is None or t2 is None or dq is None:
            return False

        result = dq_signal.calculate_linearization(
            np.array(dq_time),
            np.array(t2),
            np.array(dq),
            time_min,
            time_max,
        )
        if result is None:
            if show_warning:
                QMessageBox.warning(
                    self.parent,
                    "Invalid DQ range",
                    "Cannot fit DQ data because the selected range is invalid.",
                    QMessageBox.Ok,
                )
            return False

        if result["time_min"] != time_min:
            self.ui.DQ_DoubleSpinBox_FilterFrom.setValue(result["time_min"])
        if result["time_max"] != time_max:
            self.ui.DQ_DoubleSpinBox_FilterTo.setValue(result["time_max"])

        table = self.ui.DQ_Table_Data
        table.setColumnCount(5)
        table.setColumnWidth(4, 70)
        table.setHorizontalHeaderItem(4, QTableWidgetItem("T₂* lin"))
        table.setColumnCount(6)
        table.setColumnWidth(5, 70)
        table.setHorizontalHeaderItem(5, QTableWidgetItem("DQ Norm"))

        for row in range(table.rowCount()):
            t2_linearized = round(result["t2_linearized"][row], 4)
            dq_norm = round(result["dq_norm"][row], 4)
            table.setItem(row, 4, QTableWidgetItem(str(t2_linearized)))
            table.setItem(row, 5, QTableWidgetItem(str(dq_norm)))

        self.dq_t2_graph(show_warning=False)
        self.ui.DQ_PlotWidget_T2.plot(result["x_line"], result["y_line"], pen="r")
        self.t2_dq_graph(show_warning=False)
        return True

    def t2_dq_graph(self, show_warning=False):
        if not self._has_linearized_columns(show_warning=show_warning):
            return

        x = self._read_numeric_column(4, show_warning=show_warning)
        y = self._read_numeric_column(5, show_warning=show_warning)
        if x is None or y is None:
            return

        new_x, axis_label = dq_signal.prepare_t2_axis(
            np.array(x),
            self.ui.DQ_RadioButton_LogScale.isChecked(),
        )
        self.ui.DQ_PlotWidget_NormIntensity.getAxis("bottom").setLabel(axis_label)
        self.ui.DQ_PlotWidget_NormIntensity.clear()
        self.ui.DQ_PlotWidget_NormIntensity.plot(
            new_x,
            y,
            pen=None,
            symbol="o",
            symbolPen=None,
            symbolBrush=(255, 0, 0, 255),
            symbolSize=10,
        )
        self.highlight_selected_point_widget_2()

    def plot_fit_from_user(self):
        self.plot_fit(show_warning=True)

    def plot_fit(self, show_warning=False):
        if not self._has_table_data("No DQ data available. Load or analyze DQ files first.", show_warning):
            return

        if not self._has_linearized_columns(show_warning=False):
            if self.linearization(show_warning=show_warning):
                self.t2_dq_graph(show_warning=False)
            else:
                return

        x_values = self._read_numeric_column(4, show_warning=show_warning)
        y_values = self._read_numeric_column(5, show_warning=show_warning)
        if x_values is None or y_values is None:
            return

        fit_function = self.ui.DQ_ComboBox_FitFunction.currentText()
        if not fit_function:
            if show_warning:
                QMessageBox.warning(
                    self.parent,
                    "No DQ fit function",
                    "Select a DQ fit function first.",
                    QMessageBox.Ok,
                )
            return

        try:
            result = dq_signal.fit_distribution(
                np.array(x_values),
                np.array(y_values),
                fit_function,
                self.ui.DQ_RadioButton_LogScale.isChecked(),
            )
        except (RuntimeError, ValueError, TypeError):
            if show_warning:
                QMessageBox.warning(
                    self.parent,
                    "DQ fit error",
                    "Cannot calculate DQ fit because required columns contain non-numeric values.",
                    QMessageBox.Ok,
                )
            return

        if result is None:
            return

        self.ui.DQ_PlotWidget_NormIntensity.getAxis("bottom").setLabel(result["axis_label"])
        self.t2_dq_graph(show_warning=False)
        self.ui.DQ_PlotWidget_NormIntensity.plot(result["x_fit"], result["y_fit"], pen="r")
        self.ui.DQ_TextEdit_FitResult.setText(
            f"R²: {round(result['r_squared'], 4)} \n"
            f"X₀: {round(result['center'], 4)} \n"
            f"FWHM: {round(result['fwhm'], 4)} \n"
            f"Fraction (Lorenz): {round(result['lorenz_fraction'], 2)}"
        )

    def highlight_selected_point_widget_1(self):
        row = self.ui.DQ_Table_Data.currentRow()
        if row < 0:
            return

        x_item = self.ui.DQ_Table_Data.item(row, 0)
        y_item = self.ui.DQ_Table_Data.item(row, 3)
        if x_item is None or y_item is None:
            return

        try:
            x = float(x_item.text())
            y = float(y_item.text())
        except ValueError:
            return

        self.ui.DQ_PlotWidget_T2.plot(
            [x],
            [y],
            pen=None,
            symbol="o",
            symbolBrush=(255, 255, 0, 255),
            symbolPen="k",
            symbolSize=13,
        )

    def highlight_selected_point_widget_2(self):
        row = self.ui.DQ_Table_Data.currentRow()
        if row < 0:
            return

        x_item = self.ui.DQ_Table_Data.item(row, 4)
        y_item = self.ui.DQ_Table_Data.item(row, 5)
        if x_item is None or y_item is None:
            return

        try:
            x = float(x_item.text())
            y = float(y_item.text())
        except ValueError:
            return

        if self.ui.DQ_RadioButton_LogScale.isChecked() and x > 0:
            x = np.log10(x)

        self.ui.DQ_PlotWidget_NormIntensity.plot(
            [x],
            [y],
            pen=None,
            symbol="o",
            symbolBrush=(255, 255, 0, 255),
            symbolPen="k",
            symbolSize=13,
        )

    def _clear_plots(self):
        self.ui.DQ_PlotWidget_T2.clear()
        self.ui.DQ_PlotWidget_NormIntensity.clear()
        self.ui.DQ_TextEdit_FitResult.setText("")

    def _has_table_data(self, message, show_warning=True):
        if self.ui.DQ_Table_Data.rowCount() > 0:
            return True

        if show_warning:
            self._warn_no_data(message)
        return False

    def _has_linearized_columns(self, show_warning=False):
        has_columns = self.ui.DQ_Table_Data.columnCount() >= 6
        if has_columns:
            return True

        if show_warning:
            QMessageBox.warning(
                self.parent,
                "No DQ fit data",
                "Cannot calculate DQ fit because required columns are missing.",
                QMessageBox.Ok,
            )
        return False

    def _read_numeric_column(self, column_index, show_warning=False):
        try:
            return self.read_column_values(self.ui.DQ_Table_Data, column_index)
        except (ValueError, TypeError, IndexError):
            if show_warning:
                QMessageBox.warning(
                    self.parent,
                    "Invalid DQ data",
                    "Cannot calculate DQ fit because required columns contain non-numeric values.",
                    QMessageBox.Ok,
                )
            return None

    def _warn_no_data(self, message):
        QMessageBox.warning(self.parent, "No DQ data", message, QMessageBox.Ok)
