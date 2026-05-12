import logging
import os
import re

import numpy as np
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem

from calculations import t1t2_signal
from controllers.base_tab_controller import BaseTabController
from controllers.table_columns import T1Columns
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


class T1T2TabController(BaseTabController):
    def connect_signals(self):
        self.ui.T1T2_Button_Plot.clicked.connect(self.plot_relaxation_time_from_user)
        self.ui.T1T2_Table_Results.horizontalHeader().sectionDoubleClicked.connect(
            lambda index: self.parent.renameSection(self.ui.T1T2_Table_Results, index)
        )
        self.ui.T1T2_Table_Results.itemSelectionChanged.connect(self.plot_relaxation_time)
        self.ui.T1T2_Button_FitOneExp.clicked.connect(self.change_exponential_order)
        self.ui.T1T2_Button_FitTwoExp.clicked.connect(self.change_exponential_order)
        self.ui.T1T2_Button_FitThreeExp.clicked.connect(self.change_exponential_order)
        self.ui.T1T2_RadioButton_Seconds.clicked.connect(self.calculate_relaxation_time_from_user)
        self.ui.T1T2_RadioButton_Milliseconds.clicked.connect(self.calculate_relaxation_time_from_user)
        self.ui.T1T2_DoubleSpinBox_FitFrom.editingFinished.connect(self.calculate_relaxation_time_from_user)
        self.ui.T1T2_DoubleSpinBox_FitTo.editingFinished.connect(self.calculate_relaxation_time_from_user)
        self.ui.T1T2_DoubleSpinBox_InitialTau1.editingFinished.connect(self.calculate_relaxation_time_from_user)
        self.ui.T1T2_DoubleSpinBox_InitialTau2.editingFinished.connect(self.calculate_relaxation_time_from_user)
        self.ui.T1T2_DoubleSpinBox_InitialTau3.editingFinished.connect(self.calculate_relaxation_time_from_user)
        self.ui.T1T2_ComboBox_ChooseFile.activated.connect(
            lambda *_args: self.calculate_relaxation_time_from_user()
        )

    def update_T12_table(self):
        if self.parent.selected_T1files:
            self._status(f"Loading {len(self.parent.selected_T1files)} T1/T2 file(s)...")
        with busy_cursor():
            self._update_T12_table_impl()

    def _update_T12_table_impl(self):
        selected_files = self.parent.selected_T1files
        if not selected_files:
            self._warn_no_data("No T1/T2 data available. Select T1/T2 files first.")
            return

        legacy_files = [path for path in selected_files if path.lower().endswith(".sef")]
        if legacy_files:
            self._status("Unsupported T1/T2 file skipped.")
            QMessageBox.warning(
                self.parent,
                "Unsupported T1/T2 file",
                "Legacy FFC .sef files are no longer supported.\n\n"
                "Please convert the data to the current supported table/text/CSV format and try again.",
                QMessageBox.Ok,
            )
            self.parent.selected_T1files = [path for path in selected_files if not path.lower().endswith(".sef")]
            return

        logger.info("T1T2 loading %d files", len(selected_files))
        table = self.ui.T1T2_Table_Results
        combobox = self.ui.T1T2_ComboBox_ChooseFile
        dictionary = self.parent.tau_dictionary
        preserved_x_axis = self._preserved_x_axis(table)
        preserved_rows = self._preserved_rows(table)
        dictionary.clear()
        combobox.clear()

        extension = os.path.splitext(selected_files[0])[1].lower()
        failed_files = []

        if extension == ".csv":
            failed_files = self._load_csv_files(table, combobox, dictionary, preserved_rows, selected_files)
        elif extension == ".txt":
            failed_files = self._load_txt_files(table, combobox, dictionary, preserved_rows, selected_files)
        else:
            failed_files = self._load_default_files(table, combobox, dictionary, preserved_rows, selected_files)

        if failed_files:
            self._warn_failed_files("Some T1/T2 files could not be read. They were skipped.", failed_files)

        combobox.setCurrentIndex(-1)
        table.resizeColumnsToContents()
        loaded_count = table.rowCount()
        if loaded_count:
            self._status(f"Loaded {loaded_count} T1/T2 row(s).")
        else:
            self._status("Could not load T1/T2 files.")

    def _load_csv_files(self, table, combobox, dictionary, preserved_rows, selected_files):
        logger.info("T1T2 branch: CSV")
        table.setRowCount(0)
        csv_pattern = r"_(\-?\d+)(?=\.csv$)"
        failed_files = []

        for file_path in selected_files:
            row = table.rowCount()
            current_file = os.path.basename(file_path)
            logger.info("T1T2 parsing file %d/%d: %s", row + 1, len(selected_files), current_file)

            try:
                match = re.search(csv_pattern, current_file)
                x_axis = match.group(1) if match else row
                time_values, signal_values = self._read_csv(file_path)
            except Exception:
                logger.exception("T1T2 CSV parsing failed: %s", file_path)
                failed_files.append(current_file)
                continue

            # effective_x = preserved_x_axis.get(file_path, str(x_axis))
            # dictionary[file_path] = {"X Axis": [effective_x], "Time": time_values, "Signal": signal_values}
            effective_row = preserved_rows.get(file_path, )
            self._add_table_row(table, combobox, row, file_path, current_file, effective_row)

        return failed_files

    def _load_txt_files(self, table, combobox, dictionary, preserved_x_axis, selected_files):
        logger.info("T1T2 branch: TXT")
        table.setRowCount(0)
        pattern_all = r"T1_.*_(\s?-?\d+).txt"
        failed_files = []

        for file_path in selected_files:
            current_file = os.path.basename(file_path)
            logger.info("T1T2 parsing file: %s", current_file)

            try:
                match = re.search(pattern_all, current_file)
                x_axis = match.group(1) if match else "Variable"
                time_values, all_signal, short_signal, medium_signal, long_signal = self._read_txt_multi_signal(file_path)
            except Exception:
                logger.exception("T1T2 TXT parsing failed: %s", file_path)
                failed_files.append(current_file)
                continue

            for suffix, signal_values in (
                ("_all", all_signal),
                ("_short", short_signal),
                ("_medium", medium_signal),
                ("_long", long_signal),
            ):
                key = file_path + suffix
                dictionary[key] = {"X Axis": [], "Time": [], "Signal": []}
                t1t2_signal.add_curve(dictionary, file_path, suffix, x_axis, time_values, signal_values)
                row = table.rowCount()
                current_entry = os.path.basename(key)
                effective_x = preserved_x_axis.get(key, dictionary[key]["X Axis"][0])
                self._add_table_row(table, combobox, row, key, current_entry, effective_x)

        return failed_files

    def _load_default_files(self, table, combobox, dictionary, preserved_x_axis, selected_files):
        logger.info("T1T2 branch: default")
        table.setRowCount(0)
        pattern = r"(T1|T2)_(\s?-?\d+(\.\d+)?)((_.*)?\.(dat|txt))"
        failed_files = []

        for file_path in selected_files:
            row = table.rowCount()
            current_file = os.path.basename(file_path)
            logger.info("T1T2 parsing file %d/%d: %s", row + 1, len(selected_files), current_file)

            try:
                match = re.search(pattern, current_file)
                x_axis = match.group(2) if match else row
                time_values, signal_values = self._read_default_signal(file_path)
            except Exception:
                logger.exception("T1T2 parsing failed: %s", file_path)
                failed_files.append(current_file)
                continue

            effective_x = preserved_x_axis.get(file_path, str(x_axis))
            dictionary[file_path] = {"X Axis": [effective_x], "Time": list(time_values), "Signal": list(signal_values)}
            self._add_table_row(table, combobox, row, file_path, current_file, effective_x)

        return failed_files

    def _read_csv(self, file_path):
        time_values = []
        signal_values = []
        with open(file_path) as file:
            for line in file:
                parts = line.strip().split(",")
                time_values.append(float(parts[0]))
                signal_values.append(float(parts[1]))

        return time_values, signal_values

    def _read_txt_multi_signal(self, file_path):
        time_values = []
        all_signal = []
        short_signal = []
        medium_signal = []
        long_signal = []

        with open(file_path, "r") as data:
            lines = [t1t2_signal.clean_tabbed_line(line.rstrip("\n")) for line in data if line.strip()]

        for line in lines[1:]:
            parts = line.split("\t")
            time_values.append(float(parts[0]))
            all_signal.append(float(parts[1]))
            short_signal.append(float(parts[2]))
            medium_signal.append(float(parts[3]))
            long_signal.append(float(parts[4]))

        return time_values, all_signal, short_signal, medium_signal, long_signal

    def _read_default_signal(self, file_path):
        try:
            with open(file_path, "r") as data:
                lines = [t1t2_signal.clean_tabbed_line(line.rstrip("\n")) for line in data if line.strip()]

            time_values = []
            signal_values = []
            for line in lines[1:]:
                parts = line.split("\t")
                time_values.append(float(parts[0]))
                signal_values.append(float(parts[1]))
            return time_values, signal_values
        except Exception:
            data = np.loadtxt(file_path)
            return data[:, 0], data[:, 1]

    def _add_table_row1(self, table, combobox, row, folder, file_name, x_axis):
        table.insertRow(row)
        table.setItem(row, T1Columns.FOLDER, QTableWidgetItem(folder))
        table.setItem(row, T1Columns.FILE_NAME, QTableWidgetItem(file_name))
        table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(x_axis)))
        combobox.addItem(file_name)

    def _add_table_row(self, table, combobox, row, folder, file_name, row_data = None):
        table.insertRow(row)
        table.setItem(row, T1Columns.FOLDER, QTableWidgetItem(folder))
        table.setItem(row, T1Columns.FILE_NAME, QTableWidgetItem(file_name))

        if row_data:
            table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(row_data[0])))
            table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(row_data[1])))
            table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(row_data[2])))
            table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(row_data[3])))
            table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(row_data[4])))
            table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(row_data[5])))
        combobox.addItem(file_name)

    def _preserved_x_axis(self, table):
        preserved_x_axis = {}
        for row in range(table.rowCount()):
            file_item = table.item(row, T1Columns.FOLDER)
            x_item = table.item(row, T1Columns.X_AXIS)
            if file_item is not None and x_item is not None:
                preserved_x_axis[file_item.text()] = x_item.text()

        self.preserve_row(table)

        return preserved_x_axis

    def _preserved_rows(self, table):
        preserved_rows = {}

        for row in range(table.rowCount()):
            file_item = table.item(row, T1Columns.FOLDER)
            if file_item is None:
                continue

            file_key = file_item.text().strip()
            if not file_key:
                continue

            row_data = {}

            for col in range(table.columnCount()):
                item = table.item(row, col)
                row_data[col] = item.text() if item is not None else ""

            preserved_rows[file_key] = row_data

        return preserved_rows

    def change_exponential_order(self):
        self.calculate_relaxation_time(show_warning=True)

    def calculate_relaxation_time_from_user(self):
        self.calculate_relaxation_time(show_warning=True)

    def calculate_relaxation_time(self, show_warning=False):
        with busy_cursor():
            self._calculate_relaxation_time_impl(show_warning=show_warning)

    def _calculate_relaxation_time_impl(self, show_warning=False):
        table = self.ui.T1T2_Table_Results
        figure = self.ui.T1T2_PlotWidget_RawSignal
        idx = self.ui.T1T2_ComboBox_ChooseFile.currentIndex()
        dictionary = self.parent.tau_dictionary
        if table.rowCount() == 0 or not dictionary:
            if show_warning:
                self._warn_no_data("No T1/T2 data available. Select T1/T2 files first.")
            return

        if idx == -1:
            if show_warning:
                self._status("No rows selected.")
                QMessageBox.warning(self.parent, "No T1/T2 row selected", "No T1/T2 row selected.", QMessageBox.Ok)
            return

        item = table.item(idx, T1Columns.FOLDER)
        if item is None or item.text() not in dictionary:
            if show_warning:
                self._status("Could not fit T1/T2 data: non-numeric values.")
                QMessageBox.warning(
                    self.parent,
                    "No T1/T2 data",
                    "Cannot calculate relaxation time because required columns contain non-numeric values.",
                    QMessageBox.Ok,
                )
            return

        start = int(self.ui.T1T2_DoubleSpinBox_FitFrom.value())
        end = -(int(self.ui.T1T2_DoubleSpinBox_FitTo.value())) or None
        denominator = 1 if self.ui.T1T2_RadioButton_Seconds.isChecked() else 1000
        key = item.text()
        time_values, signal_values, fit_time, fit_signal = t1t2_signal.fit_range(
            dictionary[key]["Time"],
            dictionary[key]["Signal"],
            start,
            end,
            denominator,
        )
        tau1 = int(self.ui.T1T2_DoubleSpinBox_InitialTau1.value())
        tau2 = int(self.ui.T1T2_DoubleSpinBox_InitialTau2.value())
        tau3 = int(self.ui.T1T2_DoubleSpinBox_InitialTau3.value())

        if self.ui.T1T2_Button_FitOneExp.isChecked():
            order = 1
        elif self.ui.T1T2_Button_FitTwoExp.isChecked():
            order = 2
        else:
            order = 3

        try:
            initial_params = t1t2_signal.initial_parameters(order, fit_signal, tau1, tau2, tau3)
        except ValueError:
            if show_warning:
                self._status("Could not fit T1/T2 data: invalid range.")
                QMessageBox.warning(
                    self.parent,
                    "Invalid T1/T2 range",
                    "Cannot fit relaxation data because the selected range is invalid.",
                    QMessageBox.Ok,
                )
            return

        logger.info("T1T2 fitting selected curve: index=%d, order=%d", idx, order)
        try:
            fit_time_curve, fit_signal_curve, tau_1, tau_2, tau_3, r2, amp_1, amp_2, amp_3 = t1t2_signal.fit_relaxation(
                fit_time,
                fit_signal,
                order,
                initial_params,
            )
        except Exception:
            logger.exception("T1T2 fit failed for index=%d order=%d", idx, order)
            self._status("Could not fit T1/T2 data.")
            QMessageBox.warning(
                self.parent,
                "Fitting failed",
                f"Fitting failed for a {order}-exponential model. Decrease the coherence order and/or adjust the initial tau values.",
                QMessageBox.Ok,
            )
            return

        logger.info("T1T2 fit completed: tau1=%s tau2=%s tau3=%s r2=%s", tau_1, tau_2, tau_3, r2)
        if show_warning:
            self._status("Fit completed.")
        self.ui.T1T2_TextEdit_FitResult.setText(f"R² {r2}")
        table.setItem(idx, T1Columns.TAU_1, QTableWidgetItem(str(tau_1)))
        table.setItem(idx, T1Columns.TAU_2, QTableWidgetItem(str(tau_2)))
        table.setItem(idx, T1Columns.TAU_3, QTableWidgetItem(str(tau_3)))
        table.setItem(idx, T1Columns.A_1, QTableWidgetItem(str(amp_1)))
        table.setItem(idx, T1Columns.A_2, QTableWidgetItem(str(amp_2)))
        table.setItem(idx, T1Columns.A_3, QTableWidgetItem(str(amp_3)))
        table.resizeColumnsToContents()
        figure.clear()
        figure.plot(time_values, signal_values, pen=None, symbolPen=None, symbol="o", symbolBrush="r", symbolSize=10)
        figure.plot(fit_time_curve, fit_signal_curve, pen="b")
        dictionary[key]["T1 1"] = tau_1
        dictionary[key]["T1 2"] = tau_2
        dictionary[key]["T1 3"] = tau_3

    def plot_relaxation_time_from_user(self):
        self.plot_relaxation_time(show_warning=True)

    def plot_relaxation_time(self, show_warning=False):
        table = self.ui.T1T2_Table_Results
        graph = self.ui.T1T2_PlotWidget_RelaxationTime
        graph.clear()
        if table.rowCount() < 1:
            if show_warning:
                self._warn_no_data("Cannot update T1/T2 plot because the table is empty.")
            return

        column = self._selected_relaxation_column()
        x_axis = []
        relaxation_time = []
        number = 1
        for row in range(table.rowCount()):
            x_item = table.item(row, T1Columns.X_AXIS)
            y_item = table.item(row, column)
            try:
                x_axis.append(float(x_item.text()))
            except Exception:
                x_axis.append(number)
            try:
                relaxation_time.append(float(y_item.text()))
            except Exception:
                relaxation_time.append(0)
            number += 1

        graph.plot(x_axis, relaxation_time, pen=None, symbolPen=None, symbol="o", symbolBrush="r", symbolSize=10)

        self._plot_group_overlays()
        self.highlight_selected_relaxation_point()
        if show_warning:
            self._status("Plot updated.")

    def _plot_group_overlays(self):

        group_data = getattr(self.parent, "group_data_T1T2", {})
        colors = getattr(self.parent, "tab10_colors", [])
        if not group_data:
            return

        column = self._selected_relaxation_column()

        for i, (_group_number, group_rows) in enumerate(group_data.items()):
            group_x = []
            group_y = []
            for row_data in group_rows:
                try:
                    group_x.append(float(row_data[0]))
                    group_y.append(float(row_data[column]))
                except (ValueError, IndexError):
                    continue

            sorted_points = sorted(zip(group_x, group_y), key=lambda point: point[0])
            if len(sorted_points) > 1:
                xs, ys = zip(*sorted_points)
                color = colors[i % len(colors)] if colors else "r"
                self.ui.T1T2_PlotWidget_RelaxationTime.plot(
                    xs,
                    ys,
                    pen={"color": color, "width": 2},
                    symbol="o",
                    symbolBrush=color,
                    symbolPen=None,
                    symbolSize=8,
                )

    def highlight_selected_relaxation_point(self):
        row = self.ui.T1T2_Table_Results.currentRow()
        if row < 0:
            return

        column = self._selected_relaxation_column()
        x_item = self.ui.T1T2_Table_Results.item(row, T1Columns.X_AXIS)
        y_item = self.ui.T1T2_Table_Results.item(row, column)
        if x_item is None or y_item is None:
            return

        try:
            x = float(x_item.text())
            y = float(y_item.text())
        except ValueError:
            return

        self.ui.T1T2_PlotWidget_RelaxationTime.plot(
            [x],
            [y],
            pen=None,
            symbol="o",
            symbolBrush=(255, 255, 0, 255),
            symbolPen="k",
            symbolSize=13,
        )

    def _selected_relaxation_column(self):
        if self.ui.T1T2_RadioButton_PlotTau1.isChecked():
            return T1Columns.TAU_1
        if self.ui.T1T2_RadioButton_PlotTau2.isChecked():
            return T1Columns.TAU_2
        return T1Columns.TAU_3

    def _warn_no_data(self, message):
        self._status(message)
        QMessageBox.warning(self.parent, "No T1/T2 data", message, QMessageBox.Ok)

    def _warn_failed_files(self, message, failed_files):
        preview = "\n".join(failed_files[:5])
        if len(failed_files) > 5:
            preview += f"\n...and {len(failed_files) - 5} more"
        self._status(message)
        QMessageBox.warning(self.parent, "T1/T2 load warning", f"{message}\n\n{preview}", QMessageBox.Ok)
