import os
import re
import logging

import numpy as np
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem

import Calculator as Cal
from controllers.base_tab_controller import BaseTabController
from dialogs.save_files_dialog import SaveFilesDialog
from controllers.table_columns import T1Columns
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


class T1T2TabController(BaseTabController):
    def update_T12_table(self):
        with busy_cursor():
            self._update_T12_table_impl()

    def _update_T12_table_impl(self):
        def clean_line(line):
            while '\t\t' in line:
                line = line.replace('\t\t', '\t')
            return line.strip()

        def create_dictionary(dictionary, file, addition, x_axis, Time, Signal):
            dictionary[file + addition]["X Axis"].append(x_axis)
            dictionary[file + addition]["Time"].extend(Time)
            dictionary[file + addition]["Signal"].extend(Signal)
            return dictionary

        selected_files = self.parent.selected_T1files
        logger.info("T1T2 loading %d files", len(selected_files))
        table = self.ui.table_T1
        combobox = self.ui.T1T2_ChooseFileComboBox
        pattern = r'(T1|T2)_(\s?-?\d+(\.\d+)?)((_.*)?\.(dat|txt))'
        dictionary = self.parent.tau_dictionary
        preserved_x_axis = {}
        for row in range(table.rowCount()):
            file_item = table.item(row, T1Columns.FOLDER)
            x_item = table.item(row, T1Columns.X_AXIS)
            if file_item is not None and x_item is not None:
                preserved_x_axis[file_item.text()] = x_item.text()
        dictionary.clear()

        try:
            if os.path.splitext(selected_files[0])[1] == '.sef':
                logger.warning("T1T2 unsupported format: .sef")
                QMessageBox.warning(self.parent, "Unsupported format .sef", "Reading of NMRD data is disabled for the versions > 0.2.2.", QMessageBox.Ok)
                return
            elif os.path.splitext(selected_files[0])[1] == '.csv':
                logger.info("T1T2 branch: CSV")
                self.parent.state_bad_code = False
                table.setRowCount(len(selected_files))
                csv_pattern = r'_(\-?\d+)(?=\.csv$)'
                for row, file in zip(range(table.rowCount()), selected_files):
                    logger.info("T1T2 parsing file %d/%d: %s", row + 1, len(selected_files), os.path.basename(file))
                    current_file = os.path.basename(file)
                    try:
                        x_axis = re.search(csv_pattern, current_file).group(1)
                    except Exception:
                        x_axis = row
                    dictionary[file] = {"X Axis": [], "Time": [], "Signal": []}
                    Time, Signal = [], []
                    with open(file) as f:
                        for line in f:
                            parts = line.strip().split(",")
                            Time.append(float(parts[0]))
                            Signal.append(float(parts[1]))
                    effective_x = preserved_x_axis.get(file, str(x_axis))
                    dictionary[file]["X Axis"].append(effective_x)
                    dictionary[file]["Time"].extend(Time)
                    dictionary[file]["Signal"].extend(Signal)
                    table.setItem(row, T1Columns.FOLDER, QTableWidgetItem(file))
                    table.setItem(row, T1Columns.FILE_NAME, QTableWidgetItem(current_file))
                    table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(effective_x)))
                    combobox.addItem(f"{current_file}")

            elif os.path.splitext(selected_files[0])[1] == '.txt':
                logger.info("T1T2 branch: TXT")
                self.parent.state_bad_code = False
                table.setRowCount(len(selected_files * 4))
                pattern_all = r'T1_.*_(\s?-?\d+).txt'
                for file in selected_files:
                    logger.info("T1T2 parsing file: %s", os.path.basename(file))
                    try:
                        x_axis = re.search(pattern_all, os.path.basename(file)).group(1)
                    except Exception:
                        x_axis = 'Variable'
                    Time, Signal_all, Signal_short, Signal_med, Signal_long = [], [], [], [], []
                    dictionary[file + '_all'] = {"X Axis": [], "Time": [], "Signal": []}
                    dictionary[file + '_short'] = {"X Axis": [], "Time": [], "Signal": []}
                    dictionary[file + '_medium'] = {"X Axis": [], "Time": [], "Signal": []}
                    dictionary[file + '_long'] = {"X Axis": [], "Time": [], "Signal": []}
                    with open(file, "r") as data:
                        lines = [clean_line(line.rstrip('\n')) for line in data if line.strip()]
                        for line in lines[1:]:
                            parts = line.split('\t')
                            Time.append(float(parts[0]))
                            Signal_all.append(float(parts[1]))
                            Signal_short.append(float(parts[2]))
                            Signal_med.append(float(parts[3]))
                            Signal_long.append(float(parts[4]))
                    dictionary = create_dictionary(dictionary, file, '_all', x_axis, Time, Signal_all)
                    dictionary = create_dictionary(dictionary, file, '_short', x_axis, Time, Signal_short)
                    dictionary = create_dictionary(dictionary, file, '_medium', x_axis, Time, Signal_med)
                    dictionary = create_dictionary(dictionary, file, '_long', x_axis, Time, Signal_long)

                for row, entry in zip(range(table.rowCount()), dictionary):
                    current_file = os.path.basename(entry)
                    x_axis = preserved_x_axis.get(entry, dictionary[entry]["X Axis"][0])
                    table.setItem(row, T1Columns.FOLDER, QTableWidgetItem(entry))
                    table.setItem(row, T1Columns.FILE_NAME, QTableWidgetItem(current_file))
                    table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(preserved_x_axis.get(entry, x_axis))))
                    combobox.addItem(f"{current_file}")
            else:
                logger.info("T1T2 branch: default")
                self.parent.state_bad_code = False
                table.setRowCount(len(selected_files))
                for row, file in zip(range(table.rowCount()), selected_files):
                    logger.info("T1T2 parsing file %d/%d: %s", row + 1, len(selected_files), os.path.basename(file))
                    dictionary[file] = {"X Axis": [], "Time": [], "Signal": []}
                    current_file = os.path.basename(file)
                    try:
                        x_axis = re.search(pattern, current_file).group(2)
                    except Exception:
                        x_axis = row
                    try:
                        with open(file, "r") as data:
                            lines = [clean_line(line.rstrip('\n')) for line in data if line.strip()]
                        Time, Signal = [], []
                        for line in lines[1:]:
                            parts = line.split('\t')
                            Time.append(float(parts[0]))
                            Signal.append(float(parts[1]))
                    except Exception:
                        data = np.loadtxt(file)
                        Time, Signal = data[:, 0], data[:, 1]
                    table.setItem(row, T1Columns.FOLDER, QTableWidgetItem(file))
                    table.setItem(row, T1Columns.FILE_NAME, QTableWidgetItem(current_file))
                    effective_x = preserved_x_axis.get(file, str(x_axis))
                    table.setItem(row, T1Columns.X_AXIS, QTableWidgetItem(str(effective_x)))
                    dictionary[file]["X Axis"].append(effective_x)
                    dictionary[file]["Time"].extend(Time)
                    dictionary[file]["Signal"].extend(Signal)
                    combobox.addItem(f"{current_file}")
        except Exception:
            logger.exception("T1T2 loading failed")
            QMessageBox.warning(self.parent, "Error", "Something went wrong. Try again.", QMessageBox.Ok)
            self.parent.clear_list()

        self.ui.btn_Plot1.setEnabled(True)
        combobox.setCurrentIndex(-1)
        table.resizeColumnsToContents()

    def change_exponential_order(self):
        self.ui.DSB_ExpFitting1.setValue(100)
        self.ui.DSB_ExpFitting2.setValue(1000)
        self.ui.DSB_ExpFitting3.setValue(10000)
        self.ui.DSB_ExpFitting1.setEnabled(True)
        self.ui.DSB_ExpFitting2.setEnabled(not self.ui.T1T2_FitWith1ExpButton.isChecked())
        self.ui.DSB_ExpFitting3.setEnabled(self.ui.T1T2_FitWith3ExpButton.isChecked())
        self.calculate_relaxation_time()

    def calculate_relaxation_time(self):
        with busy_cursor():
            self._calculate_relaxation_time_impl()

    def _calculate_relaxation_time_impl(self):
        table = self.ui.table_T1
        figure = self.ui.T1_Widget_1
        idx = self.ui.T1T2_ChooseFileComboBox.currentIndex()
        dictionary = self.parent.tau_dictionary
        start = int(self.ui.T1T2_fit_from.value())
        end = -(int(self.ui.T1T2_fit_to.value())) or None
        denominator = 1 if self.ui.radioButton_16.isChecked() else 1000
        if idx == -1:
            return
        if self.parent.state_bad_code:
            key = table.item(idx, T1Columns.FOLDER).text() + ' at freq:' + table.item(idx, T1Columns.X_AXIS).text()
            denominator = 0.001
        else:
            key = table.item(idx, T1Columns.FOLDER).text()
        t0 = np.array(dictionary[key]['Time']) / denominator
        s0 = np.array(dictionary[key]['Signal'])
        t, s = t0[start:end], s0[start:end]
        t1, t2, t3 = int(self.ui.DSB_ExpFitting1.value()), int(self.ui.DSB_ExpFitting2.value()), int(self.ui.DSB_ExpFitting3.value())
        if self.ui.T1T2_FitWith1ExpButton.isChecked():
            order, p = 1, [s[0], t1, 1]
        elif self.ui.T1T2_FitWith2ExpButton.isChecked():
            order, p = 2, [s[0], t1, s[0], t2, 1]
        else:
            order, p = 3, [s[0], t1, s[0], t2, s[0], t3, 1]
        logger.info("T1T2 fitting selected curve: index=%d, order=%d", idx, order)
        try:
            tf, fit, tau1, tau2, tau3, r2, a1, a2, a3 = Cal.fit_exponent(t, s, order, p)
        except Exception:
            logger.exception("T1T2 fit failed for index=%d order=%d", idx, order)
            QMessageBox.warning(self.parent, "Fitting failed", f"Fitting failed for a {order}-exponential model. Decrease the coherence order and/or adjust the initial tau values.", QMessageBox.Ok)
            return
        logger.info("T1T2 fit completed: tau1=%s tau2=%s tau3=%s r2=%s", tau1, tau2, tau3, r2)
        self.ui.textEdit_error.setText(f"R² {r2}")
        table.setItem(idx, T1Columns.TAU_1, QTableWidgetItem(str(tau1))); table.setItem(idx, T1Columns.TAU_2, QTableWidgetItem(str(tau2))); table.setItem(idx, T1Columns.TAU_3, QTableWidgetItem(str(tau3)))
        table.setItem(idx, T1Columns.A_1, QTableWidgetItem(str(a1))); table.setItem(idx, T1Columns.A_2, QTableWidgetItem(str(a2))); table.setItem(idx, T1Columns.A_3, QTableWidgetItem(str(a3)))
        table.resizeColumnsToContents()
        figure.clear(); figure.plot(t0, s0, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10); figure.plot(tf, fit, pen='b')
        dictionary[key]['T1 1'] = tau1; dictionary[key]['T1 2'] = tau2; dictionary[key]['T1 3'] = tau3
        self.ui.btn_Plot1.setEnabled(True)

    def plot_relaxation_time(self):
        table = self.ui.table_T1
        graph = self.ui.T1_Widget_2
        column = T1Columns.TAU_1 if self.ui.T1T2_1expPlotButton.isChecked() else T1Columns.TAU_2 if self.ui.T1T2_2expPlotButton.isChecked() else T1Columns.TAU_3
        graph.clear()
        if table.rowCount() < 1:
            return
        x_axis, relaxation_time, number = [], [], 1
        for row in range(table.rowCount()):
            try: x_axis.append(float(table.item(row, T1Columns.X_AXIS).text()))
            except: x_axis.append(number)
            try: relaxation_time.append(float(table.item(row, column).text()))
            except: relaxation_time.append(0)
            number += 1
        graph.plot(x_axis, relaxation_time, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)
        self.highlight_selected_relaxation_point()

    def bad_code_makes_more_bad_code(self):
        dictionary = self.parent.tau_dictionary
        dialog = SaveFilesDialog(self.parent)
        basename = os.path.basename(self.parent.selected_T1files[0])
        save = not all(data.get('T1 1', 0) == 0 for data in dictionary.values())
        dialog.save_file_in_sef(self.parent, dictionary, 'T1 1', 1, basename, save)
        save = not all(data.get('T1 2', 0) == 0 for data in dictionary.values())
        dialog.save_file_in_sef(self.parent, dictionary, 'T1 2', 2, basename, save)
        save = not all(data.get('T1 3', 0) == 0 for data in dictionary.values())
        dialog.save_file_in_sef(self.parent, dictionary, 'T1 3', 3, basename, save)

    def highlight_selected_relaxation_point(self):
        row = self.ui.table_T1.currentRow()
        if row < 0:
            return
        column = T1Columns.TAU_1 if self.ui.T1T2_1expPlotButton.isChecked() else T1Columns.TAU_2 if self.ui.T1T2_2expPlotButton.isChecked() else T1Columns.TAU_3
        x_item = self.ui.table_T1.item(row, T1Columns.X_AXIS)
        y_item = self.ui.table_T1.item(row, column)
        if x_item is None or y_item is None:
            return
        try:
            x = float(x_item.text())
            y = float(y_item.text())
        except ValueError:
            return
        self.ui.T1_Widget_2.plot([x], [y], pen=None, symbol='o', symbolBrush=(255, 255, 0, 255), symbolPen='k', symbolSize=13)
