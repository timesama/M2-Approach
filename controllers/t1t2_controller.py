import os
import re
from itertools import islice

import numpy as np
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem

import Calculator as Cal
from controllers.base_tab_controller import BaseTabController
from dialogs.save_files_dialog import SaveFilesDialog


class T1T2TabController(BaseTabController):
    def update_T12_table(self):
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
        table = self.ui.table_T1
        combobox = self.ui.T1T2_ChooseFileComboBox
        pattern = r'(T1|T2)_(\s?-?\d+(\.\d+)?)((_.*)?\.(dat|txt))'
        dictionary = self.parent.tau_dictionary
        dictionary.clear()

        try:
            if os.path.splitext(selected_files[0])[1] == '.sef':
                if os.path.splitext(selected_files[1])[1] != '.sef':
                    QMessageBox.warning(self.parent, "Error", "Load two files: Magnetization and Profile in .sef format.", QMessageBox.Ok)
                    return
                if os.stat(selected_files[0]).st_size > os.stat(selected_files[1]).st_size:
                    filetoreadzones, filetoreaddata = selected_files[1], selected_files[0]
                else:
                    filetoreadzones, filetoreaddata = selected_files[0], selected_files[1]
                dictionary = self.process_files(filetoreadzones, filetoreaddata)
                table.setRowCount(len(dictionary))
                for row, key in zip(range(table.rowCount()), dictionary):
                    x_axis = dictionary[key]["X Axis"]
                    table.setItem(row, 0, QTableWidgetItem(filetoreaddata))
                    table.setItem(row, 1, QTableWidgetItem(os.path.basename(filetoreaddata)))
                    table.setItem(row, 2, QTableWidgetItem(x_axis))
                    combobox.addItem(f"{x_axis}")
                self.parent.tau_dictionary = dictionary
                self.parent.state_bad_code = True

            elif os.path.splitext(selected_files[0])[1] == '.csv':
                self.parent.state_bad_code = False
                table.setRowCount(len(selected_files))
                csv_pattern = r'_(\-?\d+)(?=\.csv$)'
                for row, file in zip(range(table.rowCount()), selected_files):
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
                    dictionary[file]["X Axis"].append(x_axis)
                    dictionary[file]["Time"].extend(Time)
                    dictionary[file]["Signal"].extend(Signal)
                    table.setItem(row, 0, QTableWidgetItem(file))
                    table.setItem(row, 1, QTableWidgetItem(current_file))
                    table.setItem(row, 2, QTableWidgetItem(str(x_axis)))
                    combobox.addItem(f"{current_file}")

            elif os.path.splitext(selected_files[0])[1] == '.txt':
                self.parent.state_bad_code = False
                table.setRowCount(len(selected_files * 4))
                pattern_all = r'T1_.*_(\s?-?\d+).txt'
                for file in selected_files:
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
                    x_axis = dictionary[entry]["X Axis"][0]
                    table.setItem(row, 0, QTableWidgetItem(entry))
                    table.setItem(row, 1, QTableWidgetItem(current_file))
                    table.setItem(row, 2, QTableWidgetItem(str(x_axis)))
                    combobox.addItem(f"{current_file}")
            else:
                self.parent.state_bad_code = False
                table.setRowCount(len(selected_files))
                for row, file in zip(range(table.rowCount()), selected_files):
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
                    table.setItem(row, 0, QTableWidgetItem(file))
                    table.setItem(row, 1, QTableWidgetItem(current_file))
                    table.setItem(row, 2, QTableWidgetItem(str(x_axis)))
                    dictionary[file]["X Axis"].append(x_axis)
                    dictionary[file]["Time"].extend(Time)
                    dictionary[file]["Signal"].extend(Signal)
                    combobox.addItem(f"{current_file}")
        except Exception:
            QMessageBox.warning(self.parent, "Error", "Something went wrong. Try again.", QMessageBox.Ok)
            self.parent.clear_list()

        self.ui.btn_Plot1.setEnabled(True)
        combobox.setCurrentIndex(-1)

    def process_files(self, file1, file2):
        zone_to_frequency = {}
        with open(file1, 'r') as f1:
            for line in islice(f1, 4, None):
                parts = line.split()
                zone_to_frequency[parts[-2]] = parts[0]
        dictionary = {}
        current_zone = None
        with open(file2, 'r') as f2:
            for line in f2:
                line = line.strip()
                if line.startswith("Zone"):
                    current_zone = line.split()[1]
                    continue
                if not line or line.startswith("TAU") or line.startswith("---"):
                    continue
                if current_zone in zone_to_frequency:
                    try:
                        parts = line.split()
                        time = float(parts[0]); signal = float(parts[1]); frequency = zone_to_frequency[current_zone]
                        name = file2 + ' at freq:' + str(frequency)
                        if name not in dictionary:
                            dictionary[name] = {"X Axis": frequency, "Time": [], "Signal": []}
                        dictionary[name]["Time"].append(time)
                        dictionary[name]["Signal"].append(signal)
                    except (ValueError, IndexError):
                        continue
        return dictionary

    def change_exponential_order(self):
        self.ui.DSB_ExpFitting1.setValue(100)
        self.ui.DSB_ExpFitting2.setValue(1000)
        self.ui.DSB_ExpFitting3.setValue(10000)
        self.ui.DSB_ExpFitting1.setEnabled(True)
        self.ui.DSB_ExpFitting2.setEnabled(not self.ui.T1T2_FitWith1ExpButton.isChecked())
        self.ui.DSB_ExpFitting3.setEnabled(self.ui.T1T2_FitWith3ExpButton.isChecked())
        self.calculate_relaxation_time()

    def calculate_relaxation_time(self):
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
            key = table.item(idx, 0).text() + ' at freq:' + table.item(idx, 2).text()
            denominator = 0.001
        else:
            key = table.item(idx, 0).text()
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
        try:
            tf, fit, tau1, tau2, tau3, r2, a1, a2, a3 = Cal.fit_exponent(t, s, order, p)
        except Exception:
            QMessageBox.warning(self.parent, "No covariance", f"I am sorry, I couldn't fit with {order} exponents. Decrease the order of fitting and try again.", QMessageBox.Ok)
            return
        self.ui.textEdit_error.setText(f"R² {r2}")
        self.ui.DSB_ExpFitting1.setValue(tau1); self.ui.DSB_ExpFitting2.setValue(tau2); self.ui.DSB_ExpFitting3.setValue(tau3)
        table.setItem(idx, 3, QTableWidgetItem(str(tau1))); table.setItem(idx, 5, QTableWidgetItem(str(tau2))); table.setItem(idx, 7, QTableWidgetItem(str(tau3)))
        table.setItem(idx, 4, QTableWidgetItem(str(a1))); table.setItem(idx, 6, QTableWidgetItem(str(a2))); table.setItem(idx, 8, QTableWidgetItem(str(a3)))
        figure.clear(); figure.plot(t0, s0, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10); figure.plot(tf, fit, pen='b')
        dictionary[key]['T1 1'] = tau1; dictionary[key]['T1 2'] = tau2; dictionary[key]['T1 3'] = tau3
        self.ui.btn_Plot1.setEnabled(True)

    def plot_relaxation_time(self):
        table = self.ui.table_T1
        graph = self.ui.T1_Widget_2
        column = 3 if self.ui.T1T2_1expPlotButton.isChecked() else 5 if self.ui.T1T2_2expPlotButton.isChecked() else 7
        graph.clear()
        if table.rowCount() < 1:
            return
        x_axis, relaxation_time, number = [], [], 1
        for row in range(table.rowCount()):
            try: x_axis.append(float(table.item(row, 2).text()))
            except: x_axis.append(number)
            try: relaxation_time.append(float(table.item(row, column).text()))
            except: relaxation_time.append(0)
            number += 1
        graph.plot(x_axis, relaxation_time, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)

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
