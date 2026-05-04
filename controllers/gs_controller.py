import os
import re

import numpy as np
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem

import Calculator as Cal
from controllers.base_tab_controller import BaseTabController


class GSTabController(BaseTabController):
    def update_GS_table(self):
        def clean_line(line):
            while '\t\t' in line:
                line = line.replace('\t\t', '\t')
            return line.strip()

        selected_files = self.parent.selected_GSfiles
        table = self.ui.table_GS
        combobox = self.ui.comboBox_7
        pattern = r'_([0-9]+)\.dat$'
        dictionary = self.parent.GS_dictionary
        dictionary.clear()

        try:
            table.setRowCount(len(selected_files))
            for row, file in zip(range(table.rowCount()), selected_files):
                dictionary[file] = {"X Axis": [], "sqrtTime": [], "short": [], "medium": [], "long": []}
                sqrtTime, short, medium, long = [], [], [], []
                current_file = os.path.basename(file)
                try:
                    x_axis = re.search(pattern, current_file).group(1)
                except Exception:
                    x_axis = row
                with open(file, "r") as data:
                    lines = [clean_line(line.rstrip('\n')) for line in data if line.strip()]
                for line in lines[1:]:
                    parts = line.split('\t')
                    sqrtTime.append(float(parts[0]))
                    short.append(float(parts[4]))
                    medium.append(float(parts[5]))
                    long.append(float(parts[6]))
                combobox.addItem(f"{current_file}")
                table.setItem(row, 0, QTableWidgetItem(file))
                table.setItem(row, 1, QTableWidgetItem(current_file))
                table.setItem(row, 2, QTableWidgetItem(str(x_axis)))
                dictionary[file]["X Axis"].append(x_axis)
                dictionary[file]["sqrtTime"].extend(sqrtTime)
                dictionary[file]["short"].extend(short)
                dictionary[file]["medium"].extend(medium)
                dictionary[file]["long"].extend(long)
        except Exception as e:
            QMessageBox.warning(self.parent, "Error", f"Something {e} went wrong. Try again.", QMessageBox.Ok)
            self.parent.clear_list()

    def calculate_sqrt_time(self):
        idx = self.ui.comboBox_7.currentIndex()
        if idx == -1:
            return
        table = self.ui.table_GS
        figure = self.ui.GS_Widget_1
        dictionary = self.parent.GS_dictionary
        key = table.item(idx, 0).text()
        time_original = np.array(dictionary[key]['sqrtTime']).flatten()
        if self.ui.checkBox_3.isChecked():
            time_original = np.sqrt(time_original)
        short_original = np.array(dictionary[key]['short']).flatten()
        medium_original = np.array(dictionary[key]['medium']).flatten()
        long_original = np.array(dictionary[key]['long']).flatten()
        self.ui.GS_fit_from_1.setMinimum(time_original[0]); self.ui.GS_fit_from_1.setMaximum(time_original[-1])
        self.ui.GS_fit_to_1.setMinimum(time_original[15]); self.ui.GS_fit_to_1.setMaximum(time_original[-1])
        from_val, to_val = self.ui.GS_fit_from_1.value(), self.ui.GS_fit_to_1.value()
        sp = (np.abs(time_original - from_val)).argmin(); ep = (np.abs(time_original - to_val)).argmin() + 1
        time = time_original[sp:ep]
        short = short_original[sp:ep]; medium = medium_original[sp:ep]; long = long_original[sp:ep]
        if self.ui.radioButton_short.isChecked(): signal, signal_original = short, short_original
        elif self.ui.radioButton_medium.isChecked(): signal, signal_original = medium, medium_original
        else: signal, signal_original = long, long_original
        try:
            tf, fit, sqrtT, r2v = Cal.linear_fit_GS(time, signal)
            d = Cal.calculate_domain_size(sqrtT, self.ui.GS_beta.value(), self.ui.GS_r2.value(), self.ui.GS_m2.value())
            self.ui.textEdit_error_2.setText(f"R² {r2v}")
            table.setItem(idx, 3, QTableWidgetItem(str(sqrtT)))
            table.setItem(idx, 4, QTableWidgetItem(str(d)))
            figure.clear(); figure.plot(time_original, signal_original, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10); figure.plot(tf, fit, pen='b')
            dictionary[key]['sqrtT'] = sqrtT; dictionary[key]['d'] = d
            self.ui.btn_Plot1.setEnabled(True)
        except Exception as e:
            figure.clear(); QMessageBox.warning(self.parent, "Error", f"Something {e} went wrong. Try again.", QMessageBox.Ok)

    def plot_sqrt_time(self):
        table = self.ui.table_GS
        graph = self.ui.GS_Widget_2
        graph.clear()
        if table.rowCount() < 1:
            return
        x_axis, sqrtT, n = [], [], 1
        for row in range(table.rowCount()):
            try: x_axis.append(float(table.item(row, 2).text()))
            except: x_axis.append(n)
            try: sqrtT.append(float(table.item(row, 3).text()))
            except: sqrtT.append(0)
            n += 1
        graph.plot(x_axis, sqrtT, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)
