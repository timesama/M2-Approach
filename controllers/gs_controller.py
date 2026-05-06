import os
import re
import logging

import numpy as np
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem

import Calculator as Cal
from controllers.base_tab_controller import BaseTabController
from controllers.table_columns import GSColumns
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


class GSTabController(BaseTabController):
    def update_GS_table(self):
        with busy_cursor():
            self._update_GS_table_impl()

    def _update_GS_table_impl(self):
        def clean_line(line):
            while '\t\t' in line:
                line = line.replace('\t\t', '\t')
            return line.strip()

        selected_files = self.parent.selected_GSfiles
        logger.info("GS loading %d files", len(selected_files))
        table = self.ui.table_GS
        combobox = self.ui.comboBox_7
        pattern = r'_([0-9]+)\.dat$'
        dictionary = self.parent.GS_dictionary
        dictionary.clear()

        try:
            table.setRowCount(len(selected_files))
            for row, file in zip(range(table.rowCount()), selected_files):
                logger.info("GS processing file %d/%d: %s", row + 1, len(selected_files), os.path.basename(file))
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
                table.setItem(row, GSColumns.FOLDER, QTableWidgetItem(file))
                table.setItem(row, GSColumns.FILE_NAME, QTableWidgetItem(current_file))
                table.setItem(row, GSColumns.X_AXIS, QTableWidgetItem(str(x_axis)))
                dictionary[file]["X Axis"].append(x_axis)
                dictionary[file]["sqrtTime"].extend(sqrtTime)
                dictionary[file]["short"].extend(short)
                dictionary[file]["medium"].extend(medium)
                dictionary[file]["long"].extend(long)
        except Exception as e:
            logger.exception("GS table loading failed")
            QMessageBox.warning(self.parent, "Error", f"Something {e} went wrong. Try again.", QMessageBox.Ok)
            self.parent.clear_list()
        table.resizeColumnsToContents()

    def calculate_sqrt_time(self):
        with busy_cursor():
            self._calculate_sqrt_time_impl()

    def _calculate_sqrt_time_impl(self):
        idx = self.ui.comboBox_7.currentIndex()
        if idx == -1:
            return
        logger.info("GS fit started: index=%d", idx)
        table = self.ui.table_GS
        figure = self.ui.GS_Widget_1
        dictionary = self.parent.GS_dictionary
        key = table.item(idx, GSColumns.FOLDER).text()
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
            table.setItem(idx, GSColumns.SQRT_TIME, QTableWidgetItem(str(sqrtT)))
            table.setItem(idx, GSColumns.D_NM, QTableWidgetItem(str(d)))
            table.resizeColumnsToContents()
            figure.clear(); figure.plot(time_original, signal_original, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10); figure.plot(tf, fit, pen='b')
            dictionary[key]['sqrtT'] = sqrtT; dictionary[key]['d'] = d
            logger.info("GS fit completed: sqrtT=%s d=%s", sqrtT, d)
            self.ui.btn_Plot_GS.setEnabled(True)
        except Exception as e:
            logger.exception("GS fit failed: index=%d", idx)
            figure.clear(); QMessageBox.warning(self.parent, "Error", f"Something {e} went wrong. Try again.", QMessageBox.Ok)

    def plot_sqrt_time(self):
        table = self.ui.table_GS
        graph = self.ui.GS_Widget_2
        graph.clear()
        if table.rowCount() < 1:
            return
        x_axis, sqrtT, n = [], [], 1
        for row in range(table.rowCount()):
            try: x_axis.append(float(table.item(row, GSColumns.X_AXIS).text()))
            except: x_axis.append(n)
            try: sqrtT.append(float(table.item(row, GSColumns.SQRT_TIME).text()))
            except: sqrtT.append(0)
            n += 1
        graph.plot(x_axis, sqrtT, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)
