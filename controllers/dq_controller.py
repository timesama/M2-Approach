import numpy as np
from PySide6.QtWidgets import QTableWidgetItem
from scipy.optimize import curve_fit

import Calculator as Cal
from controllers.base_tab_controller import BaseTabController


class DQTabController(BaseTabController):
    def process_processed_file(self, i, filename, x, y, z, m2, t2, file_path):
        import re
        match = re.search(r'_(\d+\.\d+)_', filename)
        dq_time = match.group(1) if match else '0'
        amplitude = Cal._calculate_amplitude(y, z)
        dq = Cal.calculate_DQ_intensity(x, amplitude)
        self.ui.table_DQ.setRowCount(i)
        self.parent.fill_table(self.ui.table_DQ, dq_time, dq, m2, t2, i)
        if self.ui.radioButton.isChecked():
            self.parent.save_figures(file_path, dq_time)
        self.update_graphs()

    def update_graphs(self):
        files = getattr(getattr(self.state, "dq", None), "files", None)
        if files is None:
            files = getattr(self.state, "dq_files", [])
        if len(files) > 1:
            self.linearization()
            self.plot_fit()
        else:
            self.dq_t2_graph()

    def dq_t2_graph(self):
        x = self.read_column_values(self.ui.table_DQ, 0)
        y = self.read_column_values(self.ui.table_DQ, 3)
        self.ui.DQ_Widget_1.clear()
        self.ui.DQ_Widget_1.plot(
            x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=10
        )

    def linearization(self):
        time_min = self.ui.dq_min.value()
        time_max = self.ui.dq_max.value()
        dq_time = np.array(self.read_column_values(self.ui.table_DQ, 0))
        t2 = np.array(self.read_column_values(self.ui.table_DQ, 3))
        dq = np.array(self.read_column_values(self.ui.table_DQ, 1))

        x = dq_time[(dq_time >= time_min) & (dq_time <= time_max)]
        if len(x) < 3:
            time_min = 0
            time_max = 20
            self.ui.dq_min.setValue(time_min)
            self.ui.dq_max.setValue(time_max)
            x = dq_time[(dq_time >= time_min) & (dq_time <= time_max)]

        y = t2[(dq_time >= time_min) & (dq_time <= time_max)]
        if len(x) <= 1 or len(y) <= 1:
            return

        coeff = np.polyfit(x, y, 1)
        integral = np.trapz(dq)
        dq_norm = dq / integral

        self.ui.table_DQ.setColumnCount(5)
        self.ui.table_DQ.setColumnWidth(4, 70)
        self.ui.table_DQ.setHorizontalHeaderItem(4, QTableWidgetItem("T₂* lin"))
        self.ui.table_DQ.setColumnCount(6)
        self.ui.table_DQ.setColumnWidth(5, 70)
        self.ui.table_DQ.setHorizontalHeaderItem(5, QTableWidgetItem("DQ Norm"))

        for row in range(self.ui.table_DQ.rowCount()):
            t2_lin = round(coeff[0] * dq_time[row] + coeff[1], 4)
            self.ui.table_DQ.setItem(row, 4, QTableWidgetItem(str(t2_lin)))
            self.ui.table_DQ.setItem(row, 5, QTableWidgetItem(str(round(dq_norm[row], 4))))

        x_line = np.arange(0, 105.1, 0.1)
        y_line = np.polyval(coeff, x_line)
        self.dq_t2_graph()
        self.ui.DQ_Widget_1.plot(x_line, y_line, pen='r')
        self.t2_dq_graph()

    def t2_dq_graph(self):
        x = self.read_column_values(self.ui.table_DQ, 4)
        y = self.read_column_values(self.ui.table_DQ, 5)

        if self.ui.radioButton_Log.isChecked():
            new_x = np.log10(x)
            self.ui.DQ_Widget_2.getAxis('bottom').setLabel("log(T₂*)")
        else:
            new_x = x
            self.ui.DQ_Widget_2.getAxis('bottom').setLabel("T₂*")

        self.ui.DQ_Widget_2.clear()
        self.ui.DQ_Widget_2.plot(
            new_x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=10
        )

    def plot_fit(self):
        _x = self.read_column_values(self.ui.table_DQ, 4)
        y = self.read_column_values(self.ui.table_DQ, 5)

        button = self.ui.radioButton_Log
        if button.isChecked():
            x = np.log10(_x)
            p = [1, 1, 1, 0]
            b = ([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1, np.inf])
            self.ui.DQ_Widget_2.getAxis('bottom').setLabel("log(T₂*)")
        else:
            x = _x
            p = [1, 10, 10, 0]
            b = ([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1, np.inf])
            self.ui.DQ_Widget_2.getAxis('bottom').setLabel("T₂*")

        text = self.ui.comboBox_FunctionDQ.currentText()
        try:
            x_fit = np.arange(0, np.max(x) + 0.001, 0.01)
        except Exception:
            x_fit = np.arange(0, 100 + 0.001, 0.01)

        b1 = ([0, 0, 0, -10], [np.inf, np.inf, np.inf, np.inf])
        if text == 'Gauss':
            params, _ = curve_fit(Cal.gaussian, x, y, p0=p, bounds=b1)
            y_fit, y_r2, cen, fwhm, w = Cal.gaussian(x_fit, *params), Cal.gaussian(x, *params), params[1], params[2], 0
            button.setEnabled(True)
        elif text == 'Lorenz':
            params, _ = curve_fit(Cal.lorenz, x, y, p0=p, bounds=b1)
            y_fit, y_r2, cen, fwhm, w = Cal.lorenz(x_fit, *params), Cal.lorenz(x, *params), params[1], params[2], 1
            button.setEnabled(True)
        elif text == 'Pseudo Voigt':
            params, _ = curve_fit(Cal.voigt, x, y, bounds=b)
            y_fit, y_r2, cen, fwhm, w = Cal.voigt(x_fit, *params), Cal.voigt(x, *params), params[1], params[2], params[3]
            button.setEnabled(True)
        else:
            button.setEnabled(False)
            return

        r_squared = Cal.calculate_r_squared(y, y_r2)
        self.t2_dq_graph()
        self.ui.DQ_Widget_2.plot(x_fit, y_fit, pen='r')
        self.ui.textEdit_4.setText(
            f"R\u00B2: {round(r_squared, 4)} \nX\u2080: {round(cen, 4)} \nFWHM: {round(fwhm, 4)} \nFraction (Lorenz): {round(w,2)}"
        )
