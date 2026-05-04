import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem
from pyqtgraph import InfiniteLine, mkPen

import Calculator as Cal
from controllers.base_tab_controller import BaseTabController


class DQMQTabController(BaseTabController):
    def _get_current_file(self):
        files = self.state.dqmq_files or []
        if not files:
            return None
        return files[0]

    def dq_mq_analysis(self):
        table = self.ui.table_DQMQ
        table.clear()
        try:
            time, dq, ref = self.plot_original()
        except Exception as e:
            QMessageBox.warning(self.parent, "Corrupt File", f"Couldn't read the file beacuse {e}", QMessageBox.Ok)
            if self.parent is not None:
                self.parent.clear_list()
            return

        table.setRowCount(len(time))
        for row in range(len(time)):
            table.setItem(row, 0, QTableWidgetItem(str(time[row])))
            table.setItem(row, 1, QTableWidgetItem(str(dq[row])))
            table.setItem(row, 2, QTableWidgetItem(str(ref[row])))

        table.resizeColumnsToContents()
        table.setHorizontalHeaderLabels(['Time', 'DQ', 'Ref', 'nDQ'])
        self.ui.pushButton_DQMQ_1.setEnabled(True)
        self.ui.pushButton_DQMQ_2.setEnabled(False)
        self.ui.pushButton_DQMQ_3.setEnabled(False)
        self.ui.pushButton_DQMQ_4.setEnabled(True)

    def plot_original(self):
        file_path = self._get_current_file()
        if file_path is None:
            return np.array([]), np.array([]), np.array([])
        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        figure.addLegend()
        try:
            time, dq, ref = Cal.read_data(file_path, 1)
        except Exception:
            time, dq, ref = Cal.read_data(file_path, 0)

        time = time + int(self.ui.DQMQtime_shift.value())
        figure.plot(time, dq, pen=mkPen('r', width=3), name='DQ')
        figure.plot(time, ref, pen=mkPen('b', width=3), name='Ref')
        return time, dq, ref

    def plot_norm(self):
        file_path = self._get_current_file()
        if file_path is None:
            return
        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        figure.addLegend()
        noise_level = self.ui.noise.value()
        time_shift = int(self.ui.DQMQtime_shift.value())
        time, dq_norm, ref_norm, _, _, _, _, _, _, _ = Cal.dqmq(file_path, 40, 100, 1, noise_level, time_shift)
        figure.plot(time, dq_norm, pen=mkPen('r', width=3), name='DQ')
        figure.plot(time, ref_norm, pen=mkPen('b', width=3), name='Ref')
        self.ui.pushButton_DQMQ_2.setEnabled(True)
        self.ui.pushButton_DQMQ_3.setEnabled(True)
        self.ui.dq_min_3.setEnabled(True)
        self.ui.dq_max_3.setEnabled(True)
        self.ui.power.setEnabled(True)

    def plot_diff(self):
        file_path = self._get_current_file()
        if file_path is None:
            return
        fit_from = self.ui.dq_min_3.value()
        fit_to = self.ui.dq_max_3.value()
        p = self.ui.power.value()
        noise_level = self.ui.noise.value()
        time_shift = int(self.ui.DQMQtime_shift.value())
        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        figure.addLegend()
        time, _, _, diff, dq_norm, ref_norm, _, _, fitted_curve, _ = Cal.dqmq(
            file_path, fit_from, fit_to, p, noise_level, time_shift
        )
        figure.plot(time, dq_norm, pen=mkPen('r', width=2), name='DQ')
        figure.plot(time, ref_norm, pen=mkPen('b', width=2), name='Ref')
        figure.plot(time, diff, pen=mkPen('k', width=3), name='Diff')
        figure.plot(time, fitted_curve, pen=mkPen('m', width=3), name='fitting')
        self.ui.pushButton_DQMQ_2.setEnabled(True)
        self.ui.pushButton_DQMQ_3.setEnabled(True)

    def plot_nDQ(self):
        file_path = self._get_current_file()
        if file_path is None:
            return
        fit_from = self.ui.dq_min_3.value()
        fit_to = self.ui.dq_max_3.value()
        p = self.ui.power.value()
        noise_level = self.ui.noise.value()
        time_shift = int(self.ui.DQMQtime_shift.value())
        smoothing = [self.ui.DQMQSmooth_from.value(), self.ui.DQMQSmooth_to.value(), int(self.ui.DQMQSmooth_window.value())]
        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        figure.addLegend()
        time, _, _, _, dq_normal, ref_normal, time0, n_dq, _, mq_normal = Cal.dqmq(
            file_path, fit_from, fit_to, p, noise_level, time_shift, smoothing
        )
        figure.addItem(InfiniteLine(pos=0.5, angle=0, pen=mkPen(color=(200, 200, 255), width=2, style=Qt.DashLine)))
        figure.plot(time, dq_normal, pen=mkPen('r', width=3), name='DQ')
        figure.plot(time, ref_normal, pen=mkPen('b', width=3), name='Ref')
        figure.plot(time, mq_normal, pen=mkPen('m', width=3), name='MQ')
        figure.plot(time0, n_dq, pen=mkPen('k', width=3), symbol='o', symbolPen='k', symbolSize=10, name='nDQ')
        for row in range(len(time)):
            self.ui.table_DQMQ.setItem(row, 3, QTableWidgetItem(str(round(n_dq[row + 1], 4))))
        self.ui.table_DQMQ.resizeColumnsToContents()

    def plot_nDQ_on_Load(self):
        table = self.ui.table_DQMQ
        table.resizeColumnsToContents()
        time, dq_normal, ref_normal, n_dq = [], [], [], []
        for row in range(table.rowCount()):
            time.append(float(table.item(row, 0).text()))
            dq_normal.append(float(table.item(row, 1).text()))
            ref_normal.append(float(table.item(row, 2).text()))
            n_dq.append(float(table.item(row, 3).text()))

        time = np.array(time)
        dq_normal = np.array(dq_normal)
        ref_normal = np.array(ref_normal)
        dq_normal, ref_normal, _ = Cal.normalize_mq(dq_normal, ref_normal, 'plus')
        n_dq = np.insert(np.array(n_dq), 0, 0)
        time0 = np.insert(time, 0, 0)

        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        figure.addLegend()
        figure.plot(time, dq_normal, pen=mkPen('r', width=3), name='DQ')
        figure.plot(time, ref_normal, pen=mkPen('b', width=3), name='Ref')
        figure.plot(time0, n_dq, pen=mkPen('k', width=3), symbol='o', symbolPen='k', symbolSize=10, name='nDQ')
