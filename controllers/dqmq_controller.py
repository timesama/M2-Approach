import numpy as np
import logging
import os
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem, QFileDialog
from pyqtgraph import InfiniteLine, mkPen
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

import Calculator as Cal
from controllers.base_tab_controller import BaseTabController
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


class DQMQTabController(BaseTabController):
    DRES_K = 0.4
    DRES_GRID = np.linspace(0, 0.94, 60000)

    def __init__(self, ui, state, parent=None):
        super().__init__(ui, state, parent)
        self.integral_sum_result = None
        self.dres_result = None

    def _get_current_file(self):
        files = self.state.dqmq_files or []
        if not files:
            return None
        return files[0]

    def dq_mq_analysis(self):
        with busy_cursor():
            table = self.ui.table_DQMQ
            logger.info("DQMQ analysis started")
        table.clear()
        table.setColumnCount(4)
        try:
            time, dq, ref = self.plot_original()
        except Exception as e:
            logger.exception("DQMQ analysis failed")
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
        logger.info("DQMQ analysis completed: %d points", len(time))

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
        logger.info("DQMQ normalization started")
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
        logger.info("DQMQ diff fit: from=%s to=%s power=%s noise=%s", fit_from, fit_to, p, noise_level)
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
        logger.info("DQMQ nDQ fit: from=%s to=%s power=%s noise=%s", fit_from, fit_to, p, noise_level)
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

    def calculate_integral_sum(self):
        file_path, _ = QFileDialog.getOpenFileName(self.parent, "Select DQ/T2 distribution CSV", "", "CSV Files (*.csv)")
        if not file_path:
            return
        data = self._read_dqt2_distribution_file(file_path)
        if data is None:
            QMessageBox.warning(self.parent, "Invalid file", "Selected file must contain exactly 6 columns.")
            return
        t, signal_norm = self._calculate_t2_summed_signal(data)
        shift = float(self.ui.DQMQSmooth_window_2.value())
        self._plot_integral_sum(t, signal_norm, shift)
        self.integral_sum_result = {
            "time": t + shift,
            "signal_norm": signal_norm,
            "source_file": file_path,
            "shift": shift,
        }

    def _read_dqt2_distribution_file(self, file_path):
        data = np.genfromtxt(file_path, delimiter=",")
        if data.ndim == 1:
            data = np.atleast_2d(data)
        if data.shape[1] != 6:
            return None
        return data

    def _calculate_t2_summed_signal(self, data, n_points=500):
        t_dq = data[:, 0]
        amp_dq = data[:, 1]
        t2 = data[:, 4]
        valid = np.isfinite(amp_dq) & np.isfinite(t2) & np.isfinite(t_dq) & (t2 > 0)
        amp_dq = amp_dq[valid]
        t2 = t2[valid]
        t_dq = t_dq[valid]
        order = np.argsort(t2)
        amp_dq = amp_dq[order]
        t2 = t2[order]
        t = np.linspace(0, np.max(t_dq), n_points)
        d_t2 = np.gradient(t2) if len(t2) > 1 else np.ones_like(t2)
        signal = np.zeros_like(t)
        for amp, t2i, delta_t2 in zip(amp_dq, t2, d_t2):
            signal += amp * np.exp(-(t / t2i) ** 2) * delta_t2
        s_min, s_max = np.min(signal), np.max(signal)
        signal_norm = np.zeros_like(signal) if s_max == s_min else (signal - s_min) / (s_max - s_min)
        return t, signal_norm

    def _plot_integral_sum(self, t, signal_norm, shift):
        self.ui.DQMQ_Widget.plot(
            t + shift,
            signal_norm,
            pen=mkPen((0, 100, 0), width=3),
            name=f"Integral sum, shift={shift:g}",
        )

    def save_integral_sum_result(self, base_file_path):
        if not self.integral_sum_result:
            return
        root, ext = os.path.splitext(base_file_path)
        save_path = f"{root}_IntegralSum{ext or '.csv'}"
        out = np.column_stack((self.integral_sum_result["time"], self.integral_sum_result["signal_norm"]))
        np.savetxt(save_path, out, delimiter=",", header="time,integral_sum_norm", comments="")

    def calculate_dres_distribution(self):
        file_path = self._get_current_file()
        if not file_path:
            QMessageBox.warning(self.parent, "No data", "Load a DQMQ file first.")
            return
        try:
            data = np.genfromtxt(file_path, delimiter=",")
            if data.ndim == 1:
                data = np.atleast_2d(data)
            tau = data[:, 0].astype(float)
            dq = data[:, 1].astype(float)
            ref = data[:, 2].astype(float)
            ndq = data[:, 3].astype(float)
            ref_max = np.max(ref)
            ndq_smoothed = savgol_filter(ndq, 3, 1) if len(ndq) >= 3 else ndq
            time0 = np.insert(tau + 1, 0, 0.0)
            ndq0 = np.insert(ndq_smoothed, 0, 0.0)
            kernel_map = {"Gauss": "gaussian", "Abragam": "abragam", "Pake": "pake", "Weibull": "weibull", "A-L": "a-l"}
            kernel = kernel_map.get(self.ui.DQMQ_Kernel_comboBox.currentText(), "gaussian")
            n_components = 2 if self.ui.radioButton_2.isChecked() else 1
            fit = self._fit_dres(time0, ndq0, kernel, n_components)
            d_plot, p_dist = self._distribution_from_fit(fit["popt"], n_components)
            fig = self.ui.DQMQ_Widget_DRes
            fig.clear()
            fig.plot(d_plot, p_dist, pen=mkPen('m', width=3), name=f"{kernel} {n_components}D")
            self.dres_result = {"dres_khz": d_plot, "p_dres": p_dist, "fit": fit}
        except Exception as exc:
            logger.exception("Dres calculation failed")
            QMessageBox.warning(self.parent, "Dres calculation failed", str(exc))

    def _fit_dres(self, time0, ndq0, kernel, n_components):
        model = self._make_fit_model(kernel, n_components)
        if n_components == 1:
            p0, lb, ub = [0.25, 1e-3, 2.0], [1e-7, 1e-7, 1e-7], [1.0, 0.8, 6.0]
        else:
            p0, lb, ub = [0.061, 0.05, 0.313, 0.033, 0.359, 0.96], [1e-3, 1e-4, 1e-5, 1e-15, 0.0, 0.1], [1.0, 1.5, 1.9, 1.1, 1.0, 4.0]
        popt, pcov = curve_fit(model, time0, ndq0, p0=p0, bounds=(lb, ub), maxfev=5000000)
        return {"popt": popt, "pcov": pcov}

    def _make_fit_model(self, kernel, n_components):
        return (lambda t, mu, sigma, beta: self._ndq_1d(t, mu, sigma, beta, kernel)) if n_components == 1 else (
            lambda t, mu1, sigma1, mu2, sigma2, frac1, beta: self._ndq_2d(t, mu1, sigma1, mu2, sigma2, frac1, beta, kernel)
        )

    def _dq_kernel(self, x, kernel, beta):
        beta = max(float(beta), 1e-12)
        if kernel == "gaussian":
            return 1.0 - np.exp(-self.DRES_K * x**2)
        if kernel == "abragam":
            return 1.0 - np.exp(-self.DRES_K * x**2) * np.sinc(x / np.pi)
        if kernel == "pake":
            return 1.0 - np.exp(-self.DRES_K * x**2) * np.cos(x)
        if kernel == "weibull":
            return 1.0 - np.exp(-self.DRES_K * x**beta)
        return 1.0 - np.exp(-self.DRES_K * x**beta) * np.cos(x)

    def _ndq_1d(self, t, mu, sigma, beta, kernel):
        p = self._p_gaussian(self.DRES_GRID, mu, sigma)
        return np.array([0.5 * np.trapz(p * self._dq_kernel(self.DRES_GRID * ti, kernel, beta), self.DRES_GRID) for ti in np.asarray(t)])

    def _ndq_2d(self, t, mu1, sigma1, mu2, sigma2, frac1, beta, kernel):
        p = self._p_gaussian_2d(self.DRES_GRID, mu1, sigma1, mu2, sigma2, frac1)
        return np.array([0.5 * np.trapz(p * self._dq_kernel(self.DRES_GRID * ti, kernel, beta), self.DRES_GRID) for ti in np.asarray(t)])

    def _p_gaussian(self, d, mu, sigma):
        sigma = max(float(sigma), 1e-9)
        p = np.exp(-(d - mu) ** 2 / (2 * sigma**2))
        area = np.trapz(p, d)
        return np.zeros_like(p) if area <= 0 else p / area

    def _p_gaussian_2d(self, d, mu1, sigma1, mu2, sigma2, frac1):
        frac1 = np.clip(frac1, 0.0, 1.0)
        return frac1 * self._p_gaussian(d, mu1, sigma1) + (1.0 - frac1) * self._p_gaussian(d, mu2, sigma2)

    def _distribution_from_fit(self, popt, n_components):
        p = self._p_gaussian(self.DRES_GRID, popt[0], popt[1]) if n_components == 1 else self._p_gaussian_2d(
            self.DRES_GRID, popt[0], popt[1], popt[2], popt[3], popt[4]
        )
        return self.DRES_GRID / (2 * np.pi) * 1000.0, p
