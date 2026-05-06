import numpy as np
import logging
import os
import importlib.util
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem, QFileDialog
from pyqtgraph import InfiniteLine, mkPen
from scipy.signal import savgol_filter

import Calculator as Cal
from calculations import dqmq_dres
from controllers.base_tab_controller import BaseTabController
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


class DQMQTabController(BaseTabController):
    def __init__(self, ui, state, parent=None):
        super().__init__(ui, state, parent)
        self.integral_sum_result = None
        self.integral_sum_curve_item = None
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
        self.integral_sum_curve_item = None
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
        self.integral_sum_curve_item = None
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
        self.integral_sum_curve_item = None
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
        self.integral_sum_curve_item = None
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
        self.integral_sum_curve_item = None
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
        try:
            with busy_cursor():
                data = self._read_dqt2_distribution_file(file_path)
                t, signal_norm = self._calculate_t2_summed_signal(data)
                shift = float(self.ui.DQMQSmooth_window_2.value())
                self._plot_integral_sum(t, signal_norm, shift)
                self.integral_sum_result = {
                    "time": t + shift,
                    "time_base": t,
                    "signal_norm": signal_norm,
                    "source_file": file_path,
                    "shift": shift,
                }
        except ValueError as exc:
            QMessageBox.warning(self.parent, "Invalid file", str(exc))
        except Exception as exc:
            logger.exception("Integral sum calculation failed")
            QMessageBox.warning(self.parent, "Integral sum calculation failed", str(exc))

    def _read_dqt2_distribution_file(self, file_path):
        data = np.genfromtxt(file_path, delimiter=",")
        if data.ndim == 0:
            raise ValueError("Selected file must contain exactly 6 columns.")
        if data.ndim == 1:
            data = np.atleast_2d(data)
        if data.shape[1] != 6:
            raise ValueError("Selected file must contain exactly 6 columns.")
        return data

    def _calculate_t2_summed_signal(self, data, n_points=500):
        t_dq = data[:, 0]
        amp_dq = data[:, 1]
        t2 = data[:, 4]
        valid = np.isfinite(amp_dq) & np.isfinite(t2) & np.isfinite(t_dq) & (t2 > 0)
        amp_dq = amp_dq[valid]
        t2 = t2[valid]
        t_dq = t_dq[valid]
        if len(t2) == 0:
            raise ValueError("Selected file does not contain valid positive T2star_lin values.")
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
        shifted_time = t + shift
        if self.integral_sum_curve_item is None:
            self.integral_sum_curve_item = self.ui.DQMQ_Widget.plot(
                shifted_time,
                signal_norm,
                pen=mkPen((0, 100, 0), width=3),
                name=f"Integral sum, shift={shift:g}",
            )
        else:
            self.integral_sum_curve_item.setData(shifted_time, signal_norm)
            self.integral_sum_curve_item.setName(f"Integral sum, shift={shift:g}")

    def update_integral_sum_shift(self):
        if not self.integral_sum_result:
            return
        shift = float(self.ui.DQMQSmooth_window_2.value())
        self.integral_sum_result["shift"] = shift
        self.integral_sum_result["time"] = self.integral_sum_result["time_base"] + shift
        self._plot_integral_sum(
            self.integral_sum_result["time_base"],
            self.integral_sum_result["signal_norm"],
            shift,
        )

    def save_integral_sum_result(self, base_file_path):
        if not self.integral_sum_result:
            return
        root, ext = os.path.splitext(base_file_path)
        save_path = f"{root}_IntegralSum{ext or '.csv'}"
        out = np.column_stack((self.integral_sum_result["time"], self.integral_sum_result["signal_norm"]))
        np.savetxt(save_path, out, delimiter=",", header="time,integral_sum_norm", comments="")

    def calculate_dres(self):
        self.calculate_dres_distribution()

    def calculate_dres_distribution(self):
        try:
            with busy_cursor():
                arrays = self._read_dqmq_table_arrays()
                kernel = self._selected_dres_kernel()
                n_components = self._selected_dres_component_count()
                fit_result = dqmq_dres.fit_selected_model(
                    arrays["Time0"],
                    arrays["nDQ0"],
                    kernel=kernel,
                    n_components=n_components,
                )
                d_plot, p_dist = dqmq_dres.build_distribution(fit_result)
                self.dres_result = {
                    "D_plot": d_plot,
                    "P": p_dist,
                    "fit_x": arrays["Time0"],
                    "fit_y": fit_result["fit"],
                    "kernel": kernel,
                    "n_components": n_components,
                    "params": fit_result["popt"],
                    "param_names": fit_result["param_names"],
                }
                self._plot_dres_result(arrays, self.dres_result)
        except ValueError as exc:
            QMessageBox.warning(self.parent, "Dres calculation", str(exc))
        except Exception as exc:
            logger.exception("Dres calculation failed")
            QMessageBox.warning(self.parent, "Dres calculation failed", str(exc))

    def _read_dqmq_table_arrays(self):
        table = self.ui.table_DQMQ
        if table.rowCount() == 0:
            raise ValueError("Run DQMQ analysis first so the table contains nDQ values.")

        values = {"tau": [], "dq": [], "ref": [], "nDQ": []}
        column_names = [("tau", 0, "Time"), ("dq", 1, "DQ"), ("ref", 2, "Ref"), ("nDQ", 3, "nDQ")]
        for row in range(table.rowCount()):
            row_values = {}
            row_has_any_value = False
            for key, column, label in column_names:
                item = table.item(row, column)
                text = item.text().strip() if item is not None else ""
                row_has_any_value = row_has_any_value or bool(text)
                if not text:
                    row_values[key] = None
                    continue
                try:
                    row_values[key] = float(text)
                except ValueError as exc:
                    raise ValueError(f"Row {row + 1} has a non-numeric {label} value.") from exc

            if not row_has_any_value:
                continue
            if any(row_values[key] is None for key, _, _ in column_names):
                raise ValueError(f"Row {row + 1} is missing one or more DQMQ table values.")
            for key in values:
                values[key].append(row_values[key])

        if not values["nDQ"]:
            raise ValueError("The DQMQ table nDQ column is empty. Run Plot nDQ before calculating Dres.")
        if len(values["nDQ"]) < 3:
            raise ValueError("Need at least 3 DQMQ table rows with nDQ values to calculate Dres.")

        tau = np.asarray(values["tau"], dtype=float)
        dq = np.asarray(values["dq"], dtype=float)
        ref = np.asarray(values["ref"], dtype=float)
        ndq_original = np.asarray(values["nDQ"], dtype=float)
        ndq = savgol_filter(ndq_original, 3, 1) if len(ndq_original) >= 3 else ndq_original

        if np.isclose(tau[0], 0.0):
            time0 = tau
            ndq0 = ndq
        else:
            time0 = np.insert(tau, 0, 0.0)
            ndq0 = np.insert(ndq, 0, 0.0)

        return {
            "tau": tau,
            "DQ": dq,
            "Ref": ref,
            "nDQ": ndq,
            "nDQ_original": ndq_original,
            "Time0": time0,
            "nDQ0": ndq0,
        }

    def _selected_dres_kernel(self):
        kernel_map = {
            "Gauss": "gaussian",
            "Abragam": "abragam",
            "Pake": "pake",
            "Weibull": "weibull",
            "A-L": "a-l",
        }
        selected = self.ui.DQMQ_Kernel_comboBox.currentText()
        if selected not in kernel_map:
            raise ValueError(f"Unknown Dres kernel: {selected}")
        return kernel_map[selected]

    def _selected_dres_component_count(self):
        if self.ui.radioButton_3.isChecked():
            return 1
        if self.ui.radioButton_2.isChecked():
            return 2
        raise ValueError("Select either 1 Dres or 2 Dres components.")

    def _plot_dres_result(self, arrays, dres_result):
        dres_figure = self.ui.DQMQ_Widget_DRes
        dres_figure.clear()
        dres_figure.addLegend()
        dres_figure.plot(
            dres_result["D_plot"],
            dres_result["P"],
            pen=mkPen('m', width=3),
            name="Dres distribution",
        )

        figure = self.ui.DQMQ_Widget
        figure.clear()
        self.integral_sum_curve_item = None
        figure.addLegend()
        figure.plot(
            arrays["Time0"],
            arrays["nDQ0"],
            pen=mkPen('k', width=3),
            symbol='o',
            symbolPen='k',
            symbolSize=10,
            name="nDQ",
        )
        figure.plot(
            dres_result["fit_x"],
            dres_result["fit_y"],
            pen=mkPen('k', width=4, style=Qt.DashLine),
            name="Dres fit",
        )

    def save_dres_result(self, base_file_path):
        if not self.dres_result:
            return
        root, _ = os.path.splitext(base_file_path)
        distribution = np.column_stack((self.dres_result["D_plot"], self.dres_result["P"]))
        fit = np.column_stack((self.dres_result["fit_x"], self.dres_result["fit_y"]))
        metadata = {
            "kernel": self.dres_result["kernel"],
            "n_components": self.dres_result["n_components"],
            **dict(zip(self.dres_result["param_names"], self.dres_result["params"])),
        }

        can_write_excel = (
            importlib.util.find_spec("pandas") is not None
            and (
                importlib.util.find_spec("openpyxl") is not None
                or importlib.util.find_spec("xlsxwriter") is not None
            )
        )
        if can_write_excel:
            import pandas as pd
            try:
                with pd.ExcelWriter(f"{root}_Dres.xlsx") as writer:
                    pd.DataFrame(distribution, columns=["Dres_over_2pi_kHz", "P_Dres"]).to_excel(
                        writer,
                        sheet_name="Dres distribution",
                        index=False,
                    )
                    pd.DataFrame(fit, columns=["time", "fitted_nDQ"]).to_excel(
                        writer,
                        sheet_name="Fit",
                        index=False,
                    )
                    pd.DataFrame(metadata.items(), columns=["parameter", "value"]).to_excel(
                        writer,
                        sheet_name="Metadata",
                        index=False,
                    )
                return
            except Exception:
                logger.exception("Excel Dres export failed; writing CSV fallback files")

        np.savetxt(
            f"{root}_Dres_distribution.csv",
            distribution,
            delimiter=",",
            header="Dres_over_2pi_kHz,P_Dres",
            comments="",
        )
        metadata_header = "\n".join([f"{key}: {value}" for key, value in metadata.items()])
        np.savetxt(
            f"{root}_Dres_fit.csv",
            fit,
            delimiter=",",
            header=f"{metadata_header}\ntime,fitted_nDQ",
            comments="# ",
        )

