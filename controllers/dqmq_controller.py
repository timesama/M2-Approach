import importlib.util
import logging
import os

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFileDialog, QMessageBox, QTableWidgetItem
from pyqtgraph import InfiniteLine, mkPen
from scipy.signal import savgol_filter

import Calculator as Cal
from calculations import dqmq_dres, dqmq_signal
from controllers.base_tab_controller import BaseTabController
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)

DQMQ_PLOT_STYLES = {
    "dq": {"color": (255, 0, 0), "width": 3},                 # red
    "ndq": {"color": (220, 20, 60), "width": 3},              # crimson

    "ref": {"color": (30, 144, 255), "width": 3},             # dodger blue
    "diff": {"color": (0, 0, 139), "width": 3},               # dark blue

    "mq": {"color": (255, 165, 0), "width": 3},               # orange
    "mq_minus_tail": {"color": (190, 70, 0), "width": 3},     # brown

    "tail_fit": {"color": (128, 128, 128), "width": 3},       # grey
    "dres_fit": {"color": (255, 105, 180), "width": 4},       # light magenta / hot pink

    "integral_sum": {"color": (0, 255, 0), "width": 3},       # bright green
}


class DQMQTabController(BaseTabController):
    def __init__(self, ui, state, parent=None):
        super().__init__(ui, state, parent)
        self.raw_data = None
        self.analysis_result = None
        self.integral_sum_result = None
        self.integral_sum_curve_item = None
        self.dres_result = None
        self.dres_is_stale = False
        self.dres_auto_calculated = False
        self.current_plot_mode = "raw"
        self.apply_plot_style_to_checkboxes()

    def apply_plot_style_to_checkboxes(self):
        checkbox_styles = {
            "DQMQ_CheckBox_DQ": "dq",
            "DQMQ_CheckBox_NDQ": "ndq",
            "DQMQ_CheckBox_Reference": "ref",
            "DQMQ_CheckBox_Difference": "diff",
            "DQMQ_CheckBox_MQ": "mq",
            "DQMQ_CheckBox_MQTailDiff": "mq_minus_tail",
            "DQMQ_CheckBox_TailFitting": "tail_fit",
            "DQMQ_CheckBox_DresFitting": "dres_fit",
            "DQMQ_CheckBox_IntegralSum": "integral_sum",
        }

        for checkbox_name, style_name in checkbox_styles.items():
            checkbox = getattr(self.ui, checkbox_name, None)
            if checkbox is None:
                continue

            color = DQMQ_PLOT_STYLES[style_name]["color"]
            checkbox.setStyleSheet(f"color: rgb{color};")

    def reset_cached_results(self, clear_raw=True):
        if clear_raw:
            self.raw_data = None
        self.analysis_result = None
        self.integral_sum_result = None
        self.integral_sum_curve_item = None
        self.dres_result = None
        self.dres_is_stale = False
        self.dres_auto_calculated = False
        self.current_plot_mode = "raw"
        if clear_raw:
            self.ui.DQMQ_PlotWidget_Signal.clear()
        self.ui.DQMQ_PlotWidget_Dres.clear()

    def _get_current_file(self):
        files = self.state.dqmq_files or []
        if not files:
            return None
        return files[0]

    def _get_widget(self, *names):
        for name in names:
            widget = getattr(self.ui, name, None)
            if widget is not None:
                return widget

        logger.warning("Missing DQMQ UI widget. Tried: %s", ", ".join(names))
        return None

    def _is_checked(self, *names, default=True):
        widget = self._get_widget(*names)
        if widget is None:
            return default
        return widget.isChecked()

    def _clear_signal_plot(self):
        figure = self.ui.DQMQ_PlotWidget_Signal
        figure.clear()
        self.integral_sum_curve_item = None
        return figure

    def _read_raw_data_from_file(self, file_path):
        try:
            time, dq, ref = Cal.read_data(file_path, 1)
        except Exception:
            time, dq, ref = Cal.read_data(file_path, 0)

        self.raw_data = {
            "time": time,
            "dq_raw": dq,
            "ref_raw": ref,
            "source_file": file_path,
        }
        return self.raw_data

    def _ensure_raw_data(self):
        file_path = self._get_current_file()
        if file_path is None:
            self._status("Select a DQMQ file first.")
            QMessageBox.warning(
                self.parent,
                "DQMQ file missing",
                "Select a DQMQ file before plotting.",
                QMessageBox.Ok,
            )
            return None

        if self.raw_data and self.raw_data.get("source_file") == file_path:
            return self.raw_data

        return self._read_raw_data_from_file(file_path)

    def _write_raw_table(self, raw_data):
        table = self.ui.DQMQ_Table_Data
        time = raw_data["time"]
        dq = raw_data["dq_raw"]
        ref = raw_data["ref_raw"]

        table.clear()
        table.setColumnCount(4)
        table.setRowCount(len(time))
        table.setHorizontalHeaderLabels(["Time", "DQ", "Ref", "nDQ"])
        for row in range(len(time)):
            table.setItem(row, 0, QTableWidgetItem(str(time[row])))
            table.setItem(row, 1, QTableWidgetItem(str(dq[row])))
            table.setItem(row, 2, QTableWidgetItem(str(ref[row])))
            table.setItem(row, 3, QTableWidgetItem(""))

        table.resizeColumnsToContents()

    def _write_ndq_to_table(self, time, n_dq):
        table = self.ui.DQMQ_Table_Data
        required_rows = len(time)
        if table.rowCount() < required_rows:
            table.setRowCount(required_rows)

        for row in range(required_rows):
            n_dq_value = n_dq[row + 1]
            if np.isfinite(n_dq_value):
                table_text = str(round(n_dq_value, 4))
            else:
                table_text = "NaN"
            table.setItem(row, 3, QTableWidgetItem(table_text))

        table.resizeColumnsToContents()

    def dq_mq_analysis(self):
        with busy_cursor():
            self._status("Loading DQMQ data...")
            logger.info("DQMQ raw load started")
            file_path = self._get_current_file()
            if file_path is None:
                self._status("Select a DQMQ file first.")
                QMessageBox.warning(
                    self.parent,
                    "DQMQ file missing",
                    "Select a DQMQ file before loading DQMQ data.",
                    QMessageBox.Ok,
                )
                return

            try:
                raw_data = self._read_raw_data_from_file(file_path)
            except Exception as exc:
                logger.exception("DQMQ raw load failed")
                self._status(f"Could not process file: {exc}")
                QMessageBox.warning(
                    self.parent,
                    "Corrupt File",
                    f"Couldn't read the file because {exc}",
                    QMessageBox.Ok,
                )
                if self.parent is not None:
                    self.parent.clear_list()
                return

            self.reset_cached_results(clear_raw=False)
            self._write_raw_table(raw_data)
            self.render_raw_plot()
            logger.info("DQMQ raw load completed: %d points", len(raw_data["time"]))
            self._status("Loaded DQMQ data.")

    def render_raw_plot(self):
        raw_data = self._ensure_raw_data()
        if raw_data is None:
            return np.array([]), np.array([]), np.array([])

        self.current_plot_mode = "raw"
        figure = self._clear_signal_plot()
        time = raw_data["time"]
        dq = raw_data["dq_raw"]
        ref = raw_data["ref_raw"]

        if self._is_checked("DQMQ_CheckBox_DQ", default=True):
            figure.plot(
                time,
                dq,
                pen=self._dqmq_pen("dq"),
                name="Raw DQ")

        if self._is_checked("DQMQ_CheckBox_Ref", "DQMQ_CheckBox_Reference", default=True):
            figure.plot(
                time,
                ref,
                pen=self._dqmq_pen("ref"),
                name="Raw Ref"
                )

        return time, dq, ref

    def plot_original(self):
        result = self.render_raw_plot()
        self._status("Plot updated.")
        return result

    def _analysis_parameters(self):
        fit_from = self.ui.DQMQ_DoubleSpinBox_FitFrom.value()
        fit_to = self.ui.DQMQ_DoubleSpinBox_FitTo.value()
        fitting_exponent = self.ui.DQMQ_DoubleSpinBox_Power.value()
        baseline_level = self.ui.DQMQ_DoubleSpinBox_Noise.value()
        time_shift = int(self.ui.DQMQ_DoubleSpinBox_TimeShift.value())
        smoothing = [
            self.ui.DQMQ_DoubleSpinBox_SmoothFrom.value(),
            self.ui.DQMQ_DoubleSpinBox_SmoothTo.value(),
            int(self.ui.DQMQ_DoubleSpinBox_SmoothWindow.value()),
        ]
        return fit_from, fit_to, fitting_exponent, baseline_level, time_shift, smoothing

    def run_full_analysis(self):
        raw_data = self._ensure_raw_data()
        if raw_data is None:
            return None
        self._status("Running DQMQ analysis...")

        if self.ui.DQMQ_Table_Data.rowCount() != len(raw_data["time"]):
            self._write_raw_table(raw_data)

        fit_from, fit_to, fitting_exponent, baseline_level, time_shift, smoothing = (
            self._analysis_parameters()
        )
        logger.info(
            "DQMQ full analysis: from=%s to=%s exponent=%s baseline=%s time_shift=%s smoothing=%s",
            fit_from,
            fit_to,
            fitting_exponent,
            baseline_level,
            time_shift,
            smoothing,
        )

        try:
            with busy_cursor():
                self.analysis_result = dqmq_signal.calculate_dqmq_analysis(
                    raw_time=raw_data["time"],
                    dq_raw=raw_data["dq_raw"],
                    ref_raw=raw_data["ref_raw"],
                    fit_from=fit_from,
                    fit_to=fit_to,
                    fitting_exponent=fitting_exponent,
                    noise_level=baseline_level,
                    time_shift=time_shift,
                    smoothing=smoothing,
                )
        except Exception as exc:
            logger.exception("DQMQ full analysis failed")
            self._status(f"Could not run DQMQ analysis: {exc}")
            QMessageBox.warning(self.parent, "DQMQ analysis failed", str(exc))
            return None

        self._write_ndq_to_table(
            self.analysis_result["time"],
            self.analysis_result["nDQ"],
        )
        self.dres_is_stale = True
        self.render_analysis_plot()
        self._status("Fit completed.")
        return self.analysis_result

    def plot_norm(self):
        return self.run_full_analysis()

    def plot_diff(self):
        return self.run_full_analysis()

    def plot_nDQ(self):
        return self.run_full_analysis()

    def on_analysis_parameter_editing_finished(self):
        if self.raw_data is None and self._get_current_file() is None:
            return
        self.run_full_analysis()

    def on_visibility_checkbox_changed(self):
        if self.current_plot_mode == "raw":
            self.render_raw_plot()
            return
        if self.current_plot_mode == "analysis":
            self.render_analysis_plot()
            self._status("Plot updated.")

    def _dqmq_pen(self, style_name, *, dashed=False):
        style = DQMQ_PLOT_STYLES[style_name]
        qt_style = Qt.DashLine if dashed else Qt.SolidLine

        return mkPen(
            style["color"],
            width=style["width"],
            style=qt_style,
        )

    def render_analysis_plot(self):
        if self.analysis_result is None:
            self.render_raw_plot()
            return

        self.current_plot_mode = "analysis"
        figure = self._clear_signal_plot()
        result = self.analysis_result
        time = result["time"]

        if self._is_checked("DQMQ_CheckBox_DQ", default=True):
            figure.plot(
                time,
                result["dq_norm"],
                pen=self._dqmq_pen("dq"),
                name="DQ"
                )

        if self._is_checked("DQMQ_CheckBox_Ref", "DQMQ_CheckBox_Reference", default=True):
            figure.plot(
                time,
                result["ref_norm"],
                pen=self._dqmq_pen("ref"),
                name="Ref"
                )

        if self._is_checked("DQMQ_CheckBox_MQ", default=True):
            figure.plot(
                time,
                result["mq_norm"],
                pen=self._dqmq_pen("mq"),
                name="MQ"
                )

        if self._is_checked("DQMQ_CheckBox_Diff", "DQMQ_CheckBox_Difference", default=True):
            figure.plot(
                time,
                result["diff"],
                pen=self._dqmq_pen("diff"),
                name="Ref - DQ",
            )

        if self._is_checked("DQMQ_CheckBox_Fit", "DQMQ_CheckBox_TailFitting", default=False):
            figure.plot(
                time,
                result["tail_fit"],
                pen=self._dqmq_pen("tail_fit", dashed=True),
                name="Tail fit",
            )

        if self._is_checked("DQMQ_CheckBox_MQTailDiff", default=False):

            figure.plot(
                time,
                result["denominator_base_norm"],
                pen=self._dqmq_pen("mq_minus_tail"),
                name="MQ - Tail",
            )

        if self._is_checked("DQMQ_CheckBox_NDQ", default=True):
            figure.addItem(
                InfiniteLine(
                    pos=0.5,
                    angle=0,
                    pen=mkPen(color=(200, 200, 255), width=2, style=Qt.DashLine),
                )
            )

            figure.plot(
                result["time0"],
                result["nDQ"],
                pen=mkPen("k", width=1),
                symbol="o",
                symbolPen=DQMQ_PLOT_STYLES["ndq"]["color"],
                symbolBrush=DQMQ_PLOT_STYLES["ndq"]["color"],
                symbolSize=7,
                name="nDQ",
            )

        if self.integral_sum_result and self._is_checked("DQMQ_CheckBox_IntegralSum", default=False):
            figure.plot(
                self.integral_sum_result["time"],
                self.integral_sum_result["signal_norm"],
                pen=self._dqmq_pen("integral_sum"),
                name=f"Integral sum, shift={self.integral_sum_result['shift']:g}",
            )

        if self.dres_result and self._is_checked("DQMQ_CheckBox_DresFit", "DQMQ_CheckBox_DresFitting", default=False):
            figure.plot(
                self.dres_result["fit_x"],
                self.dres_result["fit_y"],
                pen=self._dqmq_pen("dres_fit", dashed=True),
                name="Dres fit",
            )


    def plot_nDQ_on_Load(self):
        try:
            arrays = self._read_dqmq_table_arrays()
        except ValueError:
            self.analysis_result = None
            self.render_raw_plot()
            return

        self.analysis_result = {
            "time": arrays["tau"],
            "dq_norm": arrays["DQ"],
            "ref_norm": arrays["Ref"],
            "diff": arrays["Ref"] - arrays["DQ"],
            "mq_raw_norm": arrays["DQ"] + arrays["Ref"],
            "mq_norm": arrays["DQ"] + arrays["Ref"],
            "tail_fit": np.zeros_like(arrays["tau"]),
            "mq_tail": np.zeros_like(arrays["tau"]),
            "noise_weight": np.zeros_like(arrays["tau"]),
            "additive": np.zeros_like(arrays["tau"]),
            "denominator_base": arrays["DQ"] + arrays["Ref"],
            "denominator_base_norm": np.zeros_like(arrays["tau"]),
            "denominator": arrays["DQ"] + arrays["Ref"],
            "time0": arrays["Time0"],
            "nDQ": arrays["nDQ0"],
            "fit_from": None,
            "fit_to": None,
            "power": None,
            "noise": None,
            "time_shift": None,
            "smoothing": None,
        }
        self.render_analysis_plot()

    def calculate_integral_sum(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self.parent, "Select DQ/T2 distribution CSV", "", "CSV Files (*.csv)"
        )
        if not file_path:
            return
        try:
            with busy_cursor():
                data = self._read_dqt2_distribution_file(file_path)
                t, signal_norm = self._calculate_t2_summed_signal(data)
                shift = float(self.ui.DQMQ_DoubleSpinBox_IntegralShift.value())
                self.integral_sum_result = {
                    "time": t + shift,
                    "time_base": t,
                    "signal_norm": signal_norm,
                    "source_file": file_path,
                    "shift": shift,
                }
                if self.current_plot_mode == "analysis":
                    self.render_analysis_plot()
                self._status("Calculated integral sum.")
        except ValueError as exc:
            self._status(f"Could not process file: {exc}")
            QMessageBox.warning(self.parent, "Invalid file", str(exc))
        except Exception as exc:
            logger.exception("Integral sum calculation failed")
            self._status(f"Could not calculate integral sum: {exc}")
            QMessageBox.warning(
                self.parent, "Integral sum calculation failed", str(exc)
            )

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

        finite_values = np.isfinite(amp_dq) & np.isfinite(t2) & np.isfinite(t_dq)
        positive_t2 = t2 > 0
        valid_rows = finite_values & positive_t2
        amp_dq = amp_dq[valid_rows]
        t2 = t2[valid_rows]
        t_dq = t_dq[valid_rows]

        if len(t2) == 0:
            raise ValueError(
                "Selected file does not contain valid positive T2star_lin values."
            )

        order = np.argsort(t2)
        amp_dq = amp_dq[order]
        t2 = t2[order]

        # t = np.linspace(0, np.max(t_dq), n_points)

        t_extra = np.linspace(0, np.max(t_dq), n_points)

        t = np.unique(np.concatenate([t_dq, t_extra]))
        t.sort()
        d_t2 = np.gradient(t2) if len(t2) > 1 else np.ones_like(t2)


        signal = np.zeros_like(t)
        for amp, t2i, delta_t2 in zip(amp_dq, t2, d_t2):
            scaled_time = t / t2i
            gaussian_decay = np.exp(-(scaled_time**2))
            signal += amp * gaussian_decay * delta_t2

        s_min = np.min(signal)
        s_max = np.max(signal)
        if s_max == s_min:
            signal_norm = np.zeros_like(signal)
        else:
            signal_range = s_max - s_min
            signal_norm = (signal - s_min) / signal_range
        return t, signal_norm

    def update_integral_sum_shift(self):
        if not self.integral_sum_result:
            return
        shift = float(self.ui.DQMQ_DoubleSpinBox_IntegralShift.value())
        self.integral_sum_result["shift"] = shift
        self.integral_sum_result["time"] = self.integral_sum_result["time_base"] + shift
        if self.current_plot_mode == "analysis":
            self.render_analysis_plot()
            self._status("Integral sum shift updated.")


    def calculate_dres(self):
        self.calculate_dres_distribution(show_errors=True)

    def calculate_dres_distribution(self, show_errors=True):
        try:
            with busy_cursor():
                arrays = self._dres_input_arrays()
                kernel = self._selected_dres_kernel()
                n_components = self._selected_dres_component_count()
                k_value = self._dres_k_value()
                beta = self._weibul_beta_value()
                l_value = self._dres_l_value()
                p0 = self._dres_initial_parameters(n_components)

                fitto_value = self._dres_fitto_value()
                fitto_idx = dqmq_signal._nearest_index(arrays['Time0'], fitto_value)
                array_x  = arrays['Time0'][:fitto_idx]
                array_y  = arrays['nDQ0'][:fitto_idx]

                fit_result = dqmq_dres.fit_selected_model(
                    arrays['Time0'],
                    array_x,
                    array_y,
                    kernel=kernel,
                    n_components=n_components,
                    p0=p0,
                    k_value=k_value,
                    beta = beta,
                    l_value = l_value
                )
                d_plot, p_dist = dqmq_dres.build_distribution(fit_result)

                self.dres_result = {
                    "D_plot": d_plot,
                    "P": p_dist,
                    "fit_x": fit_result["fit_x"],
                    "fit_y": fit_result["fit_y"],
                    "kernel": kernel,
                    "n_components": n_components,
                    "params": fit_result["popt"],
                    "param_names": fit_result["param_names"],
                    "p0": p0,
                    "k_value": k_value,
                    "beta_value" : beta,
                    "l_value":l_value,
                }
                self._write_fitted_dres_parameters_to_ui(
                    fit_result["popt"],
                    n_components,
                )
                self.dres_is_stale = False
                self._plot_dres_distribution()
                if self.current_plot_mode == "analysis":
                    self.render_analysis_plot()
                self._status("Dres calculation completed.")
        except ValueError as exc:
            if show_errors:
                self._status(f"Could not calculate Dres: {exc}")
                QMessageBox.warning(self.parent, "Dres calculation", str(exc))
            else:
                logger.warning("Skipping automatic Dres calculation: %s", exc)
        except Exception as exc:
            logger.exception("Dres calculation failed")
            if show_errors:
                self._status(f"Could not calculate Dres: {exc}")
                QMessageBox.warning(self.parent, "Dres calculation failed", str(exc))

    def _dres_input_arrays(self):
        if self.analysis_result is not None:
            arrays = {
                "tau": self.analysis_result["time"],
                "DQ": self.analysis_result["dq_norm"],
                "Ref": self.analysis_result["ref_norm"],
                "nDQ": self.analysis_result["nDQ"][1:],
                "Time0": self.analysis_result["time0"],
                "nDQ0": self.analysis_result["nDQ"],
            }
        else:
            arrays = self._read_dqmq_table_arrays()

        return self._finite_dres_arrays(arrays)

    def _finite_dres_arrays(self, arrays):
        time0 = np.asarray(arrays["Time0"], dtype=float)
        ndq0 = np.asarray(arrays["nDQ0"], dtype=float)
        valid_points = np.isfinite(time0) & np.isfinite(ndq0)
        time0 = time0[valid_points]
        ndq0 = ndq0[valid_points]
        if len(ndq0) < 3:
            raise ValueError("Need at least 3 finite nDQ values to calculate Dres.")

        return {
            **arrays,
            "Time0": time0,
            "nDQ0": ndq0,
        }

    def _read_dqmq_table_arrays(self):
        table = self.ui.DQMQ_Table_Data
        if table.rowCount() == 0:
            raise ValueError(
                "Run DQMQ analysis first so the table contains nDQ values."
            )

        values = {"tau": [], "dq": [], "ref": [], "nDQ": []}
        column_names = [
            ("tau", 0, "Time"),
            ("dq", 1, "DQ"),
            ("ref", 2, "Ref"),
            ("nDQ", 3, "nDQ"),
        ]
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
                    raise ValueError(
                        f"Row {row + 1} has a non-numeric {label} value."
                    ) from exc

            if not row_has_any_value:
                continue
            if any(row_values[key] is None for key, _, _ in column_names):
                raise ValueError(
                    f"Row {row + 1} is missing one or more DQMQ table values."
                )
            for key in values:
                values[key].append(row_values[key])

        if not values["nDQ"]:
            raise ValueError(
                "The DQMQ table nDQ column is empty. Run Plot Norm before calculating Dres."
            )
        if len(values["nDQ"]) < 3:
            raise ValueError(
                "Need at least 3 DQMQ table rows with nDQ values to calculate Dres."
            )

        tau = np.asarray(values["tau"], dtype=float)
        dq = np.asarray(values["dq"], dtype=float)
        ref = np.asarray(values["ref"], dtype=float)
        ndq_original = np.asarray(values["nDQ"], dtype=float)
        ndq = (
            savgol_filter(ndq_original, 3, 1)
            if len(ndq_original) >= 3
            else ndq_original
        )

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
            "P-L": "p-l",
        }
        selected = self.ui.DQMQ_ComboBox_Kernel.currentText()
        if selected not in kernel_map:
            raise ValueError(f"Unknown Dres kernel: {selected}")
        return kernel_map[selected]

    def _selected_dres_component_count(self):
        if self.ui.DQMQ_RadioButton_OneDres.isChecked():
            return 1
        if self.ui.DQMQ_RadioButton_TwoDres.isChecked():
            return 2
        raise ValueError("Select either 1 Dres or 2 Dres components.")

    def _dres_fitto_value(self):
        return self.ui.DQMQ_DoubleSpinBox_DresFitTo.value()

    def _weibul_beta_value(self):
        return self.ui.DQMQ_DoubleSpinBox_DresWeibullBeta.value()

    def _dres_k_value(self):
        return self.ui.DQMQ_DoubleSpinBox_DresK.value()

    def _dres_l_value(self):
        return self.ui.DQMQ_DoubleSpinBox_DresL.value()

    def _dres_initial_parameters(self, n_components):
        center1 = self.ui.DQMQ_DoubleSpinBox_DresCenter1.value()
        width1 = self.ui.DQMQ_DoubleSpinBox_DresWidth1.value()

        if n_components == 1:
            return [center1, width1]

        center2 = self.ui.DQMQ_DoubleSpinBox_DresCenter2.value()
        width2 = self.ui.DQMQ_DoubleSpinBox_DresWidth2.value()
        fraction1 = self.ui.DQMQ_DoubleSpinBox_DresFraction1.value()
        return [center1, width1, center2, width2, fraction1]

    def _write_fitted_dres_parameters_to_ui(self, fitted_parameters, n_components):
        parameter_widgets = [
            self.ui.DQMQ_DoubleSpinBox_DresCenter1,
            self.ui.DQMQ_DoubleSpinBox_DresWidth1,
        ]
        parameter_values = [
            fitted_parameters[0],
            fitted_parameters[1],
        ]
        if n_components == 2:
            parameter_widgets = [
                self.ui.DQMQ_DoubleSpinBox_DresCenter1,
                self.ui.DQMQ_DoubleSpinBox_DresWidth1,
                self.ui.DQMQ_DoubleSpinBox_DresCenter2,
                self.ui.DQMQ_DoubleSpinBox_DresWidth2,
                self.ui.DQMQ_DoubleSpinBox_DresFraction1,
            ]
            parameter_values = fitted_parameters

        for widget, value in zip(parameter_widgets, parameter_values):
            previous_state = widget.blockSignals(True)
            widget.setValue(float(value))
            widget.blockSignals(previous_state)

    def mark_dres_stale(self):
        self.dres_is_stale = True
        logger.info("DQMQ Dres parameters changed; Dres result is stale")
        self._status("Dres settings changed.")

    def _plot_dres_distribution(self):
        if not self.dres_result:
            return

        dres_figure = self.ui.DQMQ_PlotWidget_Dres
        dres_figure.clear()
        dres_figure.plot(
            self.dres_result["D_plot"],
            self.dres_result["P"],
            pen=mkPen("m", width=3),
            name="Dres distribution",
        )
        fitted_parameters = self.dres_result["params"]
        center_indices = [0]
        if self.dres_result["n_components"] == 2:
            center_indices.append(2)

        for center_index in center_indices:
            center = fitted_parameters[center_index]
            center_khz = center / (2 * np.pi) * 1000.0
            center_line = InfiniteLine(
                pos=center_khz,
                angle=90,
                pen=mkPen((80, 80, 80), width=2, style=Qt.DashLine),
            )
            dres_figure.addItem(center_line)

    def save_integral_sum_result(self, base_file_path):
        if not self.integral_sum_result:
            return

        root, ext = os.path.splitext(base_file_path)
        save_path = f"{root}_IntegralSum{ext or '.csv'}"
        integral = np.column_stack(
            (self.integral_sum_result["time"], self.integral_sum_result["signal_norm"])
            )

        mq = np.column_stack(
            (self.analysis_result["time"], self.analysis_result["mq_norm"], self.analysis_result["denominator_base_norm"])
            )

        can_write_excel = importlib.util.find_spec("pandas") is not None and (
            importlib.util.find_spec("openpyxl") is not None
            or importlib.util.find_spec("xlsxwriter") is not None
        )
        if can_write_excel:
            import pandas as pd

            try:
                with pd.ExcelWriter(f"{root}_InegralSum.xlsx") as writer:
                    pd.DataFrame(
                        integral, columns=["time", "integral_sum"]
                    ).to_excel(
                        writer,
                        sheet_name="Integral Sum",
                        index=False,
                    )
                    pd.DataFrame(mq, columns=["tauDQ", "MQ", "MQ_tail"]).to_excel(
                        writer,
                        sheet_name="MQ",
                        index=False,
                    )
                return
            except Exception:
                logger.exception("Excel IntegralSum export failed; writing CSV fallback files")


    def save_dres_result(self, base_file_path):
        if not self.dres_result:
            return
        root, _ = os.path.splitext(base_file_path)
        distribution = np.column_stack(
            (self.dres_result["D_plot"], self.dres_result["P"])
        )
        fit = np.column_stack((self.dres_result["fit_x"], self.dres_result["fit_y"]))
        metadata = {
            "kernel": self.dres_result["kernel"],
            "n_components": self.dres_result["n_components"],
            "k_value": self.dres_result.get("k_value"),
            "beta": self.dres_result.get("beta"),
            **dict(zip(self.dres_result["param_names"], self.dres_result["params"])),
        }

        can_write_excel = importlib.util.find_spec("pandas") is not None and (
            importlib.util.find_spec("openpyxl") is not None
            or importlib.util.find_spec("xlsxwriter") is not None
        )
        if can_write_excel:
            import pandas as pd

            try:
                with pd.ExcelWriter(f"{root}_Dres.xlsx") as writer:
                    pd.DataFrame(
                        distribution, columns=["Dres_over_2pi_kHz", "P_Dres"]
                    ).to_excel(
                        writer,
                        sheet_name="Dres distribution",
                        index=False,
                    )
                    pd.DataFrame(fit, columns=["time", "fitted_nDQ"]).to_excel(
                        writer,
                        sheet_name="Fit",
                        index=False,
                    )
                    pd.DataFrame(
                        metadata.items(), columns=["parameter", "value"]
                    ).to_excel(
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
        metadata_header = "\n".join(
            [f"{key}: {value}" for key, value in metadata.items()]
        )
        np.savetxt(
            f"{root}_Dres_fit.csv",
            fit,
            delimiter=",",
            header=f"{metadata_header}\ntime,fitted_nDQ",
            comments="# ",
        )

    def reset_Dres_values(self):
        logger.info("Dres fitting parameters restored to defaults")
        self._status("Dres fitting parameters restored to defaults.")
        self.ui.DQMQ_DoubleSpinBox_DresCenter1.setValue(0.250)
        self.ui.DQMQ_DoubleSpinBox_DresWidth1.setValue(0.001)
        self.ui.DQMQ_DoubleSpinBox_DresCenter2.setValue(0.05)
        self.ui.DQMQ_DoubleSpinBox_DresWidth2.setValue(0.001)
        self.ui.DQMQ_DoubleSpinBox_DresFraction1.setValue(0.5)
        self.ui.DQMQ_DoubleSpinBox_DresWeibullBeta.setValue(2)
        self.ui.DQMQ_DoubleSpinBox_DresK.setValue(0.4)
        self.ui.DQMQ_DoubleSpinBox_DresL.setValue(0.4)