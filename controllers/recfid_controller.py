from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pyqtgraph.exporters
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFileDialog, QMessageBox
from pyqtgraph import mkPen

from calculations import recfid_signal as recfid
from controllers.base_tab_controller import BaseTabController
from dialogs.recfid_dialogs import RecFIDManualEchoTimeDialog, RecFIDTimeAnalysisDialog
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


@dataclass
class RecFIDState:
    fid_files: list[str] = field(default_factory=list)
    data_files: list[str] = field(default_factory=list)
    fid_empty_files: list[str] = field(default_factory=list)
    data_empty_files: list[str] = field(default_factory=list)
    manual_echo_times: list[float] = field(default_factory=list)
    manual_assignment: bool = False
    mode: str = "none"
    fid_result: dict = field(default_factory=dict)
    data_results: dict = field(default_factory=dict)
    se_max: np.ndarray = field(default_factory=lambda: np.array([]))
    se_t2: np.ndarray = field(default_factory=lambda: np.array([]))
    echo_time_array: np.ndarray = field(default_factory=lambda: np.array([]))
    extrapolation: float | None = None
    build_result: dict = field(default_factory=dict)


class RecFIDController(BaseTabController):
    """UI orchestration for the g_RecFID tab."""

    def __init__(self, ui, state, parent=None):
        super().__init__(ui, state, parent)
        self.recfid_state = RecFIDState()
        if not hasattr(self.state, "recfid"):
            self.state.recfid = self.recfid_state

    def connect_signals(self):
        self._connect("pushButton", "clicked", lambda: self.select_files("fid"))
        self._connect("pushButton_2", "clicked", lambda: self.select_files("data"))
        self._connect("pushButton_3", "clicked", lambda: self.select_files("fid_empty"))
        self._connect("pushButton_4", "clicked", lambda: self.select_files("data_empty"))
        self._connect("pushButton_5", "clicked", self.manual_assignment)
        self._connect("Btn_SaveStatistics", "clicked", self.save_statistics)
        self._connect("pushButtonExportMaxT2", "clicked", self.exportSEMaxT2)
        self._connect("pushButtonGO", "clicked", self.run)
        self._connect("pushButton_analyse", "clicked", self.time_analysis)
        for widget_name in ("doubleSpinBox_7", "doubleSpinBox_6", "comboBox", "doubleSpinBox_begin", "doubleSpinBox_finish"):
            signal = "activated" if widget_name == "comboBox" else "valueChanged"
            self._connect(widget_name, signal, lambda *_args: self.rebuild(show_warning=False))

    def initialize_plots(self):
        self._setup_graph("widgetOriginalNMRSignal", "Time, μs", "Amplitude", "FID / data")
        self._setup_graph("widgetSEMax", "Echo Time, μs", "Amplitude", "SE Max")
        self._setup_graph("widgetBUNMRSignal", "Time, μs", "Amplitude", "FID Build")
        self._setup_graph("widgetBUFFT", "Frequency, MHz", "Amplitude, a.u", "FFT Build")

    def reset(self):
        self.recfid_state = RecFIDState()
        self.state.recfid = self.recfid_state
        self.initialize_plots()

    def select_files(self, kind):
        multiple = kind in {"data", "data_empty"}
        file_mode = QFileDialog.ExistingFiles if multiple else QFileDialog.ExistingFile
        dialog = QFileDialog(self.parent)
        dialog.setFileMode(file_mode)
        dialog.setNameFilter("Data (*.dat *.txt *.csv)")
        if dialog.exec():
            files = dialog.selectedFiles()
        else:
            return
        target = {
            "fid": self.recfid_state.fid_files,
            "data": self.recfid_state.data_files,
            "fid_empty": self.recfid_state.fid_empty_files,
            "data_empty": self.recfid_state.data_empty_files,
        }[kind]
        target.clear()
        target.extend(files)
        if kind == "data":
            self._update_mode_from_data_files()
            self.recfid_state.manual_assignment = False
            self.recfid_state.manual_echo_times = []
        self._status(f"Loaded {len(files)} RecFID {kind.replace('_', ' ')} file(s).")

    def run(self):
        with busy_cursor():
            if not self._validate_required_files():
                return
            try:
                self.compare()
                self.extrapolate()
                begin = self._value("doubleSpinBox_begin", 9.0)
                finish = self._value("doubleSpinBox_finish", 20.0)
                self.build_up(begin, finish)
            except Exception as exc:
                logger.exception("RecFID analysis failed")
                self._warn("Couldn't analyse", f"Couldn't analyse the RecFID data because {exc}.")

    def rebuild(self, show_warning=True):
        if not self.recfid_state.fid_result or not self.recfid_state.data_files:
            return
        with busy_cursor():
            try:
                self.extrapolate()
                self.build_up(self._value("doubleSpinBox_begin", 9.0), self._value("doubleSpinBox_finish", 20.0))
            except Exception as exc:
                logger.exception("RecFID rebuild failed")
                if show_warning:
                    self._warn("Couldn't rebuild", f"Couldn't rebuild the RecFID result because {exc}.")

    def compare(self):
        data, data_empty, fid, fid_empty = self._choose_files_for_comparison(0)
        result = recfid.analyze_signal(data, fid, data_empty, fid_empty, self._analysis_options())
        self.recfid_state.fid_result = {
            "Time_td_fid": result.time_fid,
            "Re_td": result.signal_fid,
            "Freq": result.frequency_fid,
            "Real_fft": result.spectrum_fid,
            "M2": result.m2_fid,
            "T2": result.t2_fid,
        }
        graph_nmr = getattr(self.ui, "widgetOriginalNMRSignal", None)
        if graph_nmr is not None:
            graph_nmr.clear()
            graph_nmr.plot(result.time_data, result.signal_data, pen=mkPen("r", width=3), symbol=None)
            graph_nmr.plot(result.time_fid, result.signal_fid, pen=mkPen("b", width=3), symbol=None)
        self._status(
            f"FID T₂*={result.t2_fid}, M₂={result.m2_fid}; data T₂*={result.t2_data}, M₂={result.m2_data}"
        )
        return result

    def se_analysis(self, start, finish):
        state = self.recfid_state
        state.data_results = {}
        se_max = []
        se_t2 = []
        echo_time = self._echo_times_for_se_files()
        fid = state.fid_files[0]
        fid_empty = state.fid_empty_files[0] if state.fid_empty_files and self._checked("checkBox_subempty", True) else None
        if state.data_empty_files and len(state.data_empty_files) != len(state.data_files):
            raise ValueError("Please load the same number of data and empty data files.")
        for idx, file_name in enumerate(state.data_files):
            data_empty = state.data_empty_files[idx] if state.data_empty_files else None
            result = recfid.analyze_signal(file_name, fid, data_empty, fid_empty, self._analysis_options())
            state.data_results[file_name] = {
                "Time_td": result.time_data,
                "Re_td": result.signal_data,
                "Freq": result.frequency_data,
                "Real_fft": result.spectrum_data,
                "M2": result.m2_data,
                "T2": result.t2_data,
            }
            se_max.append(float(np.max(result.signal_data)))
            se_t2.append(result.t2_data)
        extrapolation, fitting_curve, fittingx = recfid.find_maximum_se(echo_time, se_max, start, finish)
        state.se_max = np.asarray(se_max, dtype=float)
        state.se_t2 = np.asarray(se_t2, dtype=float)
        state.echo_time_array = np.asarray(echo_time, dtype=float)
        return extrapolation, state.echo_time_array, state.se_max, fittingx, fitting_curve

    def extrapolate(self):
        state = self.recfid_state
        graph_se = getattr(self.ui, "widgetSEMax", None)
        self._update_mode_from_data_files()
        if state.mode == "se":
            start = self._value("doubleSpinBox_7", 0.0)
            end = self._value("doubleSpinBox_6", 0.0)
            extrapolation, echo_time, maximum_se, fittingx, fittingy = self.se_analysis(start, end)
            state.extrapolation = extrapolation
            if graph_se is not None:
                graph_se.clear()
                graph_se.show()
                graph_se.plot([0], [extrapolation], pen="r", symbolPen=None, symbol="o", symbolBrush="r")
                graph_se.plot(fittingx, fittingy, pen=mkPen("k", width=3, style=Qt.DashLine), symbol=None)
                graph_se.plot(echo_time, maximum_se, pen=None, symbolPen=None, symbol="o", symbolBrush="b")
        else:
            result = self.compare()
            extrapolation = float(np.max(result.signal_data))
            if state.mode == "mse":
                extrapolation /= max(self._value("doubleSpinBox_8", 2.0), 1.0)
            state.extrapolation = extrapolation
            state.se_max = np.array([float(np.max(result.signal_data))])
            state.se_t2 = np.array([result.t2_data])
            state.echo_time_array = np.array([0.0])
            state.data_results[state.mode.upper()] = {
                "Time_td": result.time_data,
                "Re_td": result.signal_data,
                "Freq": result.frequency_data,
                "Real_fft": result.spectrum_data,
                "M2": result.m2_data,
                "T2": result.t2_data,
            }
            if graph_se is not None:
                graph_se.clear()
                graph_se.hide()
        self._status(f"RecFID max = {round(state.extrapolation, 2)}")
        return state.extrapolation

    def build_up(self, begin, finish):
        state = self.recfid_state
        if state.extrapolation is None:
            raise ValueError("Run extrapolation before build-up.")
        fid = state.fid_result
        if not fid:
            raise ValueError("Run FID/data analysis before build-up.")
        time_build, data_build, data_fit, freq_build, spectrum_build, m2_build, t2_build = recfid.analyze_build_up(
            fid["Time_td_fid"],
            fid["Re_td"],
            state.extrapolation,
            self._text("comboBox", "Gaussian"),
            begin,
            finish,
            self._value("doubleSpinBox_5", 100.0),
        )
        state.build_result = {
            "Time": time_build,
            "Re": data_build,
            "Fit": data_fit,
            "Freq": freq_build,
            "Real_fft": spectrum_build,
            "M2": float(m2_build),
            "T2": float(t2_build),
        }
        graph_build_fid = getattr(self.ui, "widgetBUNMRSignal", None)
        if graph_build_fid is not None:
            graph_build_fid.clear()
            graph_build_fid.plot(time_build, data_build, pen=mkPen("b", width=3), symbol=None)
            graph_build_fid.plot(fid["Time_td_fid"], fid["Re_td"], pen=mkPen("r", width=3), symbol=None)
        graph_build_fft = getattr(self.ui, "widgetBUFFT", None)
        if graph_build_fft is not None:
            graph_build_fft.clear()
            graph_build_fft.plot(freq_build, spectrum_build, pen=mkPen("b", width=3), symbol=None)
            graph_build_fft.plot(fid["Freq"], fid["Real_fft"], pen=mkPen("r", width=3), symbol=None)
        self._status(
            f"Build-up T₂*={round(t2_build, 5)}, M₂={round(m2_build, 5)}; FID T₂*={fid['T2']}, M₂={fid['M2']}"
        )
        return m2_build, t2_build

    def manual_assignment(self):
        if not self.recfid_state.data_files:
            self._warn("No SE files", "Load (M)SE files before assigning echo times manually.")
            return
        dialog = RecFIDManualEchoTimeDialog(
            self.recfid_state.data_files, self.recfid_state.manual_echo_times, self.parent
        )
        if dialog.exec():
            try:
                self.recfid_state.manual_echo_times = dialog.echo_times()
                self.recfid_state.manual_assignment = True
            except ValueError as exc:
                self._warn("Invalid echo time", str(exc))

    def time_analysis(self):
        if not self.recfid_state.fid_result:
            self._warn("No build-up data", "Run RecFID analysis before analysing time ranges.")
            return
        with busy_cursor():
            start_range, finish_range, m2_values, t2_values = self._calculate_time_analysis(self._text("comboBox", "Gaussian"))
        dialog = RecFIDTimeAnalysisDialog(self.parent)
        dialog.plot_data(start_range, finish_range, t2_values, m2_values)
        dialog.setWindowTitle("RecFID T2 Map")
        dialog.resize(800, 600)
        dialog.show()
        self._plot_window = dialog

    def save_statistics(self):
        if not self.recfid_state.fid_result:
            self._warn("No build-up data", "Run RecFID analysis before exporting a T2 map.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self.parent,
            "Save Frequency and T2 Data",
            "statistical_analysis.txt",
            "Text Files (*.txt);;All Files (*)",
            options=QFileDialog.DontConfirmOverwrite,
        )
        if not file_path:
            return
        with busy_cursor():
            for index in range(getattr(self.ui, "comboBox").count()):
                self.ui.comboBox.setCurrentIndex(index)
                start_range, finish_range, _m2, t2 = self._calculate_time_analysis(self.ui.comboBox.currentText())
                ranges = [(begin, finish) for begin in start_range for finish in finish_range]
                self._write_frequencies_t2(ranges, t2.flatten(), file_path, self.ui.comboBox.currentText())
        self._status(f"Saved RecFID statistics to {file_path}")

    def exportSEMaxT2(self):
        state = self.recfid_state
        if state.extrapolation is None or state.se_max.size == 0:
            self._warn("No SE maxima", "Run RecFID analysis before exporting SE maxima and T2 values.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self.parent,
            "Save Max Amplitudes of SE and T2*",
            "MaxSeT2.txt",
            "Text Files (*.txt);;All Files (*)",
            options=QFileDialog.DontConfirmOverwrite,
        )
        if not file_path:
            return
        with open(file_path, "a", encoding="utf-8") as handle:
            handle.write(f"Max Extrapolation: {state.extrapolation}\n")
            handle.write("Echo Time\tT2\tRatio\tEcho Amp\n")
            for echo_time, maximum, t2 in zip(state.echo_time_array, state.se_max, state.se_t2):
                ratio = state.extrapolation / maximum if maximum else 0
                handle.write(f"{echo_time}\t{t2}\t{ratio}\t{maximum}\t\n")
            handle.write("\n")
        self._status(f"Saved RecFID SE maxima and T2 values to {file_path}")

    def save_results(self):
        if not self.recfid_state.build_result:
            self._warn("No RecFID result", "Run RecFID analysis before saving reconstructed FID results.")
            return
        directory = QFileDialog.getExistingDirectory(self.parent, "Select RecFID export directory")
        if not directory:
            return
        base = Path(directory)
        build = self.recfid_state.build_result
        np.savetxt(base / "recfid_built_fid.csv", np.column_stack((build["Time"], build["Re"])), delimiter=",", header="Time,Re", comments="")
        np.savetxt(base / "recfid_built_fft.csv", np.column_stack((build["Freq"], build["Real_fft"])), delimiter=",", header="Frequency,Amplitude", comments="")
        for widget_name, file_name in (("widgetBUNMRSignal", "recfid_built_fid.png"), ("widgetBUFFT", "recfid_built_fft.png")):
            widget = getattr(self.ui, widget_name, None)
            if widget is not None:
                exporter = pyqtgraph.exporters.ImageExporter(widget.plotItem)
                exporter.export(str(base / file_name))
        self._status(f"Saved RecFID results to {directory}")

    def _calculate_time_analysis(self, function_name):
        fid = self.recfid_state.fid_result
        start_range, finish_range = recfid.time_range_grid(fid["Time_td_fid"], self._value("spinBoxTimeRange", 24))
        m2 = []
        t2 = []
        for begin in start_range:
            for finish in finish_range:
                if finish < begin + 3:
                    m2_build = 0
                    t2_build = 0
                else:
                    try:
                        _tb, _db, _df, freq, spectrum, m2_build, t2_build = recfid.analyze_build_up(
                            fid["Time_td_fid"], fid["Re_td"], self.recfid_state.extrapolation, function_name, begin, finish, self._value("doubleSpinBox_5", 100.0)
                        )
                    except Exception:
                        logger.exception("RecFID time-analysis point failed: begin=%s finish=%s", begin, finish)
                        m2_build = 0
                        t2_build = 0
                m2.append(m2_build)
                t2.append(t2_build)
        return (
            start_range,
            finish_range,
            np.asarray(m2).reshape(len(start_range), len(finish_range)),
            np.asarray(t2).reshape(len(start_range), len(finish_range)),
        )

    def _write_frequencies_t2(self, ranges, t2s, file_path, function_name):
        with open(file_path, "a", encoding="utf-8") as handle:
            handle.write(f"Name: {self.recfid_state.fid_files}\n")
            handle.write(f"Mode: {self.recfid_state.mode}\n")
            handle.write(f"Function: {function_name}\n\n")
            for range_bf, t2 in zip(ranges, t2s):
                handle.write(f"{range_bf[0]}\t{range_bf[1]}\t{t2}\t\n")
            handle.write("\n")

    def _analysis_options(self):
        return recfid.AnalysisOptions(
            subtract_empty=self._checked("checkBox_subempty", True),
            cut_beginning=self._checked("checkBox_cutbegin", True),
            normalize_to_fid=self._checked("checkBox_normtofid", True),
            normalize_from=self._value("doubleSpinBox", 70.0),
            normalize_to=self._value("doubleSpinBox_2", 90.0),
            long_component=self._checked("checkBox_longcomp", False),
            long_component_from=self._value("doubleSpinBox_4", 55.0),
            apodize_time_domain=self._checked("checkBox_td_apodization", True),
            apodization_time=self._value("doubleSpinBox_5", 100.0),
            adjust_frequency_phase=True,
            adjust_fid_zero=self._checked("checkBox_adjust_zero", False),
            fid_zero_shift=self._value("doubleSpinBox_3", 0.0),
            smooth=self._checked("checkBoxSmooth", False),
            smooth_order=int(self._value("spinBoxOrder", 1)),
            smooth_window=int(self._value("spinBoxWindow", 5)),
        )

    def _choose_files_for_comparison(self, number):
        state = self.recfid_state
        if state.mode == "se":
            data = state.data_files[number]
            data_empty = state.data_empty_files[number] if state.data_empty_files else None
        else:
            data = state.data_files[0]
            data_empty = state.data_empty_files[0] if state.data_empty_files else None
        return data, data_empty, state.fid_files[0], state.fid_empty_files[0] if state.fid_empty_files else None

    def _echo_times_for_se_files(self):
        state = self.recfid_state
        if state.manual_assignment:
            return state.manual_echo_times
        echo_time = []
        try:
            for file_name in state.data_files:
                echo_time.append(recfid.extract_echo_time(file_name))
        except ValueError:
            self.manual_assignment()
            if not state.manual_assignment:
                raise
            echo_time = state.manual_echo_times
        return echo_time

    def _update_mode_from_data_files(self):
        count = len(self.recfid_state.data_files)
        if count > 4:
            self.recfid_state.mode = "se"
        elif count == 1:
            self.recfid_state.mode = "mse"
        elif count > 1:
            self.recfid_state.mode = "single_se"
        else:
            self.recfid_state.mode = "none"

    def _validate_required_files(self):
        if not self.recfid_state.fid_files:
            self._warn("No FID file", "Load a FID file first.")
            return False
        if not self.recfid_state.data_files:
            self._warn("No (M)SE file", "Load (M)SE data first.")
            return False
        return True

    def _connect(self, widget_name, signal_name, slot):
        widget = getattr(self.ui, widget_name, None)
        if widget is None:
            logger.warning("RecFID widget %s is missing; signal not connected", widget_name)
            return
        signal = getattr(widget, signal_name, None)
        if signal is None:
            logger.warning("RecFID widget %s has no %s signal", widget_name, signal_name)
            return
        signal.connect(slot)

    def _setup_graph(self, widget_name, xlabel, ylabel, title):
        widget = getattr(self.ui, widget_name, None)
        if widget is None:
            return
        widget.clear()
        widget.getAxis("left").setLabel(ylabel)
        widget.getAxis("bottom").setLabel(xlabel)
        widget.setTitle(title)

    def _checked(self, widget_name, default=False):
        widget = getattr(self.ui, widget_name, None)
        return widget.isChecked() if widget is not None else default

    def _value(self, widget_name, default=0.0):
        widget = getattr(self.ui, widget_name, None)
        return widget.value() if widget is not None else default

    def _text(self, widget_name, default=""):
        widget = getattr(self.ui, widget_name, None)
        return widget.currentText() if widget is not None else default

    def _warn(self, title, message):
        QMessageBox.warning(self.parent, title, message, QMessageBox.Ok)

    def _status(self, message):
        if self.parent is not None and hasattr(self.parent, "statusBar"):
            self.parent.statusBar().showMessage(message, 8000)
        logger.info(message)
