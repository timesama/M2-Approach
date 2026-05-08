from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QFileDialog, QMessageBox
from pyqtgraph import mkPen

from calculations import recfid_signal as recfid
from controllers.base_tab_controller import BaseTabController
from dialogs.recfid_dialogs import RecFIDTimeAnalysisDialog
from dialogs.recfid_echo_time_dialog import RecFIDEchoTimeDialog
from dialogs.open_files_dialog import OpenFilesDialog
import dialogs.open_files_dialog as open_files_dialog_module
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


class RecFIDController(BaseTabController):
    """UI orchestration for the g_RecFID tab."""

    def __init__(self, ui, state, parent=None):
        super().__init__(ui, state, parent)
        self.reset_tab(clear_plots=False)

    def reset_tab(self, clear_plots=True):
        self.recfid_mode = "none"
        self.selected_fid_files = []
        self.selected_data_files = []
        self.selected_fid_empty_files = []
        self.selected_data_empty_files = []
        self.echo_time_values = []
        self.manual_echo_time_enabled = False
        self.fid_result = {}
        self.data_results = {}
        self.se_max_result = {}
        self.extrapolation = None
        self.build_result = {}
        self.time_analysis_result = {}
        self._clear_original_signal_gradient_label()
        self._close_time_analysis_dialog()
        if clear_plots:
            self.initialize_plots()
            self._clear_result_text_widgets()

    def connect_signals(self):
        self._connect("RecFID_Button_LoadFID", "clicked", self.load_fid_files)
        self._connect("RecFID_Button_LoadData", "clicked", self.load_data_files)
        self._connect("RecFID_Button_LoadFIDEmpty", "clicked", self.load_empty_fid_files)
        self._connect("RecFID_Button_LoadDataEmpty", "clicked", self.load_empty_data_files)
        self._connect("RecFID_Button_Clear", "clicked", self.reset_tab)
        self._connect("RecFID_Button_ManualEchoTimes", "clicked", self.open_manual_echo_time_dialog)
        self._connect("RecFID_Button_SaveStatistics", "clicked", self.save_statistics)
        self._connect("RecFID_Button_ExportSEMaxT2", "clicked", self.export_se_max_t2)
        self._connect("RecFID_Button_Run", "clicked", self.run_full_pipeline)
        self._connect("RecFID_Button_TimeAnalysis", "clicked", self.run_time_analysis)
        for widget_name in (
            "RecFID_DoubleSpinBox_SEFitFrom",
            "RecFID_DoubleSpinBox_SEFitTo",
            "RecFID_DoubleSpinBox_BuildFrom",
            "RecFID_DoubleSpinBox_BuildTo",
            "RecFID_DoubleSpinBox_MseDivider",
        ):
            self._connect(widget_name, "editingFinished", lambda *_args: self.rebuild(show_warning=False))
        self._connect("RecFID_ComboBox_BuildFunction", "activated", lambda *_args: self.rebuild(show_warning=False))

    def initialize_plots(self):
        self._setup_graph("RecFID_PlotWidget_OriginalNMRSignal", "Time, μs", "Amplitude", "FID / data")
        self._setup_graph("RecFID_PlotWidget_SEMax", "Echo Time, μs", "Amplitude", "SE Max")
        self._setup_graph("RecFID_PlotWidget_BuildUpNMRSignal", "Time, μs", "Amplitude", "FID Build")
        self._setup_graph("RecFID_PlotWidget_BuildUpSpectra", "Frequency, MHz", "Amplitude, a.u", "FFT Build")

    def load_fid_files(self):
        files = self._select_files(multiple=False)
        if files is None:
            return
        self.selected_fid_files = files
        self._status(f"Loaded {len(files)} RecFID FID file(s).")

    def load_data_files(self):
        files = self._select_files(multiple=True)
        if files is None:
            return
        self.selected_data_files = files
        self.recfid_mode = "mse" if len(files) == 1 else "se" if len(files) > 1 else "none"
        self.echo_time_values = []
        self.manual_echo_time_enabled = False
        self._status(f"Loaded {len(files)} RecFID data file(s); mode={self.recfid_mode}.")

    def load_empty_fid_files(self):
        files = self._select_files(multiple=False)
        if files is None:
            return
        self.selected_fid_empty_files = files
        self._status(f"Loaded {len(files)} RecFID empty FID file(s).")

    def load_empty_data_files(self):
        files = self._select_files(multiple=True)
        if files is None:
            return
        self.selected_data_empty_files = files
        self._status(f"Loaded {len(files)} RecFID empty data file(s).")

    def run_full_pipeline(self):
        with busy_cursor():
            if not self._validate_required_files() or not self._validate_empty_data_files():
                return
            try:
                if self.run_basic_data_analysis() is None:
                    return
                if self.extrapolate() is None:
                    return
                self.build_up()
            except Exception as exc:
                logger.exception("RecFID analysis failed")
                self._warn("Couldn't analyse", f"Couldn't analyse the RecFID data because {exc}.")

    def rebuild(self, show_warning=True):
        if not self.fid_result or not self.selected_data_files:
            return
        with busy_cursor():
            try:
                if self.extrapolate() is None:
                    return
                self.build_up()
            except Exception as exc:
                logger.exception("RecFID rebuild failed")
                if show_warning:
                    self._warn("Couldn't rebuild", f"Couldn't rebuild the RecFID result because {exc}.")

    def run_basic_data_analysis(self):
        data, data_empty, fid, fid_empty = self._choose_files_for_comparison(0)
        result = recfid.analyze_signal(data, fid, data_empty, fid_empty, self._analysis_options())
        result = self._apply_mse_divider_to_analysis_result(result)
        if result is None:
            return None
        self.fid_result = {
            "Time_td_fid": result.time_fid,
            "Re_td": result.signal_fid,
            "Freq": result.frequency_fid,
            "Real_fft": result.spectrum_fid,
            "M2": result.m2_fid,
            "T2": result.t2_fid,
        }
        self.data_results = {
            data: self._data_result_from_analysis(result),
        }
        self._plot_all_processed_signals()
        self._status(
            f"FID T₂*={result.t2_fid}, M₂={result.m2_fid}; data T₂*={result.t2_data}, M₂={result.m2_data}"
        )
        return result

    def run_se_analysis(self):
        self.data_results = {}
        se_max = []
        se_t2 = []
        echo_time = self._echo_times_for_se_files()
        fid = self.selected_fid_files[0]
        fid_empty = self.selected_fid_empty_files[0] if self.selected_fid_empty_files and self._checked("RecFID_CheckBox_SubtractEmpty", True) else None

        for index, file_name in enumerate(self.selected_data_files):
            data_empty = self.selected_data_empty_files[index] if self.selected_data_empty_files else None
            result = recfid.analyze_signal(file_name, fid, data_empty, fid_empty, self._analysis_options())
            self.data_results[file_name] = self._data_result_from_analysis(result)
            se_max.append(float(np.max(result.signal_data)))
            se_t2.append(result.t2_data)

        extrapolation, fitting_curve, fittingx = recfid.find_maximum_se(
            echo_time,
            se_max,
            self._value("RecFID_DoubleSpinBox_SEFitFrom", 0.0),
            self._value("RecFID_DoubleSpinBox_SEFitTo", 0.0),
        )
        self.se_max_result = {
            "echo_time": np.asarray(echo_time, dtype=float),
            "maximum": np.asarray(se_max, dtype=float),
            "t2": np.asarray(se_t2, dtype=float),
            "fit_x": fittingx,
            "fit_y": fitting_curve,
        }
        self._plot_all_processed_signals()
        return extrapolation

    def extrapolate(self):
        if self.recfid_mode == "se":
            extrapolation = self.run_se_analysis()
            self.extrapolation = extrapolation
            self._plot_se_maximum()
        elif self.recfid_mode == "mse":
            result = self.run_basic_data_analysis()
            if result is None:
                self.extrapolation = None
                self.se_max_result = {}
                self.build_result = {}
                return None
            raw_maximum = float(np.max(result.signal_data))
            self.extrapolation = raw_maximum
            self.se_max_result = {
                "echo_time": np.array([0.0]),
                "maximum": np.array([raw_maximum]),
                "t2": np.array([result.t2_data]),
            }
            self._clear_or_hide("RecFID_PlotWidget_SEMax", hide=True)
        else:
            raise ValueError("Load one data file for MSE mode or multiple data files for SE mode before running.")

        self._status(f"RecFID max = {round(self.extrapolation, 2)}")
        return self.extrapolation

    def build_up(self):
        if self.extrapolation is None:
            raise ValueError("Run extrapolation before build-up.")
        if not self.fid_result:
            raise ValueError("Run FID/data analysis before build-up.")
        begin = self._value("RecFID_DoubleSpinBox_BuildFrom", 9.0)
        finish = self._value("RecFID_DoubleSpinBox_BuildTo", 20.0)
        time_build, data_build, data_fit, freq_build, spectrum_build, m2_build, t2_build = recfid.analyze_build_up(
            self.fid_result["Time_td_fid"],
            self.fid_result["Re_td"],
            self.extrapolation,
            self._text("RecFID_ComboBox_BuildFunction", "Gaussian"),
            begin,
            finish,
            self._value("RecFID_DoubleSpinBox_ApodizationSigma", 100.0),
        )
        self.build_result = {
            "Time": time_build,
            "Re": data_build,
            "Fit": data_fit,
            "Freq": freq_build,
            "Real_fft": spectrum_build,
            "M2": float(m2_build),
            "T2": float(t2_build),
        }
        self._plot_build_up(freq_build, spectrum_build, time_build, data_build)
        self._status(
            f"Build-up T₂*={round(t2_build, 5)}, M₂={round(m2_build, 5)}; "
            f"FID T₂*={self.fid_result['T2']}, M₂={self.fid_result['M2']}"
        )
        return m2_build, t2_build

    def run_time_analysis(self):
        if not self.fid_result or self.extrapolation is None:
            self._warn("No build-up data", "Run RecFID analysis before analysing time ranges.")
            return
        with busy_cursor():
            start_range, finish_range, m2_values, t2_values = self._calculate_time_analysis(
                self._text("RecFID_ComboBox_BuildFunction", "Gaussian")
            )
        dialog = RecFIDTimeAnalysisDialog(self.parent)
        dialog.plot_data(start_range, finish_range, t2_values, m2_values)
        dialog.setWindowTitle("RecFID T2 Map")
        dialog.resize(800, 600)
        dialog.show()
        self._plot_window = dialog
        self.time_analysis_result = {
            "start_range": start_range,
            "finish_range": finish_range,
            "m2": m2_values,
            "t2": t2_values,
        }

    def save_statistics(self):
        if not self.fid_result or self.extrapolation is None:
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
        combo = getattr(self.ui, "RecFID_ComboBox_BuildFunction")
        with busy_cursor():
            for index in range(combo.count()):
                combo.setCurrentIndex(index)
                start_range, finish_range, _m2, t2 = self._calculate_time_analysis(combo.currentText())
                ranges = [(begin, finish) for begin in start_range for finish in finish_range]
                self._write_frequencies_t2(ranges, t2.flatten(), file_path, combo.currentText())
        self._status(f"Saved RecFID statistics to {file_path}")

    def export_se_max_t2(self):
        if self.extrapolation is None or not self.se_max_result:
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
            handle.write(f"Max Extrapolation: {self.extrapolation}\n")
            handle.write("Echo Time\tT2\tRatio\tEcho Amp\n")
            for echo_time, maximum, t2 in zip(
                self.se_max_result.get("echo_time", []),
                self.se_max_result.get("maximum", []),
                self.se_max_result.get("t2", []),
            ):
                ratio = self.extrapolation / maximum if maximum else 0
                handle.write(f"{echo_time}\t{t2}\t{ratio}\t{maximum}\t\n")
            handle.write("\n")
        self._status(f"Saved RecFID SE maxima and T2 values to {file_path}")

    def open_manual_echo_time_dialog(self):
        if not self.selected_data_files:
            self._warn("No SE files", "Load (M)SE files before assigning echo times manually.")
            return
        dialog = RecFIDEchoTimeDialog(self.selected_data_files, self.echo_time_values, self.parent)
        if dialog.exec():
            try:
                self.echo_time_values = dialog.echo_times()
                self.manual_echo_time_enabled = True
            except ValueError as exc:
                self._warn("Invalid echo time", str(exc))

    def save_results(self):
        if not self.build_result:
            self._warn("No RecFID result", "Run RecFID analysis before saving reconstructed FID results.")
            return
        directory = QFileDialog.getExistingDirectory(self.parent, "Select RecFID export directory")
        if not directory:
            return
        base = Path(directory)
        np.savetxt(
            base / "recfid_built_fid.csv",
            np.column_stack((self.build_result["Time"], self.build_result["Re"])),
            delimiter=",",
            header="Time,Re",
            comments="",
        )
        np.savetxt(
            base / "recfid_built_fft.csv",
            np.column_stack((self.build_result["Freq"], self.build_result["Real_fft"])),
            delimiter=",",
            header="Frequency,Amplitude",
            comments="",
        )
        for widget_name, file_name in (
            ("RecFID_PlotWidget_BuildUpNMRSignal", "recfid_built_fid.png"),
            ("RecFID_PlotWidget_BuildUpSpectra", "recfid_built_fft.png"),
        ):
            widget = getattr(self.ui, widget_name, None)
            if widget is not None:
                exporter = pyqtgraph.exporters.ImageExporter(widget.plotItem)
                exporter.export(str(base / file_name))
        self._status(f"Saved RecFID results to {directory}")

    # Backward-compatible method aliases used by previous controller wiring/tests.
    reset = reset_tab
    clear_all = reset_tab
    run = run_full_pipeline
    compare = run_basic_data_analysis
    se_analysis = run_se_analysis
    manual_assignment = open_manual_echo_time_dialog
    time_analysis = run_time_analysis
    exportSEMaxT2 = export_se_max_t2

    def _calculate_time_analysis(self, function_name):
        start_range, finish_range = recfid.time_range_grid(
            self.fid_result["Time_td_fid"], self._value("RecFID_SpinBox_TimeAnalysisRange", 24)
        )
        m2 = []
        t2 = []
        for begin in start_range:
            for finish in finish_range:
                if finish < begin + 3:
                    m2_build = 0
                    t2_build = 0
                else:
                    try:
                        _tb, _db, _df, _freq, _spectrum, m2_build, t2_build = recfid.analyze_build_up(
                            self.fid_result["Time_td_fid"],
                            self.fid_result["Re_td"],
                            self.extrapolation,
                            function_name,
                            begin,
                            finish,
                            self._value("RecFID_DoubleSpinBox_ApodizationSigma", 100.0),
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
            handle.write(f"Name: {self.selected_fid_files}\n")
            handle.write(f"Mode: {self.recfid_mode}\n")
            handle.write(f"Function: {function_name}\n\n")
            for range_bf, t2 in zip(ranges, t2s):
                handle.write(f"{range_bf[0]}\t{range_bf[1]}\t{t2}\t\n")
            handle.write("\n")

    def _analysis_options(self):
        return recfid.AnalysisOptions(
            subtract_empty=self._checked("RecFID_CheckBox_SubtractEmpty", True),
            cut_beginning=self._checked("RecFID_CheckBox_CutBeginning", True),
            normalize_to_fid=self._checked("RecFID_CheckBox_NormalizeToFID", True),
            normalize_from=self._value("RecFID_DoubleSpinBox_NormalizeFrom", 70.0),
            normalize_to=self._value("RecFID_DoubleSpinBox_NormalizeTo", 90.0),
            long_component=self._checked("RecFID_CheckBox_LongComponent", False),
            long_component_from=self._value("RecFID_DoubleSpinBox_LongComponentEnd", 55.0),
            apodize_time_domain=self._checked("RecFID_CheckBox_TimeDomainApodization", True),
            apodization_time=self._value("RecFID_DoubleSpinBox_ApodizationSigma", 100.0),
            adjust_frequency_phase=True,
            adjust_fid_zero=self._checked("RecFID_CheckBox_AdjustZero", False),
            fid_zero_shift=self._value("RecFID_DoubleSpinBox_ZeroShift", 0.0),
            smooth=self._checked("RecFID_CheckBox_Smooth", False),
            smooth_order=int(self._value("RecFID_SpinBox_SmoothOrder", 1)),
            smooth_window=int(self._value("RecFID_SpinBox_SmoothWindow", 5)),
        )

    def _choose_files_for_comparison(self, number):
        if self.recfid_mode == "se":
            data = self.selected_data_files[number]
            data_empty = self.selected_data_empty_files[number] if self.selected_data_empty_files else None
        else:
            data = self.selected_data_files[0]
            data_empty = self.selected_data_empty_files[0] if self.selected_data_empty_files else None
        fid_empty = self.selected_fid_empty_files[0] if self.selected_fid_empty_files else None
        return data, data_empty, self.selected_fid_files[0], fid_empty

    def _echo_times_for_se_files(self):
        if self.manual_echo_time_enabled:
            return self.echo_time_values
        echo_times = []
        try:
            for file_name in self.selected_data_files:
                echo_times.append(recfid.extract_echo_time(file_name))
        except ValueError:
            self._warn(
                "Could not parse echo time",
                "Could not parse echo time from filenames. Please enter echo times manually.",
            )
            self.open_manual_echo_time_dialog()
            if not self.manual_echo_time_enabled:
                raise
            echo_times = self.echo_time_values
        return echo_times

    def _validate_required_files(self):
        if not self.selected_fid_files:
            self._warn("No FID file", "Load a FID file first.")
            return False
        if not self.selected_data_files:
            self._warn("No (M)SE file", "Load (M)SE data first.")
            return False
        if self.recfid_mode not in {"mse", "se"}:
            self.recfid_mode = "mse" if len(self.selected_data_files) == 1 else "se" if len(self.selected_data_files) > 1 else "none"
        return True

    def _validate_empty_data_files(self):
        if not self._checked("RecFID_CheckBox_SubtractEmpty", True) or not self.selected_data_empty_files:
            return True
        if self.recfid_mode == "mse" and len(self.selected_data_empty_files) == 1:
            return True
        if self.recfid_mode == "se" and len(self.selected_data_empty_files) == len(self.selected_data_files):
            return True
        self._warn(
            "Different amount of files",
            "The number of empty data files must match the number of data files when subtracting empty data. Empty data files were cleared.",
        )
        self.selected_data_empty_files = []
        return True


    def _apply_mse_divider_to_analysis_result(self, result):
        if self.recfid_mode != "mse":
            return result
        divider = self._validated_mse_divider()
        if divider is None:
            return None
        # MSE divider applies only to the MSE data amplitude, not the FID reference.
        return recfid.SignalAnalysisResult(
            time_data=result.time_data,
            signal_data=result.signal_data / divider,
            time_fid=result.time_fid,
            signal_fid=result.signal_fid,
            frequency_data=result.frequency_data,
            spectrum_data=result.spectrum_data / divider,
            frequency_fid=result.frequency_fid,
            spectrum_fid=result.spectrum_fid,
            m2_data=result.m2_data,
            t2_data=result.t2_data,
            m2_fid=result.m2_fid,
            t2_fid=result.t2_fid,
        )

    def _data_result_from_analysis(self, result):
        return {
            "Time_td": result.time_data,
            "Re_td": result.signal_data,
            "Freq": result.frequency_data,
            "Real_fft": result.spectrum_data,
            "M2": result.m2_data,
            "T2": result.t2_data,
        }

    def _validated_mse_divider(self):
        divider = self._value("RecFID_DoubleSpinBox_MseDivider", 2.0)
        if not np.isfinite(divider) or divider <= 0:
            self._warn("Invalid MSE divider", "MSE divider must be greater than zero.")
            return None
        return divider

    def _plot_all_processed_signals(self):
        if not self.fid_result:
            return
        graph_nmr = getattr(self.ui, "RecFID_PlotWidget_OriginalNMRSignal", None)
        if graph_nmr is None:
            return
        self._prepare_plot(graph_nmr)
        graph_nmr.plot(
            self.fid_result["Time_td_fid"],
            self.fid_result["Re_td"],
            pen=mkPen(QColor(255, 0, 0), width=4),
            symbol=None,
        )
        for index, result in enumerate(self.data_results.values()):
            graph_nmr.plot(
                result["Time_td"],
                result["Re_td"],
                pen=mkPen(self._winter_color(index, max(len(self.data_results), 1)), width=3),
                symbol=None,
            )
        # Gradient label is the compact visible color-scale indicator for data curves.
        self._add_original_signal_gradient_label(graph_nmr)
        self._status("Original NMR plot: black = FID; data curves use a winter color scale by echo time/file order.")

    def _add_original_signal_gradient_label(self, graph):
        if not self.data_results:
            return

        self._clear_original_signal_gradient_label()

        scale_text = self._original_signal_color_scale_text()
        html = (
            '<div style="background-color:rgba(255,255,255,210); padding:4px; border:1px solid #444;">'
            '<span style="color:#000000; font-weight:bold;">FID: black</span><br>'
            '<span style="color:#0033ff; font-weight:bold;">■</span>'
            '<span style="color:#007fbf; font-weight:bold;">■</span>'
            '<span style="color:#00bf80; font-weight:bold;">■</span>'
            '<span style="color:#00ff7f; font-weight:bold;">■</span>'
            f' Data winter scale<br>{scale_text}'
            '</div>'
        )

        label = pg.TextItem(html=html, anchor=(1, 0))
        view_box = graph.getViewBox()
        label.setParentItem(view_box)
        label.setZValue(10_000)

        def update_position():
            rect = view_box.sceneBoundingRect()
            label.setPos(rect.width() - 10, 10)

        update_position()
        view_box.sigResized.connect(update_position)

        self.original_signal_gradient_label = label
        self.original_signal_gradient_label_update = update_position

    def _original_signal_color_scale_text(self):
        if self.recfid_mode == "se" and self.se_max_result.get("echo_time") is not None:
            echo_times = np.asarray(self.se_max_result.get("echo_time", []), dtype=float)
            if echo_times.size > 1 and np.all(np.isfinite(echo_times)):
                return f"Echo time: {np.min(echo_times):g} → {np.max(echo_times):g}"
        if len(self.data_results) > 1:
            return f"File order: 1 → {len(self.data_results)}"
        return "Data"

    def _clear_original_signal_gradient_label(self):
        label = getattr(self, "original_signal_gradient_label", None)
        if label is not None:
            graph = getattr(self.ui, "RecFID_PlotWidget_OriginalNMRSignal", None)
            if graph is not None:
                try:
                    graph.removeItem(label)
                except (RuntimeError, ValueError):
                    logger.debug("RecFID gradient label was already removed", exc_info=True)
            self.original_signal_gradient_label = None

    def _prepare_plot(self, graph):
        graph.clear()
        graph.show()
        legend = getattr(graph.plotItem, "legend", None)
        if legend is not None:
            legend.clear()

    def _winter_color(self, index, total):
        fraction = 0.0 if total <= 1 else index / (total - 1)
        green = int(round(255 * fraction))
        blue = int(round(255 * (1 - 0.5 * fraction)))
        return QColor(0, green, blue)

    def _plot_se_maximum(self):
        graph = getattr(self.ui, "RecFID_PlotWidget_SEMax", None)
        if graph is None:
            return
        graph.clear()
        graph.show()
        graph.plot([0], [self.extrapolation], pen="r", symbolPen=None, symbol="o", symbolBrush="r")
        graph.plot(self.se_max_result["fit_x"], self.se_max_result["fit_y"], pen=mkPen("k", width=3, style=Qt.DashLine), symbol=None)
        graph.plot(
            self.se_max_result["echo_time"],
            self.se_max_result["maximum"],
            pen=None,
            symbolPen=None,
            symbol="o",
            symbolBrush="b",
        )

    def _plot_build_up(self, freq_build, spectrum_build, time_build, data_build):
        graph_build_fid = getattr(self.ui, "RecFID_PlotWidget_BuildUpNMRSignal", None)
        if graph_build_fid is not None:
            graph_build_fid.clear()
            graph_build_fid.plot(time_build, data_build, pen=mkPen("b", width=4), symbol=None)
            graph_build_fid.plot(self.fid_result["Time_td_fid"], self.fid_result["Re_td"], pen=mkPen("r", width=3), symbol=None)
        graph_build_fft = getattr(self.ui, "RecFID_PlotWidget_BuildUpSpectra", None)
        if graph_build_fft is not None:
            graph_build_fft.clear()
            graph_build_fft.plot(freq_build, spectrum_build, pen=mkPen("b", width=4), symbol=None)
            graph_build_fft.plot(self.fid_result["Freq"], self.fid_result["Real_fft"], pen=mkPen("r", width=3), symbol=None)

    def _clear_or_hide(self, widget_name, hide=False):
        widget = getattr(self.ui, widget_name, None)
        if widget is not None:
            widget.clear()
            if hide:
                widget.hide()

    def _select_files(self, multiple):
        open_files_dialog_module.State_multiple_files = multiple
        dialog = OpenFilesDialog(self.parent)
        if dialog.exec():
            return dialog.selectedFiles()
        return None

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
        if widget_name == "RecFID_PlotWidget_OriginalNMRSignal":
            self._clear_original_signal_gradient_label()
        widget.clear()
        legend = getattr(widget.plotItem, "legend", None)
        if legend is not None:
            legend.clear()
        widget.show()
        widget.getAxis("left").setLabel(ylabel)
        widget.getAxis("bottom").setLabel(xlabel)
        widget.setTitle(title)

    def _clear_result_text_widgets(self):
        for widget_name in (
            "RecFID_TextEdit_OriginalResults",
            "RecFID_TextEdit_SEMaxResult",
            "RecFID_TextEdit_BuildUpResults",
        ):
            widget = getattr(self.ui, widget_name, None)
            if widget is not None:
                widget.clear()

    def _close_time_analysis_dialog(self):
        dialog = getattr(self, "_plot_window", None)
        if dialog is not None:
            dialog.close()
            self._plot_window = None

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
