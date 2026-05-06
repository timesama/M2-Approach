import os
import logging
import numpy as np
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem
from scipy.signal import savgol_filter
import Calculator as Cal
from controllers.base_tab_controller import BaseTabController
from dialogs.open_files_dialog import OpenFilesDialog
import dialogs.phasing_manual as phasing_manual_module
from dialogs.phasing_manual import PhasingManual
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


class GeneralSEDQController(BaseTabController):
    def __init__(self, ui, state, parent=None, se_controller=None, dq_controller=None):
        super().__init__(ui, state, parent)
        self.se_controller = se_controller
        self.dq_controller = dq_controller

    def analysis(self):
        with busy_cursor():
            mw = self.parent
            while self.ui.comboBox_4.count() > 0:
                self.ui.comboBox_4.removeItem(0)

            if mw.tab == "SE":
                files = mw.selected_files
                self.ui.SE_PlotWidget_Main.clear()
                self.ui.SE_ComboBox_YAxis.setCurrentIndex(-1)
            else:
                files = mw.selected_files_DQ_single
                self.ui.DQ_PlotWidget_T2.clear()
                self.ui.DQ_PlotWidget_NormIntensity.clear()
                self.ui.DQ_TextEdit_FitResult.setText("")
                self.ui.DQ_ComboBox_FitFunction.setCurrentIndex(-1)

            if len(files) == 0:
                return

            logger.info("%s analysis started: %d files", mw.tab, len(files))
            logger.info(
                "Options - glycerol:%s baseline:%s long_component:%s smoothing:%s",
                self.ui.checkBox_glycerol.isChecked(),
                self.ui.checkBox_baseline.isChecked(),
                self.ui.checkBox_long_component.isChecked(),
                self.ui.checkBox_Smooth.isChecked(),
            )

            if self.ui.checkBox_glycerol.isChecked() and mw.selected_files_gly == []:
                self.open_select_dialog_glycerol()

            if self.ui.checkBox_baseline.isChecked() and mw.selected_files_empty == []:
                self.open_select_dialog_baseline()

            mw.disable_buttons()
            self.ui.btn_SelectFiles.setEnabled(False)
            self.ui.btn_Load.setEnabled(False)
            self.ui.radioButton.setEnabled(False)
            self.ui.comboBox_4.setCurrentIndex(-1)
            if self.ui.checkBox_Smooth.isChecked():
                mw.window_array = np.linspace(
                    self.ui.SmoothWindowFrom.value(),
                    self.ui.SmoothWindowTo.value(),
                    len(files),
                    dtype=np.int32,
                )
            else:
                mw.window_array = np.array([])

            for i, file_path in enumerate(files, start=1):
                logger.info(
                    "%s processing file %d/%d: %s",
                    mw.tab,
                    i,
                    len(files),
                    os.path.basename(file_path),
                )

                try:
                    file_path_gly = mw.selected_files_gly[i - 1] if self.ui.checkBox_glycerol.isChecked() else []
                    file_path_empty = mw.selected_files_empty[i - 1] if self.ui.checkBox_baseline.isChecked() else []
                    self.process_file_data(file_path, file_path_gly, file_path_empty, i)
                except Exception:
                    self.analysis_error(file_path, files)

            self.update_legends_and_dq_graphs()
            self.ui.btn_Start.setStyleSheet("background-color: none")
            logger.info("%s analysis completed: %d files", mw.tab, len(files))

    def analysis_error(self, file_path, files):
        logger.exception("Failed to process file %s", file_path)

        mw = self.parent
        QMessageBox.warning(
            mw,
            "Invalid Data",
            f"Couldn't read the {os.path.basename(file_path)}. Deleting the file.",
            QMessageBox.Ok,
        )

        files.remove(file_path)
        self.ui.btn_SelectFiles.setEnabled(True)
        self.ui.btn_Start.setStyleSheet("background-color: none")
        self.ui.FidWidget.clear()
        self.ui.FFTWidget.clear()
        self.ui.btn_Phasing.setEnabled(False)
        mw.enable_buttons()

    def update_legends_and_dq_graphs(self):
        self.parent.enable_buttons()

        if self.parent.tab == "DQ":
            self.dq_controller.update_graphs()

    def process_file_data(self, file_path, file_path_gly, file_path_empty, i):
        mw = self.parent
        filename = os.path.basename(file_path)
        subtract = self.ui.checkBox_baseline.isChecked()
        data = np.loadtxt(file_path)
        x, y, z = data[:, 0], data[:, 1], data[:, 2]
        time, re_signal, im_signal = Cal.analysis_time_domain(file_path, file_path_empty, subtract)

        if self.ui.checkBox_glycerol.isChecked():
            time_reference, re_reference, im_reference = Cal.analysis_time_domain(file_path_gly, [], False)
            time, re_signal, im_signal = Cal.magnet_inhomogenity_correction(
                time,
                time_reference,
                re_signal,
                re_reference,
                im_signal,
                im_reference,
            )

        if self.ui.checkBox_long_component.isChecked():
            time, re_signal, im_signal = Cal.subtract_long_component(time, re_signal, im_signal)

        amplitude = Cal._calculate_amplitude(re_signal, im_signal)
        mw.update_graphs(time, amplitude, re_signal, im_signal, self.ui.FidWidget)

        time_fid, fid = Cal.final_analysis_time_domain(time, re_signal, im_signal, 2**16)
        frequency = Cal._calculate_frequency_scale(time_fid)
        fft = np.fft.fftshift(np.fft.fft(fid))

        if mw.window_array.size != 0:
            window = mw.window_array[i - 1]
            real_part = savgol_filter(np.real(fft), window, 1)
            imaginary_part = savgol_filter(np.imag(fft), window, 1)
            fft = np.array(real_part + 1j * imaginary_part)

        logger.info("SE/DQ frequency-domain processing: calculating apodization and M2/T2")
        amp_spectra, re_spectra, im_spectra = Cal._simple_baseline_correction(fft)

        if mw.tab == "SE":
            phased_store = mw.phased_spectra_SE
        else:
            phased_store = mw.phased_spectra_DQ

        phased_record = phased_store.get(filename)
        if phased_record:
            re_spectra = np.array(phased_record.get("re", re_spectra))
            im_spectra = np.array(phased_record.get("im", im_spectra))
            amp_spectra = Cal._calculate_amplitude(re_spectra, im_spectra)

        real_apod = Cal._calculate_apodization(re_spectra, frequency)
        m2, t2 = Cal._calculate_M2(real_apod, frequency)

        self.state.spectrum.frequency = frequency
        self.state.spectrum.re_spectra = re_spectra
        self.state.spectrum.im_spectra = im_spectra
        mw.update_graphs(frequency, amp_spectra, re_spectra, im_spectra, self.ui.FFTWidget)

        if self.ui.comboBox_4.findText(filename) == -1:
            self.ui.comboBox_4.addItem(filename)

        if self.ui.comboBox_4.currentIndex() == -1:
            if mw.tab == "SE":
                self.se_controller.process_processed_file(i, filename, amplitude, m2, t2, file_path)
            elif mw.tab == "DQ":
                self.dq_controller.process_processed_file(i, filename, x, y, z, m2, t2, file_path)

    def after_phasing(self):
        mw = self.parent
        i = self.ui.comboBox_4.currentIndex()

        if i < 0:
            return

        frequency = self.state.spectrum.frequency
        re_spectra = self.state.spectrum.re_spectra
        im_spectra = self.state.spectrum.im_spectra
        phased_window = getattr(mw, "phasing_manual_window", None)

        if phased_window is not None and getattr(phased_window, "Real_freq_phased", None) is not None:
            re_spectra = phased_window.Real_freq_phased
            self.state.spectrum.re_spectra = re_spectra
            self.state.spectrum.im_spectra = np.zeros_like(re_spectra)
            im_spectra = self.state.spectrum.im_spectra

        filename = self.ui.comboBox_4.currentText()

        if mw.tab == "SE" and filename:
            mw.phased_spectra_SE[filename] = {"re": re_spectra.tolist(), "im": im_spectra.tolist()}
        elif mw.tab == "DQ" and filename:
            mw.phased_spectra_DQ[filename] = {"re": re_spectra.tolist(), "im": im_spectra.tolist()}

        real_apod = Cal._calculate_apodization(re_spectra, frequency)
        amp_spectra = Cal._calculate_amplitude(re_spectra, im_spectra)

        mw.update_graphs(frequency, amp_spectra, re_spectra, im_spectra, self.ui.FFTWidget)

        m2, t2 = Cal._calculate_M2(real_apod, frequency)

        if mw.tab == "SE":
            table = self.ui.SE_Table_Data
        else:
            table = self.ui.DQ_Table_Data

        table.setItem(i, 2, QTableWidgetItem(str(round(m2, 6))))
        table.setItem(i, 3, QTableWidgetItem(str(round(t2, 3))))

        if mw.tab == "SE":
            self.se_controller.update_graphs()
        elif mw.tab == "DQ":
            self.dq_controller.update_graphs()

    def open_phasing_manual(self):
        mw = self.parent
        phasing_manual_module.Frequency = self.state.spectrum.frequency
        phasing_manual_module.Re_spectra = self.state.spectrum.re_spectra
        phasing_manual_module.Im_spectra = self.state.spectrum.im_spectra
        mw.phasing_manual_window = PhasingManual()
        mw.phasing_manual_window.read_data()
        mw.phasing_manual_window.show()
        mw.phasing_manual_window.closed.connect(self.after_phasing)

    def open_select_dialog_glycerol(self):
        dlg = OpenFilesDialog(self.parent)
        dlg.setWindowTitle("Select Glycerol Files")
        if dlg.exec():
            self.parent.selected_files_gly.extend(dlg.selectedFiles())
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)

    def open_select_dialog_baseline(self):
        dlg = OpenFilesDialog(self.parent)
        dlg.setWindowTitle("Select Baseline Files")
        if dlg.exec():
            self.parent.selected_files_empty.extend(dlg.selectedFiles())
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)
