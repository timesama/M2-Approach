# This Python file uses the following encoding: utf-8
import sys, os, re

from PySide6.QtWidgets import QPushButton, QCheckBox, QHBoxLayout, QVBoxLayout, QLabel, QSlider, QApplication, QMainWindow, QScrollArea, QFileDialog, QMessageBox, QVBoxLayout, QDialog, QTableWidgetItem
from PySide6.QtCore import Qt
import pyqtgraph as pg
from pyqtgraph import mkPen
import pyqtgraph.exporters
import numpy as np
# from pyqtgraph import setConfigOptions
# setConfigOptions(useOpenGL=False)

from ui_form import Ui_Main
import Math as math
from ui_manual import Ui_Form

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

pg.CONFIG_OPTIONS['background'] = 'w'
pg.CONFIG_OPTIONS['foreground'] = 'k'

class Main(QMainWindow):
    State_multiple_files = False

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_Main()
        self.ui.setupUi(self)


        # Set window geometry
        screen = QApplication.primaryScreen()

        if screen:
            self.showMaximized()

        self.scroll = QScrollArea(self)
        self.scroll.setWidget(self.ui.mainwidget)
        self.scroll.setWidgetResizable(True)
        self.setCentralWidget(self.scroll)
        # Maximize minimize
        self.setWindowFlags(self.windowFlags() | self.windowFlags().WindowMaximizeButtonHint)

        # Prepare the widnow - clean, initialize, hide, etc
        self.start()
        self.ui.pushButton_manual.clicked.connect(self.manual_assignment)

        # Events load data
        self.ui.pushButtonFID.clicked.connect(lambda: self.load_data("fid"))
        self.ui.pushButtonData.clicked.connect(lambda: self.load_data("none"))
        self.ui.pushButtonFIDEmpty.clicked.connect(lambda: self.load_data("fid_empty"))
        self.ui.pushButtonDataEmpty.clicked.connect(lambda: self.load_data("data_empty"))

        self.ui.pushButtonClear.clicked.connect(self.start)

        # Events analysis
        self.ui.pushButtonGO.clicked.connect(self.run)
        self.ui.pushButton_analyse.clicked.connect(self.time_analysis)
        self.ui.radioButtonSE.clicked.connect(self.change_state) #TODO
        self.ui.radioButtonMSE.clicked.connect(self.change_state) #TODO

        # Save for statisticss
        self.ui.Btn_SaveStatistics.clicked.connect(self.save_statistics)
        self.ui.pushButtonExportMaxT2.clicked.connect(self.exportSEMaxT2)

        # Rerun analysis on
        self.ui.doubleSpinBox_7.valueChanged.connect(self.rebuild)
        self.ui.doubleSpinBox_6.valueChanged.connect(self.rebuild)
        self.ui.comboBox.activated.connect(self.rebuild)
        self.ui.doubleSpinBox_begin.valueChanged.connect(self.rebuild)
        self.ui.doubleSpinBox_finish.valueChanged.connect(self.rebuild)

    def start(self):
        # Initializae array
        self.selected_files_fid = []
        self.selected_files_se = []
        self.selected_files_mse = []
        self.selected_files_fid_empty = []
        self.selected_files_se_empty = []
        self.selected_files_mse_empty = []
        self.data_echo_time = []

        #Initialize dictionary
        self.fid_dict = {}
        self.data_dict = {}
        self.se_maxT2_dict = {}

        # Initialize state se or mse
        self.state = 'none'
        self.manual_assignment_state = False

        # Hide what should be hidden
        # Disable what should be disabled
        self.ui.label.hide()
        self.ui.groupBox.setEnabled(True)
        self.ui.widgetSEMax.show()
        self.ui.doubleSpinBox_6.setEnabled(True)
        self.ui.doubleSpinBox_7.setEnabled(True)
        self.ui.groupBox_6.setEnabled(False)

        # Clear graphs
        self.setup_graph(self.ui.widgetOriginalSpectra, "Frequency, MHz", "Amplitude, a.u", "FFT")
        self.setup_graph(self.ui.widgetOriginalNMRSignal, "Time, μs", "Amplitude", "FID")
        self.setup_graph(self.ui.widgetBUSpectra, "Frequency, MHz", "Amplitude, a.u", "FFT Build")
        self.setup_graph(self.ui.widgetSEMax, "Echo Time, μs", "Amplitude", "SE Max")
        self.setup_graph(self.ui.widgetBUNMRSignal, "Time, μs", "Amplitude", "FID Build")

        #Clear text
        self.ui.textEdit.setText("")
        self.ui.textEdit_2.setText("")
        self.ui.textEdit_3.setText("")

    def setup_graph(self, graph_widget, xlabel="", ylabel="", title=""):
        graph_widget.clear()
        graph_widget.getAxis('left').setLabel(ylabel)
        graph_widget.getAxis('bottom').setLabel(xlabel)
        graph_widget.setTitle(title)

    def change_state(self):
        if self.ui.radioButtonSE.isChecked():
            self.state = "se"
        elif self.ui.radioButtonMSE.isChecked():
            self.state = 'mse'
        else:
            return
        
        print(self.state)

    def load_data(self, state):
        if self.ui.radioButtonSE.isChecked():
            self.state = "se"
        elif self.ui.radioButtonMSE.isChecked():
            self.state = 'mse'
        else:
            return

        if state == 'fid':
            Main.State_multiple_files = False
            selected_data=self.selected_files_fid
        elif state == 'none' and self.state == 'se':
            state == 'se'
            Main.State_multiple_files = True
            selected_data=self.selected_files_se
        elif state == 'none' and self.state == 'mse':
            state == 'mse'
            Main.State_multiple_files = False
            selected_data=self.selected_files_mse
        elif state == 'fid_empty':
            Main.State_multiple_files = False
            selected_data=self.selected_files_fid_empty
        elif state == 'data_empty' and self.state == 'se':
            Main.State_multiple_files = True
            selected_data=self.selected_files_se_empty
        elif state == 'data_empty' and self.state == 'mse':
            Main.State_multiple_files = False
            selected_data=self.selected_files_mse_empty
        else:
            return

        dlg = OpenFilesDialog(self)
        if dlg.exec():
            Files = dlg.selectedFiles()
        selected_data.clear()
        selected_data.extend(Files)

        if len(self.selected_files_fid)>0 and self.state!='none':
            self.ui.pushButtonGO.setEnabled(True)

    def checkstate_analysis(self):
        # Read analysis params
        time_a = self.ui.doubleSpinBox_5.value()
        n_start = self.ui.doubleSpinBox.value()
        n_end = self.ui.doubleSpinBox_2.value()
        end = self.ui.doubleSpinBox_4.value()
        zero = self.ui.doubleSpinBox_3.value()
        SmoothOrder = self.ui.spinBoxOrder.value()
        SmoothWindow = self.ui.spinBoxWindow.value()

        Subtraction = self.ui.checkBox_subempty.isChecked()
        Cut_beginning = self.ui.checkBox_cutbegin.isChecked()
        Normalize_to_fid = self.ui.checkBox_normtofid.isChecked()
        Long_component = self.ui.checkBox_longcomp.isChecked()
        Adjust = self.ui.checkBox_adjust_freq_phase.isChecked()
        Apodize = self.ui.checkBox_td_apodization.isChecked()
        AdjustZero = self.ui.checkBox_adjust_zero.isChecked()
        Smooth = self.ui.checkBoxSmooth.isChecked()


        # if self.ui.checkBox_subempty.isChecked():
        #     Subtraction = True
        # else:
        #     Subtraction = False

        # if self.ui.checkBox_cutbegin.isChecked():
        #     Cut_beginning = True
        # else:
        #     Cut_beginning = False

        # if self.ui.checkBox_normtofid.isChecked():
        #     Normalize_to_fid = True
        # else:
        #     Normalize_to_fid = False

        # if self.ui.checkBox_longcomp.isChecked():
        #     Long_component = True
        # else:
        #     Long_component = False

        # if self.ui.checkBox_adjust_freq_phase.isChecked():
        #     Adjust = True
        # else:
        #     Adjust = False

        # if self.ui.checkBox_td_apodization.isChecked():
        #     Apodize = True
        # else:
        #     Apodize = False

        # if self.ui.checkBox_adjust_zero.isChecked():
        #     AdjustZero = True
        # else:
        #     AdjustZero = False

        return  Subtraction, Cut_beginning, Normalize_to_fid, n_start, n_end,  Long_component, end, Apodize, time_a, Adjust, AdjustZero, zero, Smooth, SmoothOrder, SmoothWindow

    def choose_files_for_comparison(self, number):
        # Read FID and 1! File from Se or MSE
        if self.state == 'mse':
            data = self.selected_files_mse[0]
            if self.selected_files_mse_empty!=[]:
                data_empty = self.selected_files_mse_empty[0]
            else:
                data_empty = self.selected_files_mse_empty
        elif self.state == 'se' or self.state == 'single_se':
            data = self.selected_files_se[number]
            if self.selected_files_se_empty!=[]:
                data_empty = self.selected_files_se_empty[number]
            else:
                data_empty = self.selected_files_se_empty

        data_fid = self.selected_files_fid[0]
        if self.selected_files_fid_empty!=[]:
            data_fid_empty = self.selected_files_fid_empty[0]
        else:
            data_fid_empty = self.selected_files_fid_empty

        return data, data_empty, data_fid, data_fid_empty

    def basic_data_analysis(self):
        Subtraction, Cut_beginning, Normalize_to_fid, n_start, n_end,  Long_component, end, Apodize, time_a, Adjust, AdjustZero, zero, Smooth, SmoothOrder, SmoothWindow = self.checkstate_analysis()
        data, data_empty, data_fid, data_fid_empty = self.choose_files_for_comparison(0)

        Time_td, Re_td, Time_td_fid, Re_td_fid = math.nmr_signal_correction(data, data_fid, data_empty, data_fid_empty, Subtraction, Cut_beginning, Normalize_to_fid, n_start, n_end,  Long_component, end, AdjustZero, zero, Smooth, SmoothOrder, SmoothWindow)

        # Try for the sake of beauty
        Time_td, Re_td, Time_td_fid, Re_td_fid = math.for_the_sake_of_beauty(Time_td, Re_td, Time_td_fid, Re_td_fid, Apodize, time_a)

        Frequency, Real_fft = math.freq_domain_correction(Time_td, Re_td, 0, Apodize, time_a, Adjust)
        Frequency_fid, Real_fft_fid = math.freq_domain_correction(Time_td_fid, Re_td_fid, 0, Apodize, time_a, Adjust)

        M2, T2          = math.calculate_M2(Real_fft, Frequency)
        M2_fid, T2_fid  = math.calculate_M2(Real_fft_fid, Frequency_fid)

        M2 = round(M2, 5)
        T2 = round(T2, 5)
        M2_fid = round(M2_fid, 5)
        T2_fid = round(T2_fid, 5)

        self.fid_dict = {
            'Time_td_fid': Time_td_fid,
            'Re_td': Re_td_fid,
            'Freq': Frequency_fid,
            'Real_fft': Real_fft_fid,
            'M2': M2_fid,
            'T2': T2_fid
        }

        return Time_td, Re_td,Time_td_fid, Re_td_fid,Frequency, Real_fft,Frequency_fid, Real_fft_fid,M2,T2,M2_fid,T2_fid

    def se_analysis(self, start, finish):
        self.data_dict = {} #???? Why am I giving here the second time?
        self.se_max = []
        self.se_T2 = []
        maximum_se = []
        self.echo_time_array  = []
        echo_time = []
        Extrapolation = []
        self.Manual_echo_time = False

        data_fid = self.selected_files_fid[0]
        if self.selected_files_fid_empty!=[] and self.ui.checkBox_subempty.isChecked():
            data_fid_empty = self.selected_files_fid_empty[0]
        else:
            data_fid_empty = self.selected_files_fid_empty

        for i, filename in enumerate(self.selected_files_se):
            # Subtraction, Cut_beginning, Normalize_to_fid, n_start, n_end,  Long_component, end, _, _, _, AdjustZero, zero = self.checkstate_analysis()

            Subtraction, Cut_beginning, Normalize_to_fid, n_start, n_end,  Long_component, end, Apodize, time_a, Adjust, AdjustZero, zero, Smooth, SmoothOrder, SmoothWindow = self.checkstate_analysis()

            data = filename

            if self.selected_files_se_empty != []:
                # I know it can be without !=, but this will help me later

                if len(self.selected_files_se_empty) == len(self.selected_files_se):
                    data_empty = self.selected_files_se_empty[i]
                else:
                    QMessageBox.warning(self, "Different amount of files", f"Please upload the same amount of data and empty files.\nThere are {len(self.selected_files_se)} data files.\nClearing the empty files.", QMessageBox.Ok)
                    self.selected_files_se_empty = []
                    return
            else:
                data_empty = self.selected_files_se_empty
                # = []

            Time_td, Re_td, _, _ = math.nmr_signal_correction(data, data_fid, data_empty, data_fid_empty, Subtraction, Cut_beginning, Normalize_to_fid, n_start, n_end,  Long_component, end, AdjustZero, zero, Smooth, SmoothOrder, SmoothWindow)

            # Subtraction, Cut_beginning, Normalize_to_fid, n_start, n_end,  Long_component, end, Apodize, time_a, Adjust, AdjustZero, zero = self.checkstate_analysis()
            FrequencySE, Real_fftSE = math.freq_domain_correction(Time_td, Re_td, 0, Apodize, time_a, Adjust)
            M2_SE, T2_SE    = math.calculate_M2(Real_fftSE, FrequencySE)

            # Just in case
            self.data_dict[filename] = {
                'Time_td': Time_td,
                'Re_td': Re_td,

            }

            if self.manual_assignment_state == False:
                try:
                    pattern_echo_time = re.compile(r'.*_\s*(\d+)_c\.dat')
                    match = pattern_echo_time.search(filename)
                    file_key = match.group(1)
                    echo_time.append(file_key)
                except:
                    QMessageBox.warning(self, "Couldn't read echo time", f"Couldn't read the echo time automatically.\nPlease insert the values manually in the following window.\nYou can name your files in such a manner to avoid this warning in the future:\nFileName_EchoTimeValue_c.dat", QMessageBox.Ok)
                    self.manual_assignment()

            elif self.manual_assignment_state == True:
                echo_time = self.data_echo_time

            maximum_se.append(max(Re_td))
            self.se_max.append(max(Re_td))
            self.se_T2.append(T2_SE)

        try:
            Extrapolation, fitting_curve, fittingx = math.find_maximum_se(echo_time, maximum_se, start, finish)
        except Exception as e:
            QMessageBox.warning(self, "Couldn't find maximum", f"Couldn't find maximum because {e}.\nChange the fitting parameters and try again", QMessageBox.Ok)
            return

        self.se_max = np.array(self.se_max)
        self.se_T2 = np.array(self.se_T2)
        self.echo_time_array = np.array(echo_time)

        return Extrapolation, echo_time, maximum_se, fittingx, fitting_curve

    def exportSEMaxT2(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontConfirmOverwrite

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Max Amplitudes of SE and T2*",
            "MaxSeT2.txt",
            "Text Files (*.txt);;All Files (*)",
            options=options
        )

        with open(file_path, 'a', encoding='utf-8') as f:
            f.write(f'Max Extrapolation: {self.Extrapolation}\n')
            f.write(f'Echo Time\tT2\tRatio\tEcho Amp\n')
            for et, max, T2 in zip( self.echo_time_array, self.se_max, self.se_T2):
                f.write(f"{et}\t{T2}\t{self.Extrapolation/max}\t{max}\t\n")
            f.write("\n")

    def build_up(self, begin, finish):
        graph_build_fid = self.ui.widgetBUNMRSignal
        graph_build_fft = self.ui.widgetBUSpectra

        Frequency_fid = self.fid_dict['Freq']
        Real_fft_fid = self.fid_dict['Real_fft']
        T2_fid = self.fid_dict['T2']
        M2_fid = self.fid_dict['M2']
        Time_td_fid = self.fid_dict['Time_td_fid']
        Re_td_fid = self.fid_dict['Re_td']

        try:
            function_to_fit = self.ui.comboBox.currentText()
            Time_build_full, Data_build_full, _ = math.build_up_fid(Time_td_fid, Re_td_fid, self.Extrapolation, function_to_fit, begin, finish)

            # PLO
            graph_build_fid.clear()
            graph_build_fid.plot(Time_build_full, Data_build_full, pen=mkPen('b', width=5), symbol = None)
            graph_build_fid.plot(Time_td_fid, Re_td_fid, pen=mkPen('r', width=5), symbol = None)

        except Exception as e:
            QMessageBox.warning(self, "Couldn't build", f"Couldn't build up FID because {e}.\nChange something and try again", QMessageBox.Ok)
            graph_build_fid.clear()
            return

        try:
            time_a = self.ui.doubleSpinBox_5.value()
            Frequency_buildupfid, Real_buildupfid  = math.freq_domain_correction(Time_build_full, Data_build_full, 0, True, time_a, True)

            # PLOT
            graph_build_fft.clear()
            graph_build_fft.plot(Frequency_buildupfid, Real_buildupfid, pen=mkPen('b', width=5), symbol = None)
            graph_build_fft.plot(Frequency_fid, Real_fft_fid, pen=mkPen('r', width=5), symbol = None)

            M2_bu, T2_bu    = math.calculate_M2(Real_buildupfid, Frequency_buildupfid)

            self.ui.textEdit_3.setText(f"FID:\nT₂* = {T2_fid}\nM₂ = {M2_fid}\nBuild up:\nT₂* = {round(T2_bu, 5)}\nM₂ = {round(M2_bu, 5)}")
            return M2_bu, T2_bu

        except Exception as e:
            QMessageBox.warning(self, "Couldn't build", f"Couldn't build up FFT because {e}.\nChange something and try again", QMessageBox.Ok)
            return

    def extrapolate(self):
        graph_se = self.ui.widgetSEMax
        try:
            if self.state == 'se' or self.state == 'single_se':
                if len(self.selected_files_se) > 4:
                    self.ui.doubleSpinBox_6.setEnabled(True)
                    self.ui.doubleSpinBox_7.setEnabled(True)
                    start   = self.ui.doubleSpinBox_7.value()
                    end     = self.ui.doubleSpinBox_6.value()
                    self.Extrapolation, echo_time, maximum_se, fittingx, fittingy = self.se_analysis(start, end)
                    echo_time = np.array(echo_time, dtype=float)
                    maximum_se = np.array(maximum_se, dtype=float)
                else:
                    self.state = 'single_se'
                    Time_td, Re_td,_, _,Frequency, Real_fft,_, _,M2,T2,_,_ = self.basic_data_analysis()
                    self.Extrapolation = max(Re_td)
                    self.data_dict['Single SE'] = {
                        'Time_td': Time_td,
                        'Re_td': Re_td,
                        'Freq': Frequency,
                        'Real_fft': Real_fft,
                        'M2': M2,
                        'T2': T2
                    }
            elif self.state == 'mse':
                Time_td, Re_td,_, _,Frequency, Real_fft,_, _,M2,T2,_,_ = self.basic_data_analysis()
                self.Extrapolation = max(Re_td)

                self.data_dict['MSE'] = {
                    'Time_td': Time_td,
                    'Re_td': Re_td,
                    'Freq': Frequency,
                    'Real_fft': Real_fft,
                    'M2': M2,
                    'T2': T2
                }

            show_maximum = round(self.Extrapolation, 2)
            match self.state:
                case "se":
                    graph_se.clear()
                    graph_se.show()
                    graph_se.plot([0], [self.Extrapolation], pen='r', symbolPen=None, symbol='o', symbolBrush='r')
                    graph_se.plot(fittingx, fittingy, pen=mkPen('k', width=3, style=Qt.DashLine), symbol = None)
                    graph_se.plot(echo_time, maximum_se,pen=None, symbolPen=None, symbol='o', symbolBrush='b')
                    self.ui.textEdit_2.setText(f"Max = {show_maximum}")
                case "mse":
                    self.Extrapolation = self.Extrapolation/2 # 4 phases, derived by 2, logical, yeah 
                    show_maximum = round(self.Extrapolation, 2)
                    graph_se.clear()
                    graph_se.hide()
                    self.ui.textEdit_2.setText(f"Max = {show_maximum}")
                    self.ui.doubleSpinBox_6.setEnabled(False)
                    self.ui.doubleSpinBox_7.setEnabled(False)
                case "single_se":
                    graph_se.clear()
                    graph_se.hide()
                    self.ui.textEdit_2.setText(f"Max = {show_maximum}")
                    self.ui.doubleSpinBox_6.setEnabled(False)
                    self.ui.doubleSpinBox_7.setEnabled(False)

        except Exception as e:
            QMessageBox.warning(self, "Couldn't extrapolate", f"Couldn't extrapolate because {e}.\nChange something and try again", QMessageBox.Ok)
            graph_se.clear()
            return

    def compare(self):

        graph_nmr = self.ui.widgetOriginalNMRSignal
        graph_fft = self.ui.widgetOriginalSpectra

        try:
            self.ui.groupBox.setEnabled(False)
            self.fid_dict = {"Time_td_fid": [], "Re_td": [], "Freq": [], "Real_fft": [], 'M2': [], 'T2': []}

            Time_td, Re_td,Time_td_fid, Re_td_fid,Frequency, Real_fft,Frequency_fid, Real_fft_fid,M2,T2,M2_fid,T2_fid = self.basic_data_analysis()

            # Plot and write text for tab 1
            graph_nmr.clear()
            graph_nmr.plot(Time_td, Re_td, pen=mkPen('r', width=5), symbol = None)
            graph_nmr.plot(Time_td_fid, Re_td_fid, pen=mkPen('b', width=5), symbol = None)

            graph_fft.clear()
            graph_fft.plot(Frequency, Real_fft, pen=mkPen('r', width=5), symbol = None)
            graph_fft.plot(Frequency_fid, Real_fft_fid, pen=mkPen('b', width=5), symbol = None)

            self.ui.label.show()

            self.ui.textEdit.setText(f"FID:\nT₂* = {T2_fid}\nM₂ = {M2_fid}\nData:\nT₂* = {T2}\nM₂ = {M2}")

        except Exception as e:
            QMessageBox.warning(self, "Couldn't analyse", f"Couldn't analyse the data because {e}.\nI have cleared the files, try again", QMessageBox.Ok)
            self.start()
            return

    def run(self):
        self.compare()
        self.extrapolate()

        Time_range = self.fid_dict['Time_td_fid']
        minimum = min(Time_range)
        self.ui.doubleSpinBox_finish.setMinimum(minimum)
        # Does it work? It doesn't work TODO

        begin = self.ui.doubleSpinBox_begin.value()
        finish = self.ui.doubleSpinBox_finish.value()
        _, _ = self.build_up(begin, finish)
        self.ui.groupBox_6.setEnabled(True)

    def rebuild(self):
        self.extrapolate()

        begin = self.ui.doubleSpinBox_begin.value()
        finish = self.ui.doubleSpinBox_finish.value()
        _, _ = self.build_up(begin, finish)

    def manual_assignment(self):
        data = {
            'se' : self.selected_files_se,
        }

        array = self.data_echo_time
        dialog = EchoTimeRangeDialog(self, data, array)
        dialog.exec()
        echo_time, self.manual_assignment_state  = dialog.form_array_echo_time()
        self.data_echo_time = echo_time

        return echo_time

    def time_analysis(self):
        M2 = []
        T2 = []
        range_bf = []

        Time = self.fid_dict["Time_td_fid"]
        minimum = min(Time)
        maximum = self.ui.spinBoxTimeRange.value()
        start_range = np.arange(minimum, maximum, 0.5)
        finish_range = np.arange(minimum+5, maximum+5, 0.5)
        for begin in start_range:
            for finish in finish_range:
                if finish < begin + 3:
                    M2_bu = 0
                    T2_bu = 0
                else:
                    try:
                        M2_bu, T2_bu = self.build_up(begin, finish)
                    except:
                        M2_bu = 0
                        T2_bu = 0

                M2.append(M2_bu)
                T2.append(T2_bu)
                range_bf.append([begin, finish])


        T2_reshaped = np.array(T2).reshape(len(start_range), len(finish_range))
        M2_reshaped = np.array(M2).reshape(len(start_range), len(finish_range))

        self.plot_window = MatplotlibWindow(self)
        self.plot_window.plot_data(start_range, finish_range, T2_reshaped, M2_reshaped)
        self.plot_window.setWindowTitle("2D Visualization")
        self.plot_window.resize(800, 600)
        self.plot_window.show()

    def save_statistics(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontConfirmOverwrite

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Frequency and T2 Data",
            "statistical_analysis.txt",
            "Text Files (*.txt);;All Files (*)",
            options=options
        )

        self.popup = SavingPopup()
        self.popup.show()
        QApplication.processEvents()


        if not file_path:
            print("Saving cancelled.")
            return

        for comboBoxIndex in range(self.ui.comboBox.count()):
            self.ui.comboBox.setCurrentIndex(comboBoxIndex)

            M2 = []
            T2 = []
            range_bf = []

            Time = self.fid_dict["Time_td_fid"]
            minimum = min(Time)
            maximum = self.ui.spinBoxTimeRange.value()

            start_range = np.arange(minimum, maximum, 0.5)
            finish_range = np.arange(minimum+5, maximum+5, 0.5)
            for begin in start_range:
                for finish in finish_range:
                    if finish < begin + 3:
                        M2_bu = 0
                        T2_bu = 0
                    else:
                        try:
                            M2_bu, T2_bu = self.build_up(begin, finish)
                        except:
                            M2_bu = 0
                            T2_bu = 0

                    M2.append(M2_bu)
                    T2.append(T2_bu)
                    range_bf.append([begin, finish])

            self.write_freqencies_T2(range_bf, T2, file_path)

        self.popup.close()

    def write_freqencies_T2(self, ranges, T2s, file_path):
        with open(file_path, 'a', encoding='utf-8') as f:
            f.write(f'Name: {self.selected_files_fid}\n')
            f.write(f'Mode: {self.state}\n')
            f.write(f'Function: {self.ui.comboBox.currentText()}\n\n')
            for range_bf, T2 in zip(ranges, T2s):
                f.write(f"{range_bf[0]}\t{range_bf[1]}\t{T2}\t\n")
            f.write("\n")

class MatplotlibWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.canvas = FigureCanvas(plt.figure())

        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(100)
        self.slider.setValue(100)
        self.slider.valueChanged.connect(self.update_plot)

        self.slider_label = QLabel("T2 Threshold: 100%")

        self.rainbow_checkbox = QCheckBox("Rainbow")
        self.rainbow_checkbox.stateChanged.connect(self.update_plot)

        self.save_button = QPushButton("Save Data")
        self.save_button.clicked.connect(self.save_data)

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)

        slider_layout = QHBoxLayout()
        slider_layout.addWidget(self.slider_label)
        slider_layout.addWidget(self.slider)
        slider_layout.addWidget(self.rainbow_checkbox)
        slider_layout.addWidget(self.save_button)

        layout.addLayout(slider_layout)
        self.setLayout(layout)

        self.start_range = None
        self.finish_range = None
        self.T2 = None
        self.M2 = None

    def plot_data(self, start_range, finish_range, T2, M2):
        self.start_range = start_range
        self.finish_range = finish_range
        self.T2 = T2
        self.M2 = M2
        self.update_plot()

    def update_plot(self):
        if self.T2 is None and self.M2 is not None:
            return

        threshold = self.slider.value() / 100.0 * np.max(self.T2)
        self.slider_label.setText(f"T₂ Threshold: {threshold:.2f}")
        T2_masked = np.ma.masked_where((self.T2 < 5) | (self.T2 > threshold), self.T2)

        cmap = "rainbow" if self.rainbow_checkbox.isChecked() else "RdYlGn"

        self.canvas.figure.clf()

        gs = self.canvas.figure.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.1)
        ax_main = self.canvas.figure.add_subplot(gs[0])
        ax_hist = self.canvas.figure.add_subplot(gs[1])

        c = ax_main.pcolormesh(self.finish_range, self.start_range, T2_masked, shading="auto", cmap=cmap)
        self.canvas.figure.colorbar(c, ax=ax_main, label="T₂ Value")
        ax_main.set_xlabel("Finish Range")
        ax_main.set_ylabel("Start Range")
        ax_main.set_xlim(self.finish_range[0], self.finish_range[-1])
        ax_main.set_ylim(self.start_range[0], self.start_range[-1])

        ax_hist.hist(T2_masked.flatten(), bins=20, color="blue")
        ax_hist.set_ylabel("Count")
        ax_hist.yaxis.set_label_position("right")
        ax_hist.yaxis.tick_right()
        ax_hist.set_xlabel("T₂ value")
        ax_main.set_title("T₂ Distribution")
        ax_hist.set_xlim(np.min(T2_masked), np.max(T2_masked))

        # self.canvas.figure.tight_layout()
        self.canvas.draw()

    def save_data(self):
        if self.start_range is None or self.finish_range is None or self.T2 is None or self.M2 is None:
            return

        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Data", "", "CSV Files (*.csv);;Text Files (*.txt)", options=options)
        if not file_path:
            return

        data = np.column_stack((self.start_range.repeat(len(self.finish_range)),
                                np.tile(self.finish_range, len(self.start_range)),
                                self.M2.flatten(), self.T2.flatten()))
        header = "Start Range,Finish Range,M2Value,T2 Value"
        np.savetxt(file_path, data, delimiter=",", header=header, comments="")

    def closeEvent(self, event):
        self.deleteLater()
        event.accept()

class OpenFilesDialog(QFileDialog):
    def __init__(self, parent=None):

        super().__init__(parent)

        if Main.State_multiple_files:
            self.setFileMode(QFileDialog.ExistingFiles)  # Allow selecting multiple files
        else:
            pass

        self.setNameFilter(str("Data (*.dat *.txt *.csv)"))

        initial_directory_file = "selected_folder.txt"

        try:
            with open(initial_directory_file, 'r') as file:
                directory = file.read().strip()
        except Exception as e:
            print(f"Couldn't read the initial directory: {e}")
            directory = os.path.dirname(sys.argv[0])

        self.setDirectory(directory)

        self.selected_files = []  # dictionary to store selected file paths

    def on_file_selected(self):
        options = QFileDialog.Options()
        files, _ = QFileDialog.getOpenFileNames(self, "Load Files", "", "Data Files (*.dat *.txt *.csv)", options=options)
        if files:
            self.selected_files.extend(files)

class EchoTimeRangeDialog(QDialog):
    def __init__(self, parent=None, data=None, array = None):
        super().__init__(parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.ui.pushButton.setEnabled(False)

        if array is not None:
            self.echo_time = array

        if (data is not None) and (len(data['se'])>1):
            self.ui.pushButton.setEnabled(True)
            self.populate_table(data)

        self.ui.pushButton.clicked.connect(self.form_array_echo_time)

    def populate_table(self, selected_data):
        table = self.ui.tableWidget

        # idiotic way, but meh
        total_length = 0
        for key in selected_data.values():
            total_length += len(key)

        files_array = list(selected_data.values())
        files_array = [item for item in files_array if item]

        table.setRowCount(total_length)
        row = 0
        for value in files_array:
            for file in value:
                filename = os.path.basename(file)
                item = QTableWidgetItem(str(filename))
                table.setItem(row, 0, item)

                if self.echo_time == []:
                    try:
                        pattern_echo_time = re.compile(r'.*_\s*(\d+)_c\.dat')
                        match = pattern_echo_time.search(filename)
                        file_key = match.group(1)
                        item = QTableWidgetItem(str(file_key))
                    except:
                        item = QTableWidgetItem('')
                else:
                    item = QTableWidgetItem(self.echo_time[row])
                table.setItem(row, 1, item)

                row += 1
        self.ui.tableWidget.resizeColumnsToContents()
        self.adjust_window_size()

    def adjust_window_size(self):
        table_width = self.ui.tableWidget.horizontalHeader().length()
        self.resize(table_width + 100, 300)

    def form_array_echo_time(self):
        echo_time = []
        table = self.ui.tableWidget

        for row in range(table.rowCount()):
            item = table.item(row, 1)

            if item is not None and item.text() != '' and item.text().isdigit():
                echo_time.append((item.text()))
                manual_state = True
            else:
                QMessageBox.warning(self, "Couldn't write echo time", f"Couldn't write the entered echo time in row {row+1}.\nPlease insert the numerical value.", QMessageBox.Ok)
                manual_state = False
                return

        self.echo_time = echo_time

        return echo_time, manual_state

class SavingPopup(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Please Wait")
        self.setFixedSize(250, 100)

        layout = QVBoxLayout()

        # Simple, readable label
        label = QLabel("Wait for the results to save", self)
        label.setStyleSheet("font-size: 14px; padding: 10px; color: black;")

        layout.addWidget(label)
        self.setLayout(layout)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = Main()
    widget.show()
    sys.exit(app.exec())
