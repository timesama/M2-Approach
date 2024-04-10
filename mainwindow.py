# This Python file uses the following encoding: utf-8
import sys, os, re
import logging
from PySide6.QtWidgets import QApplication, QMainWindow, QFileDialog, QTableWidgetItem, QDialog, QMessageBox, QPushButton
from PySide6.QtCore import QCoreApplication, Signal, SIGNAL
from PySide6.QtGui import QColor
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import pyqtgraph as pg
from pyqtgraph.exporters import ImageExporter
import matplotlib.pyplot as plt
from ui_Form import Ui_NMR
from ui_Notification import Ui_Note
from ui_Error import Ui_Error
from ui_PhasingManual import Ui_Phasing as Ui_PhasingManual

pg.CONFIG_OPTIONS['background'] = 'w'
pg.CONFIG_OPTIONS['foreground'] = 'k'

# Global 
Frequency = []
Re_spectra = []
Im_spectra = []

# Configure logging settings
logging.basicConfig(filename='app.log', level=logging.DEBUG)
logger = logging.getLogger()

class LoggedButton(QPushButton):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def mousePressEvent(self, event):
        logger.debug(f'Button pressed: {self.text()}')
        super().mousePressEvent(event)

# Math procedures
def analysis_time_domain(file_path):
    # 1. Read data
    Time, Real, Imag = read_data(file_path)
    # 2. Crop time below zero
    T_cr, R_cr, I_cr = crop_time_zero(Time, Real, Imag)

    # 3. Phase the data
    R_ph, I_ph = time_domain_phase(R_cr, I_cr)

    # 4. Adjust Frequency
    # 4.1 Calculate Freq
    Frequency = calculate_frequency_scale(T_cr)
    # 4.2 Shift Freq
    R_sh, I_sh = adjust_frequency(Frequency, R_ph, I_ph)

    return T_cr, R_sh, I_sh

def reference_long_component(Time, Component, Amplitude_gly, coeff):
    # 2. Normalize (reference) components to Amplitude of the reference
    Component_n = Component/Amplitude_gly

    # 3. Cut the ranges for fitting
    minimum = find_nearest(Time, 50)
    maximum = find_nearest(Time, 250)

    Time_range = Time[minimum:maximum]
    Component_n_range = Component_n[minimum:maximum]

    # 7. Fit data to exponential decay
    popt, _      = curve_fit(decaying_exponential, Time_range, Component_n_range, p0=coeff)
    
    # 8. Set the ranges for subtraction
    Time_cropped  = Time[0:maximum]
    Component_c   = Component_n[0:maximum]

    # 9. Calculate the curves fitted to data within the desired range
    Component_f = decaying_exponential(Time_cropped, *popt)

    # 10. Subtract
    Component_sub = Component_c - Component_f

    return Component_sub, Time_cropped

def long_component(Time_s, Time_r, Re_s, Re_r, Im_s, Im_r):
    # r stands for reference, s stands for sample
    Amp_r = calculate_amplitude(Re_r, Im_r)
    Amp_s = calculate_amplitude(Re_s, Im_s)

    # 1. Crop the arrays together (they should be of the same length, but I know, I know...)
    if len(Time_s) > len(Time_r):
        Time  =   Time_s[:len(Time_r)]
        Amp_s   =   Amp_s[:len(Time_r)]
        Re_s    =   Re_s[:len(Time_r)]
        Im_s    =   Im_s[:len(Time_r)]
    else:  
        Time  =   Time_r[:len(Time_s)]
        Amp_r   =   Amp_r[:len(Time_s)]

    coeff_re = [0.9, 400, 0.1]
    coeff_im = [1, 400, 0]

    Real_subtracted, Time_cropped   = reference_long_component(Time, Re_s, Amp_r, coeff_re)
    Im_subtracted, _     = reference_long_component(Time, Im_s, Amp_r, coeff_im)
    
    # 11. Normalize
    Re_n, Im_n = normalize(Real_subtracted, Im_subtracted)

    return Time_cropped, Re_n, Im_n

def final_analysis_time_domain(Time, Real, Imaginary):
    # 5. Apodize the time-domain
    Re_ap, Im_ap = apodization(Time, Real, Imaginary)
    
    # 6. Add zeros
    Tim, Fid = add_zeros(Time, Re_ap, Im_ap, 16383)

    #stophere
    return Tim, Fid

def frequency_domain_analysis(FFT, Frequency):

    # 8. Simple baseline
    _, Re, _ = simple_baseline_correction(FFT)

    # 9. Apodization
    Real_apod = calculate_apodization(Re, Frequency)

    # 10. M2 & T2
    M2, T2 = calculate_M2(Real_apod, Frequency)

    return M2, T2

def read_data(file_path):
    data = np.loadtxt(file_path)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    return x, y, z

def crop_time_zero(Time, Real, Imaginary):
    if Time[0] < 0:
        Time_start = 0
        Time_crop_idx = np.where(Time >= Time_start)[0][0]
        Time_cropped = Time[Time_crop_idx:]
        Real_cropped = Real[Time_crop_idx:]
        Imaginary_cropped = Imaginary[Time_crop_idx:]
        return Time_cropped, Real_cropped, Imaginary_cropped
    else:
        return Time, Real, Imaginary

def time_domain_phase(Real, Imaginary):
    delta = np.zeros(360)
    
    for phi in range(360):
        Re_phased = Real * np.cos(np.deg2rad(phi)) - Imaginary * np.sin(np.deg2rad(phi))
        Im_phased = Real * np.sin(np.deg2rad(phi)) + Imaginary * np.cos(np.deg2rad(phi))
        Magnitude_phased = calculate_amplitude(Re_phased, Im_phased)
        
        Re_cut = Re_phased[:5]
        Ma_cut = Magnitude_phased[:5]
        
        delta[phi] = np.mean(Ma_cut - Re_cut)
    
    idx = np.argmin(delta)
    #print(idx)

    Re = Real * np.cos(np.deg2rad(idx)) - Imaginary * np.sin(np.deg2rad(idx))
    Im = Real * np.sin(np.deg2rad(idx)) + Imaginary * np.cos(np.deg2rad(idx))

    return Re, Im

def adjust_frequency(Frequency, Re, Im):
    # Create complex FID
    Fid_unshifted = np.array(Re + 1j * Im)

    # FFT
    FFT = np.fft.fftshift(np.fft.fft(Fid_unshifted))

    # Check the length of FFT and Frequency (it is always the same, this is just in case)
    if len(Frequency) != len(FFT):
        Frequency = np.linspace(Frequency[0], Frequency[-1], len(FFT))

    # Find index of max spectrum (amplitude)
    index_max = np.argmax(FFT)

    # Find index of zero (frequency)
    index_zero = find_nearest(Frequency, 0)

    # Find difference
    delta_index = index_max - index_zero

    # Shift the spectra (amplitude) by the difference in indices
    FFT_shifted = np.concatenate((FFT[delta_index:], FFT[:delta_index]))

    # iFFT
    Fid_shifted = np.fft.ifft(np.fft.fftshift(FFT_shifted))

    # Define Real, Imaginary and Amplitude
    Re_shifted = np.real(Fid_shifted)
    Im_shifted = np.imag(Fid_shifted)

    return Re_shifted, Im_shifted

def normalize(Real, Imaginary):
    Amplitude = np.sqrt(Real ** 2 + Imaginary ** 2)
    Amplitude_max = np.max(Amplitude)
    Amp = Amplitude/Amplitude_max
    Re = Real/Amplitude_max
    Im = Imaginary/Amplitude_max
    return Re, Im

def apodization(Time, Real, Imaginary):
    Amplitude = calculate_amplitude(Real, Imaginary)
    coeffs = np.polyfit(Time, Amplitude, 1)  # Fit an exponential decay function
    c = np.polyval(coeffs, Time)
    d = np.argmin(np.abs(c - 1e-5))
    sigma = Time[d]
    if sigma == 0:
        sigma = 1000
    apodization_function = np.exp(-(Time / sigma) ** 4)
    Re_ap = Real * apodization_function
    Im_ap = Imaginary * apodization_function
    return Re_ap, Im_ap

def add_zeros(Time, Real, Imaginary, number_of_points):
    length_diff = number_of_points - len(Time)
    amount_to_add = np.zeros(length_diff+1)

    Re_zero = np.concatenate((Real, amount_to_add))
    Im_zero = np.concatenate((Imaginary, amount_to_add))

    dt = Time[1] - Time[0]
    Time_to_add = Time[-1] + np.arange(1, length_diff + 1) * dt

    Time = np.concatenate((Time, Time_to_add))
    Fid = np.array(Re_zero + 1j * Im_zero)
    Fid = Fid[:-1]

    return Time, Fid

def simple_baseline_correction(FFT):
    twentyperc = int(round(len(FFT) * 0.02))
    Baseline = np.mean(np.real(FFT[:twentyperc]))
    FFT_corrected = FFT - Baseline
    Re = np.real(FFT_corrected)
    Im = np.imag(FFT_corrected)
    Amp = calculate_amplitude(Re, Im)
    return Amp, Re, Im

def calculate_apodization(Real, Freq):
    # Find sigma at 2% from the max amplitude of the spectra
    Maximum = np.max(np.abs(Real))
    idx_max = np.argmax(np.abs(Real))
    ten_percent = Maximum * 0.02

    b = np.argmin(np.abs(Real[idx_max:] - ten_percent))
    sigma_ap = Freq[idx_max + b]

    apodization_function_s = np.exp(-(Freq / sigma_ap) ** 6)

    Real_apod = Real * apodization_function_s
    
    return Real_apod

def calculate_amplitude(Real, Imaginary):
    # NoClass
    Amp = np.sqrt(Real ** 2 + Imaginary ** 2)
    return Amp

def calculate_frequency_scale(Time):
    numberp = len(Time)

    dt = Time[1] - Time[0]
    f_range = 1 / dt
    f_nyquist = f_range / 2
    df = 2 * (f_nyquist / numberp)
    Freq = np.arange(-f_nyquist, f_nyquist + df, df)
    Freq = Freq[:-1]

    return Freq

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def calculate_M2(FFT_real, Frequency):
    # NoClass
    # Take the integral of the REAL PART OF FFT by counts
    Integral = np.trapz(np.real(FFT_real))
    
    # Normalize FFT to the Integral value
    Fur_normalized = np.real(FFT_real) / Integral
    
    # Calculate the integral of normalized FFT to receive 1
    Integral_one = np.trapz(Fur_normalized)
    
    # Multiplication (the power ^n will give the nth moment (here it is n=2)
    Multiplication = (Frequency ** 2) * Fur_normalized
    
    # Calculate the integral of multiplication - the nth moment
    # The (2pi)^2 are the units to transform from rad/sec to Hz
    # ppbly it should be (2pi)^n for generalized moment calculation
    M2 = (np.trapz(Multiplication)) * 4 * np.pi ** 2
    
    # Check the validity
    if np.abs(np.mean(Multiplication[0:10])) > 10 ** (-6):
        print('Apodization is wrong!')

    if M2 < 0:
        M2 = 0
        T2 = 0
    else:
        T2 = np.sqrt(2/M2)
    
    return M2, T2

def calculate_SFC(Amplitude):
    S = np.mean(Amplitude[1:4])
    L = np.mean(Amplitude[50:70])
    SFC = (S-L)/S
    return SFC

def calculate_DQ_intensity(Time, Amplitude):
    idx_time = np.argmin(np.abs(Time - 4))
    DQ = np.mean(Amplitude[:idx_time])
    return DQ

def decaying_exponential(x, a, b, c):
    return a * np.exp(-x/b) + c

def decaying_2exponential(x, a1, b1, a2, b2, c):
    return a1 * np.exp(-x/b1) + a2 * np.exp(-x/b2) + c

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_NMR()
        self.ui.setupUi(self)

        # Set window geometry
        screen = QApplication.primaryScreen()

        if screen:
            available_geometry = screen.availableGeometry()
            left = top = 10
            width = available_geometry.width() - left
            height = available_geometry.height() - 3 * top
            self.setMaximumSize(width, height)
            self.showMaximized()
        

        self.selected_files = []
        self.selected_folders = []
        self.selected_files_gly = []
        self.selected_DQfiles = []
        self.dq_t2 = {}
        self.t1_dictionary ={}

        # Connect buttons to their respective slots
        self.ui.btn_SelectFiles.clicked.connect(self.clear_list)
        self.ui.btn_SelectFiles.clicked.connect(self.open_select_dialog)
        self.ui.btn_Add.clicked.connect(self.open_select_dialog)
        self.ui.btn_Start.clicked.connect(self.analysis)
        self.ui.btn_Save.clicked.connect(self.save_data)
        self.ui.btn_Load.clicked.connect(self.load_data)
        self.ui.btn_Phasing.clicked.connect(self.open_phasing_manual)
        self.ui.radioButton_Log.clicked.connect(self.plot_fit)
        self.ui.tabWidget.currentChanged.connect(self.groupBox_status)
        self.ui.btn_SelectFilesDQ.clicked.connect(self.open_select_DQdialog)
        self.ui.btn_ClearTable_2.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable.clicked.connect(self.clear_list)
        self.ui.btn_Launch.clicked.connect(self.launch)
        self.ui.btn_SelectFolders_T1.clicked.connect(self.open_folder_dialog)
        self.ui.btn_Plot1.clicked.connect(self.state_1exp)
        self.ui.btn_Plot2.clicked.connect(self.state_2exp)
        #self.ui.checkBox.clicked.connect(self.disable_buttons2)

        # Graph setup
        self.setup_graph(self.ui.FFTWidget, "Frequency, MHz", "Amplitude, a.u", "FFT")
        self.setup_graph(self.ui.FidWidget, "Time, ms", "Amplitude", "FID")
        self.setup_graph(self.ui.SEWidget, "Temperature, °C", "Choose", "")
        self.setup_graph(self.ui.DQ_Widget_1, "DQ Filtering Time", "T₂*", "")
        self.setup_graph(self.ui.DQ_Widget_2, "", "Norm. DQ Intensity", "")
        self.setup_graph(self.ui.DQ_Widget_3, "DQ Filtering Time", "T₂*", "")
        self.setup_graph(self.ui.DQ_Widget_4, "", "Norm. DQ Intensity", "")
        self.setup_graph(self.ui.DQ_Widget_5, "Name", "Center", "")
        self.setup_graph(self.ui.T1_Widget_1, "Time, μs", "Signal", "")
        self.setup_graph(self.ui.T1_Widget_2, "Temperature, °C", "T₁, μs", "")

        # Table setup
        self.setup_table(self.ui.table_SE)
        self.setup_table(self.ui.table_DQ)

        # Connect table signals to slots
        self.ui.table_DQ.currentItemChanged.connect(self.update_dq_graphs)
        self.ui.table_SE.currentItemChanged.connect(self.update_yaxis)
        self.ui.radioButton_2.clicked.connect(self.calculate_T1)
        self.ui.radioButton_3.clicked.connect(self.calculate_T1)


        # Connect combobox signals to slots
        self.ui.comboBox.currentIndexChanged.connect(self.update_yaxis)
        self.ui.comboBox_3.currentIndexChanged.connect(self.update_yaxis)
        self.ui.comboBox_2.currentIndexChanged.connect(self.plot_fit)
        self.ui.comboBox_6.currentIndexChanged.connect(self.calculate_T1)

        # Connect change events
        self.ui.dq_min.valueChanged.connect(self.update_dq_graphs)
        self.ui.dq_max.valueChanged.connect(self.update_dq_graphs)

        self.ui.dq_min_2.valueChanged.connect(self.update_DQ_comparison_plot)
        self.ui.dq_max_2.valueChanged.connect(self.update_DQ_comparison_plot)

        self.ui.radioButton_Log_2.clicked.connect(self.update_DQ_comparison_plot) # this is a bad coding
        self.ui.comboBox_5.currentIndexChanged.connect(self.update_DQ_comparison_plot)
        self.ui.comboBox_4.currentIndexChanged.connect(self.update_file)

        # Disable buttons initially
        self.disable_buttons()

    def update_file(self):
        i = self.ui.comboBox_4.currentIndex() + 1
        current_tab_index =  self.ui.tabWidget.currentIndex()
        try:
            file_path = self.selected_files[i-1]
        except:
            return
        
        
        if self.ui.checkBox_2.isChecked():
            try:
                file_path_gly = self.selected_files_gly[i-1]
                self.process_file_data(file_path, file_path_gly, current_tab_index, i)
            except:
                return
        else:
            self.process_file_data(file_path, [], current_tab_index, i)

            
        # Update general figures
        if current_tab_index == 1:
            self.highlight_row(self.ui.table_DQ, i) 
            self.update_dq_graphs()
        elif current_tab_index == 0:
            self.highlight_row(self.ui.table_SE, i)
            self.update_yaxis()

        self.ui.btn_Phasing.setEnabled(True)

        #TODO sometime I should add the highlight of the certain point on graph, but I am too lazy
            
    def clear_list(self):
        self.selected_files = []
        self.selected_files_gly = []
        self.selected_DQfiles = []
        self.selected_folders = []
        self.ui.table_SE.setRowCount(0)
        self.ui.table_DQ.setRowCount(0)
        self.ui.table_DQ_2.setRowCount(0)
        self.ui.table_T1.setRowCount(0)

        self.ui.T1_Widget_1.clear()
        self.ui.T1_Widget_2.clear()
        self.ui.comboBox_6.currentIndexChanged.disconnect(self.calculate_T1)
        while self.ui.comboBox_6.count()>0:          
            self.ui.comboBox_6.removeItem(0)
        self.ui.comboBox_6.currentIndexChanged.connect(self.calculate_T1)

    def terminate(self):
        self.disable_buttons()
        self.selected_files = []
        self.selected_files_gly = []
        self.selected_DQfiles = []
        self.ui.table_SE.setRowCount(0)
        self.ui.table_DQ.setRowCount(0)
        self.ui.table_DQ_2.setRowCount(0)
        self.ui.DQ_Widget_1.clear()
        self.ui.DQ_Widget_2.clear()
        self.ui.DQ_Widget_3.clear()
        self.ui.DQ_Widget_4.clear()
        self.ui.DQ_Widget_5.clear()
        self.ui.SEWidget.clear()
        self.ui.FFTWidget.clear()
        self.ui.FidWidget.clear()
            
    def highlight_row(self, table, row_selected):

        #for row in range(table.rowCount()):
        #table.selectRow(row_selected-1)


        for col in range(table.columnCount()):
            for row in range(table.rowCount()):
                item = table.item(row, col)
                if item is not None:
                    item.setBackground(QColor(255, 255, 255))
            item_selected = table.item(row_selected-1, col)
            if item_selected is not None:
                item_selected.setBackground(QColor(255, 255, 0))

        # self.ui.table_SE.selectRow(5)
        # self.ui.table_SE.currentRow()
        #TODO normally
  
    def setup_graph(self, graph_widget, xlabel="", ylabel="", title=""):
        graph_widget.getAxis('left').setLabel(ylabel)
        graph_widget.getAxis('bottom').setLabel(xlabel)
        graph_widget.setTitle(title)

    def setup_table(self, table_widget):
        column_widths = [70] * 4
        for i, width in enumerate(column_widths):
            table_widget.setColumnWidth(i, width)

    def open_select_DQdialog(self):
        dlg = OpenFilesDialog(self)
        if dlg.exec():
            DQfileNames = dlg.selectedFiles()
            self.selected_DQfiles.extend(DQfileNames)
            self.update_DQ_comparison()
       
    def open_select_dialog(self):
        dlg = OpenFilesDialog(self)
        if dlg.exec():
            fileNames = dlg.selectedFiles()
            self.selected_files.extend(fileNames)
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)

    def open_folder_dialog(self):
        options = QFileDialog.Options()
        initial_directory = "C:/Mega/NMR/003_Temperature"
        
        while True:
            folder_path = QFileDialog.getExistingDirectory(self, "Select Folder", initial_directory, options=options)
            if folder_path:
                self.selected_folders.append(folder_path)
                initial_directory = os.path.dirname(folder_path)
            else:
                break 
        # Clear Combobox
        while self.ui.comboBox_6.count()>0:
            self.ui.comboBox_6.removeItem(0)
        self.update_T1_table()
           
    def open_select_dialog_glycerol(self):
        dlg = OpenFilesDialog(self)
        if dlg.exec():
            fileNames_gly = dlg.selectedFiles()
            self.selected_files_gly.extend(fileNames_gly)
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)

    def open_phasing_manual(self):
        self.phasing_manual_window = PhasingManual()
        self.phasing_manual_window.read_data()
        self.phasing_manual_window.show()

        self.phasing_manual_window.closed.connect(self.after_phasing)

    def groupBox_status(self):
        current_tab_index =  self.ui.tabWidget.currentIndex()
        if current_tab_index == 2 or current_tab_index == 3:
            self.ui.groupBox.setEnabled(False)
        else:
            self.ui.groupBox.setEnabled(True)

    def disable_buttons(self):
        self.ui.btn_Start.setEnabled(False)
        self.ui.btn_Save.setEnabled(False)
        self.ui.btn_Phasing.setEnabled(False)
        self.ui.dq_min.setEnabled(False)
        self.ui.dq_max.setEnabled(False)
        self.ui.comboBox.setEnabled(False)
        self.ui.comboBox_2.setEnabled(False)
        self.ui.radioButton_Log.setEnabled(False)
        #self.ui.checkBox.setEnabled(False)
        self.ui.btn_Add.setEnabled(False)
        self.ui.radioButton_Log_2.setEnabled(False)
        self.ui.comboBox_5.setEnabled(False)
        self.ui.dq_min_2.setEnabled(False)
        self.ui.dq_max_2.setEnabled(False)
        self.ui.btn_Launch.setEnabled(False)
        self.ui.btn_Plot1.setEnabled(False)
        self.ui.btn_Plot2.setEnabled(False)

    def enable_buttons(self):
        self.ui.btn_SelectFiles.setEnabled(True)
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Save.setEnabled(True)
        self.ui.radioButton.setEnabled(True)
        self.ui.btn_Load.setEnabled(True)
        self.ui.dq_min.setEnabled(True)
        self.ui.dq_max.setEnabled(True)
        self.ui.comboBox.setEnabled(True)
        self.ui.comboBox_2.setEnabled(True)
        #self.ui.checkBox.setEnabled(True)
        self.ui.comboBox_4.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)
    
    def state_1exp(self):
        self.State_exp1 = True
        self.plot_t1_temperature()
    
    def state_2exp(self):
        self.State_exp1 = False
        self.plot_t1_temperature()

    # All these functions refer to the general analysis where FFT and FID are produced
    def analysis(self):
        # Clear Combobox
        while self.ui.comboBox_4.count()>0:
            self.ui.comboBox_4.removeItem(0)

        # Clear Graphs
        self.ui.SEWidget.clear()
        self.ui.DQ_Widget_1.clear()
        self.ui.DQ_Widget_2.clear()
        self.ui.DQ_Widget_3.clear()
        self.ui.DQ_Widget_4.clear()
        self.ui.comboBox.setCurrentIndex(-1)
        self.ui.comboBox_2.setCurrentIndex(-1)
        self.ui.textEdit_4.setText("")

        current_tab_index = self.ui.tabWidget.currentIndex()

        legend = pg.LegendItem(offset=(300, 10), parent=self.ui.FidWidget.graphicsItem())

        if not self.load_data_and_check_validity((self.selected_files[0]), current_tab_index):
            print('This is not statement')
            return
        
        if self.ui.checkBox_2.isChecked():
                self.open_select_dialog_glycerol()

        for i, file_path in enumerate(self.selected_files, start=1):
            filename = os.path.basename(file_path)

            self.disable_buttons()
            self.ui.btn_SelectFiles.setEnabled(False)
            self.ui.btn_Load.setEnabled(False)
            self.ui.radioButton.setEnabled(False)

            self.ui.comboBox_4.addItem(f"{filename}")
            self.ui.comboBox_4.setCurrentIndex(-1)
            
            self.ui.textEdit_6.setText(f"Analysing file {i} out of {len(self.selected_files)}")                
            
            if self.ui.checkBox_2.isChecked():
                if len(self.selected_files_gly) != len(self.selected_files):
                    QMessageBox.warning(self, "Invalid Data", f"The amount of reference files is not the same as sample files. Terminate analysis.", QMessageBox.Ok)
                    self.terminate()
                    self.ui.btn_SelectFiles.setEnabled(True)
                    return
                else: 
                    file_path_gly = self.selected_files_gly[i-1]
                    try:
                        self.process_file_data(file_path, file_path_gly, current_tab_index, i)
                    except:
                        self.analysis_error(file_path)
                        return
            else:
                try:
                    self.process_file_data(file_path, [], current_tab_index, i)
                except:
                    self.analysis_error(file_path)
                    return
                    
        self.finalize_analysis(legend, current_tab_index)

    def analysis_error(self, file_path):
        QMessageBox.warning(self, "Invalid Data", f"Error in {file_path}. Restarting analysis", QMessageBox.Ok)
        self.selected_files.remove(file_path)
        self.ui.btn_SelectFiles.setEnabled(True)
        self.analysis()

    def finalize_analysis(self, legend, current_tab_index):
        self.ui.textEdit_6.setText(f"Finished")
        legend.removeItem('Amplitude')
        legend.removeItem('Re')
        legend.removeItem('Im')
        legend.addItem(self.ui.FidWidget.plotItem.listDataItems()[0], name='Amp')
        legend.addItem(self.ui.FidWidget.plotItem.listDataItems()[1], name='Re')
        legend.addItem(self.ui.FidWidget.plotItem.listDataItems()[2], name='Im')
        self.enable_buttons()

        if current_tab_index == 1:
            self.update_dq_graphs()

    def process_file_data(self, file_path, file_path_gly, current_tab_index, i):
        global Frequency, Re_spectra, Im_spectra

        # Read data
        data = np.loadtxt(file_path)
        x, y, z = data[:, 0], data[:, 1], data[:, 2]

        # Read name of filename
        filename = os.path.basename(file_path)

        if self.ui.checkBox_2.isChecked():
            # Read data
            data = np.loadtxt(file_path)
            x, y, z = data[:, 0], data[:, 1], data[:, 2]
            # longcomponent
            Time_r, Re_r, Im_r = analysis_time_domain(file_path_gly)
            Time_s, Re_s, Im_s = analysis_time_domain(file_path)
            Time, Re, Im = long_component(Time_s, Time_r, Re_s, Re_r, Im_s, Im_r)
            
        else:
            Time, Re, Im = analysis_time_domain(file_path)
        
        Amp = calculate_amplitude(Re, Im)
        self.update_graphs(Time, Amp, Re, Im, self.ui.FidWidget)

        Time_fid, Fid =  final_analysis_time_domain(Time, Re, Im)

        Frequency = calculate_frequency_scale(Time_fid)
        # if self.ui.checkBox.isChecked():
        #     FFT = self.FFT_handmade(Fid, Time_fid, Frequency)  #(math procedure)
        # else:
        #     FFT = np.fft.fftshift(np.fft.fft(Fid))
        FFT = np.fft.fftshift(np.fft.fft(Fid))

        # 8. Simple baseline
        Amp_spectra, Re_spectra, Im_spectra = simple_baseline_correction(FFT)
        # 9. Apodization
        Real_apod = calculate_apodization(Re_spectra, Frequency)

        # Update FFT graphs
        self.update_graphs(Frequency, Amp_spectra, Re_spectra, Im_spectra, self.ui.FFTWidget)

        if self.ui.comboBox_4.currentIndex() == -1:
            M2, T2 = calculate_M2(Real_apod, Frequency)

            if current_tab_index == 0:
                match = re.search(r'.*_(-?\s*\d+\.?\d*).*.dat', filename)
                temperature = self.extract_info(match)

                SFC = calculate_SFC(Amp)
                self.ui.table_SE.setRowCount(len(self.selected_files))
                self.fill_table(self.ui.table_SE, temperature, SFC, M2, T2, i)

                if self.ui.radioButton.isChecked():
                    self.save_figures(file_path, temperature)

            elif current_tab_index == 1:
                match = re.search(r'_(\d+\.\d+)_', filename)
                dq_time = self.extract_info(match)

                Amplitude = calculate_amplitude(y, z)
                DQ = calculate_DQ_intensity(x, Amplitude)     
                self.ui.table_DQ.setRowCount(len(self.selected_files))       
                self.fill_table(self.ui.table_DQ, dq_time, DQ, M2, T2, i)

                if self.ui.radioButton.isChecked():
                    self.save_figures(file_path, dq_time)
        else:
            pass
     
    def after_phasing(self):
        global Frequency, Re_spectra, Im_spectra

        i = self.ui.comboBox_4.currentIndex()
        current_tab_index =  self.ui.tabWidget.currentIndex()

        Real_apod   = calculate_apodization(Re_spectra, Frequency) #(math procedure)
        Amp_spectra = calculate_amplitude(Re_spectra, Im_spectra)

        # Update FFT graph
        self.update_graphs(Frequency, Amp_spectra, Re_spectra, Im_spectra, self.ui.FFTWidget)

        M2, T2 = calculate_M2(Real_apod, Frequency)
        M2_r = round(M2, 6)
        T2_r = round(T2, 6)

        if current_tab_index == 0:
            table = self.ui.table_SE

            table.setItem(i, 2, QTableWidgetItem(str(M2_r)))
            table.setItem(i, 3, QTableWidgetItem(str(T2_r)))

            self.update_yaxis()

        elif current_tab_index == 1:
            table = self.ui.table_DQ

            table.setItem(i, 2, QTableWidgetItem(str(M2_r)))
            table.setItem(i, 3, QTableWidgetItem(str(T2_r)))

            #self.update_dq_graphs()

        # update table

    def extract_info(self, pattern):
        if pattern:
            info = pattern.group(1)
        else:
            info = '0'
        return info

    # Working with graphs
    def update_graphs(self, x, y1, y2, y3, graph):
        graph.clear()
        graph.plot(x, y1, pen='k')
        graph.plot(x, y2, pen='r')
        graph.plot(x, y3, pen='b')

    # Working with SE graphs
    def update_yaxis(self):

        axis_x = self.ui.comboBox_3.currentText()

        if axis_x == "T, C":
            self.ui.SEWidget.getAxis('bottom').setLabel("Temperature, °C")
            self.ui.table_SE.setHorizontalHeaderLabels(["Temp.", "SFC", "M₂", "T₂*"])
        elif axis_x == "XS, %":
            self.ui.SEWidget.getAxis('bottom').setLabel("XS, %")
            self.ui.table_SE.setHorizontalHeaderLabels(["XS", "SFC", "M₂", "T₂*"])
        else:
            return
        
        x = self.read_column_values(self.ui.table_SE, 0)

        text = self.ui.comboBox.currentText()

        if text == "SFC":  
            y =  self.read_column_values(self.ui.table_SE, 1)
            self.ui.SEWidget.getAxis('left').setLabel("SFC")
            
        elif text == "M₂":  
            y =  self.read_column_values(self.ui.table_SE, 2)
            self.ui.SEWidget.getAxis('left').setLabel("M₂")
            
        elif text == "T₂*":  
            y =  self.read_column_values(self.ui.table_SE, 3)
            self.ui.SEWidget.getAxis('left').setLabel("T₂*")

        else:  # Set y
            y =  [0]
            x = [0]
            self.ui.SEWidget.getAxis('left').setLabel("Not Set")
            return            
            
        self.ui.SEWidget.clear()
        self.ui.SEWidget.plot(x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=5)

    # Working with DQ graphs
    def update_dq_graphs(self):
        if len(self.selected_files) > 1:
            self.linearization()
            self.plot_fit()
        else:
            return

    def dq_t2_graph(self):
        x = self.read_column_values(self.ui.table_DQ, 0)
        y= self.read_column_values(self.ui.table_DQ, 3)
        self.ui.DQ_Widget_1.clear()
        self.ui.DQ_Widget_1.plot(x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=5)

    def linearization(self):
        # Calculate line
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
        
        Integral = np.trapz(dq)
        DQ_norm = dq/Integral

        self.ui.table_DQ.setColumnCount(5)
        self.ui.table_DQ.setColumnWidth(4,70)
        self.ui.table_DQ.setHorizontalHeaderItem(4, QTableWidgetItem("T₂* lin"))

        self.ui.table_DQ.setColumnCount(6)
        self.ui.table_DQ.setColumnWidth(5,70)
        self.ui.table_DQ.setHorizontalHeaderItem(5, QTableWidgetItem("DQ Norm"))

        for row in range(self.ui.table_DQ.rowCount()):
            T2_lin = coeff[0] * dq_time[row] + coeff[1]
            T2_lin = round(T2_lin, 4)
            item = QTableWidgetItem(str(T2_lin))
            self.ui.table_DQ.setItem(row, 4, item)
            item2 = QTableWidgetItem(str(round(DQ_norm[row], 4)))
            self.ui.table_DQ.setItem(row, 5, item2)


        # Draw line
        self.graph_line = self.ui.DQ_Widget_1
        
        if coeff is not None:
            # Generate x values for the line
            x_line = np.arange(0, 105.1, 0.1)
            
            # Calculate y values using the coefficients
            y_line = np.polyval(coeff, x_line)
            
            # Plot the graph and the line
            self.dq_t2_graph()
            self.graph_line.plot(x_line, y_line, pen='r')  
            self.t2_dq_graph()

    def t2_dq_graph(self):
        x = self.read_column_values(self.ui.table_DQ, 4)
        y= self.read_column_values(self.ui.table_DQ, 5)

        graph_DQ_distr = self.ui.DQ_Widget_2
        button = self.ui.radioButton_Log
        if button.isChecked():
            new_x = np.log10(x)
            graph_DQ_distr.getAxis('bottom').setLabel("log(T₂*)")
        else:
            new_x = x
            graph_DQ_distr.getAxis('bottom').setLabel("T₂*")

        graph_DQ_distr.clear()
        graph_DQ_distr.plot(new_x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=5)

    def plot_fit(self):        
        _x = self.read_column_values(self.ui.table_DQ, 4)
        y = self.read_column_values(self.ui.table_DQ, 5)

        button = self.ui.radioButton_Log
        if button.isChecked():
            x = np.log10(_x)
            p = [1, 1, 1, 0]
            b=([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1, np.inf])
            self.ui.DQ_Widget_2.getAxis('bottom').setLabel("log(T₂*)")
        else:
            x = _x
            p = [1, 5, 5, 0]
            b=([0, 0, 0, 0, 0], [ np.inf, np.inf, np.inf, 1, np.inf])
            self.ui.DQ_Widget_2.getAxis('bottom').setLabel("T₂*")

        text = self.ui.comboBox_2.currentText()
        
        x_fit = np.arange(0, np.max(x) + 0.001, 0.01)

        b1=([0, 0, 0, -10], [np.inf, np.inf, np.inf, np.inf])

        if text == 'Gauss':
            params, _ = curve_fit(self.gaussian, x, y, p0=p, bounds=b1)
            y_fit = self.gaussian(x_fit, *params)
            y_r2 = self.gaussian(x, *params)
            cen = params[1]
            fwhm = params[2]
            w = 0
            button.setEnabled(True)
        elif text == 'Lorenz':
            params, _ = curve_fit(self.lorenz, x, y, p0=p, bounds=b1)
            y_fit = self.lorenz(x_fit, *params)
            y_r2 = self.lorenz(x, *params)
            cen = params[1]
            fwhm = params[2]
            w = 1
            button.setEnabled(True)
        elif text == 'Pseudo Voigt':
            params, _ = curve_fit(self.voigt, x, y,  bounds = b)
            y_fit = self.voigt(x_fit, *params)
            y_r2 = self.voigt(x, *params)
            cen = params[1]
            fwhm = params[2]
            w = params[3]
            button.setEnabled(True)
        else:
            button.setEnabled(False)
            return


        #R2
        R = self.calculate_r_squared(y, y_r2)
        # Plot the line
        self.t2_dq_graph()
        self.ui.DQ_Widget_2.plot(x_fit, y_fit, pen='r')

        # Display R2, Xo and FWHM
        self.ui.textEdit_4.setText(f"R\u00B2: {round(R, 4)} \nX\u2080: {round(cen, 4)} \nFWHM: {round(fwhm, 4)} \nFraction (Lorenz): {round(w,2)}")

    def gaussian(self, x, amp, cen, wid, y0):
        # NoClass
        return amp * np.exp(-(x - cen)**2 / (2 * wid**2)) + y0

    def lorenz(self, x, amp, cen, wid, y0):
        # NoClass
        return (amp * (wid**2)) / ((x - cen)**2 + (wid**2)) + y0
    
    def voigt(self, x, amp, cen, wid, frac, y0):
        # NoClass
        lorentzian = amp *(2 * wid) / (np.pi * (4 * (x - cen)**2 + wid**2))
        gaussian = amp *(np.exp((-4 * np.log(2) * (x - cen)**2) / wid**2)) / (wid * np.sqrt(np.pi / (4 * np.log(2))))
        return  (frac * lorentzian + (1 - frac) * gaussian) + y0

    def calculate_r_squared(self, y_true, y_pred):
        # NoClass

        y_true = np.array(y_true)
        y_pred = np.array(y_pred)

        # Calculate the mean of the observed values
        y_mean = np.mean(y_true)

        # Calculate the total sum of squares (TSS)
        ss_tot = np.sum((y_true - y_mean) ** 2)

        # Calculate the residual sum of squares (RSS)
        ss_res = np.sum((y_true - y_pred) ** 2)

        # Calculate R-squared
        r_squared = 1 - (ss_res / ss_tot)

        return r_squared
    
    # T1 section
    def update_T1_table1(self):
        table = self.ui.table_T1
        table.setRowCount(len(self.selected_folders))
        self.ui.btn_Plot1.setEnabled(True)
        self.ui.btn_Plot2.setEnabled(True)

        

        for row, parent_folder in enumerate(self.selected_folders, start=0):
            filename = os.path.basename(parent_folder)
            files_in_folder = os.listdir(parent_folder)

            t1_files = [file for file in files_in_folder if file.startswith("T1")]

            if len(t1_files) != 1:
                QMessageBox.warning(self, "Invalid Data", f"The folder {filename} does not have one T1 file and will be removed from the table and file list.", QMessageBox.Ok)
                table.removeRow(row)
                del self.selected_folders[row]
            
        for row, parent_folder in enumerate(self.selected_folders, start=0):
            foldername = os.path.dirname(parent_folder)
            samplename = os.path.split(foldername)[1]
            filename = os.path.basename(parent_folder)
            files_in_folder = os.listdir(parent_folder)


            t1_files = [file for file in files_in_folder if file.startswith("T1")]


            try:
                full_path_to_data = os.path.join(parent_folder, t1_files[0])
                with open(full_path_to_data, "r") as file:
                    lines = [line.replace('\t\t\t', '').rstrip('\n') for line in file if not (line.rstrip('\n').endswith('\t\t\t\t'))]
                
                Time = []
                Signal = []

                for line in lines[1:]:  # Skip the first line !!!
                    parts = line.split('\t')
                    if len(parts) != 2:
                        raise ValueError("Data has more than 2 columns")
                    time_value = float(parts[0])
                    signal_value = float(parts[1])
                    Time.append(time_value)
                    Signal.append(signal_value)

                if filename not in self.t1_dictionary:
                    self.t1_dictionary[filename] = {"Time": [], "Signal": []}
                    

                self.t1_dictionary[filename]["Time"].extend(Time)
                self.t1_dictionary[filename]["Signal"].extend(Signal)

            except ValueError as e:
                QMessageBox.warning(self, "Invalid Data", f"I couldn't read {filename} due to: {str(e)}, removing file from the table and file list.", QMessageBox.Ok)
                table.removeRow(row)
                del self.selected_folders[row]
            except Exception as e:
                QMessageBox.warning(self, "Invalid Data", f"I couldn't read {filename}, removing file from the table and file list.", QMessageBox.Ok)
                table.removeRow(row)
                del self.selected_folders[row]

                
            
            Sample = QTableWidgetItem(samplename)
            Temperature = QTableWidgetItem(filename)
            
            table.setItem(row, 0, Sample)
            table.setItem(row, 1, Temperature)

            self.ui.comboBox_6.addItem(f"{filename}")
            self.ui.comboBox_6.setCurrentIndex(-1)


        # print("Time values for temperature", filename, ":", t1_dictionary[filename]["Time"])
        # print("Signal values for temperature", filename, ":", t1_dictionary[filename]["Signal"])
    
    def update_T1_table(self):

        for parent_folder in self.selected_folders:

            #foldername = os.path.dirname(parent_folder)
            #samplename = os.path.split(foldername)[1]
            filename = os.path.basename(parent_folder)
            files_in_folder = os.listdir(parent_folder)

            t1_files = [file for file in files_in_folder if file.startswith("T1")]

            if len(t1_files) != 1:
                QMessageBox.warning(self, "Invalid Data", f"The folder {filename} does not have one T1 file and will be removed from the table and file list.", QMessageBox.Ok)
                self.selected_folders = [folder for folder in self.selected_folders if folder != parent_folder]

            try:
                full_path_to_data = os.path.join(parent_folder, t1_files[0])
                with open(full_path_to_data, "r") as file:
                    lines = [line.replace('\t\t\t', '').rstrip('\n') for line in file if not (line.rstrip('\n').endswith('\t\t\t\t'))]
                
                Time = []
                Signal = []

                for line in lines[1:]:  # Skip the first line !!!
                    parts = line.split('\t')
                    if len(parts) != 2:
                        raise ValueError("Data has more than 2 columns")
                    time_value = float(parts[0])
                    signal_value = float(parts[1])
                    Time.append(time_value)
                    Signal.append(signal_value)

                if filename not in self.t1_dictionary:
                    self.t1_dictionary[filename] = {"Time": [], "Signal": []}
                    

                self.t1_dictionary[filename]["Time"].extend(Time)
                self.t1_dictionary[filename]["Signal"].extend(Signal)

            except ValueError as e:
                QMessageBox.warning(self, "Invalid Data", f"I couldn't read {filename} due to: {str(e)}, removing file from the table and file list.", QMessageBox.Ok)
                self.selected_folders = [folder for folder in self.selected_folders if folder != parent_folder]
            except Exception as e:
                QMessageBox.warning(self, "Invalid Data", f"I couldn't read {filename}, removing file from the table and file list.", QMessageBox.Ok)
                self.selected_folders = [folder for folder in self.selected_folders if folder != parent_folder]
            
        for row, parent_folder in enumerate(self.selected_folders, start=0):
            filename = os.path.basename(parent_folder)
            foldername = os.path.dirname(parent_folder)
            samplename = os.path.split(foldername)[1]


            table = self.ui.table_T1
            table.setRowCount(len(self.selected_folders))
            self.ui.btn_Plot1.setEnabled(True)
            self.ui.btn_Plot2.setEnabled(True)                
            
            Sample = QTableWidgetItem(samplename)
            Temperature = QTableWidgetItem(filename)
            
            table.setItem(row, 0, Sample)
            table.setItem(row, 1, Temperature)

            self.ui.comboBox_6.addItem(f"{filename}")
            self.ui.comboBox_6.setCurrentIndex(-1)

    def calculate_T1(self):
        selected_file_idx = self.ui.comboBox_6.currentIndex()

        if selected_file_idx == -1:
            return
        
        value_from_row = self.ui.table_T1.item(selected_file_idx, 1).text()
        Time = np.array(self.t1_dictionary[value_from_row]['Time'])/1000
        Signal = np.array(self.t1_dictionary[value_from_row]['Signal'])

        Time_fit = np.arange(min(Time), max(Time) + 1, 1)
        if self.ui.radioButton_2.isChecked():
            p = [-10, 200, 15]
            b=([-np.inf, 0, -np.inf], [np.inf, 50000, np.inf])
            popt, pcov = curve_fit(decaying_exponential, Time, Signal, p0 = p,bounds = b, maxfev=100000)
            fitted_curve = decaying_exponential(Time_fit, *popt)
            tau = round(popt[1],1)
            tau_str = str(tau)
            self.ui.textEdit_T1.setText(f"T1: {tau}")
            item = QTableWidgetItem(tau_str)
            self.ui.table_T1.setItem(selected_file_idx,2,item)
            item2 = QTableWidgetItem('0')
            self.ui.table_T1.setItem(selected_file_idx,3,item2)

        else:
            try:
                p = [-10, 200, -10, 200, 15]
                b=([-np.inf, 0, -np.inf, 0, -np.inf], [np.inf, 50000, np.inf, 50000, np.inf])
                popt, pcov = curve_fit(decaying_2exponential, Time, Signal, p0 = p, bounds = b, maxfev=100000)
                fitted_curve = decaying_2exponential(Time_fit, *popt)
                tau = round(popt[1],1)
                tau2 = round(popt[3],1)
                self.ui.textEdit_T1.setText(f"T1: {tau} \n T1: {tau2}")
                tau_str = str(tau)
                tau_str2 = str(tau2)
                item = QTableWidgetItem(tau_str)
                item2 = QTableWidgetItem(tau_str2)
                self.ui.table_T1.setItem(selected_file_idx,2,item)
                self.ui.table_T1.setItem(selected_file_idx,3,item2)
            except:
                QMessageBox.warning(self, "No covariance", f"I am sorry, I couldn't fit with two exponents. Fitting with one.", QMessageBox.Ok)
                p = [-10, 200, 15]
                b=([-np.inf, 0, -np.inf], [np.inf, 50000, np.inf])
                popt, pcov = curve_fit(decaying_exponential, Time, Signal, p0 = p,bounds=b, maxfev=100000)
                fitted_curve = decaying_exponential(Time_fit, *popt)
                tau = round(popt[1],1)
                tau_str = str(tau)
                self.ui.textEdit_T1.setText(f"T1: {tau}")
                item = QTableWidgetItem(tau_str)
                self.ui.table_T1.setItem(selected_file_idx,2,item)
                item2 = QTableWidgetItem('0')
                self.ui.table_T1.setItem(selected_file_idx,3,item2)

        
        figure = self.ui.T1_Widget_1
        figure.clear()
        figure.plot(Time, Signal, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=5)
        figure.plot(Time_fit, fitted_curve, pen='r')

    def plot_t1_temperature(self):
        table = self.ui.table_T1
        graph = self.ui.T1_Widget_2
        graph.clear()

        if table.rowCount() < 1:
            return

        Temperature = []
        T1 =[]
        for row in range(table.rowCount()):
            if table.item(row,1) and table.item(row,2) is not None:
                C = float(table.item(row, 1).text())
                if self.State_exp1 == True:
                    T = float(table.item(row, 2).text())
                else:
                    T = float(table.item(row, 3).text())
                Temperature.append(C)
                T1.append(T)
            else:
                print('oops')

        graph.plot(Temperature,T1, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=5)


    # Working with tables
    def update_DQ_comparison(self):
        table = self.ui.table_DQ_2
        table.setRowCount(len(self.selected_DQfiles))
        self.ui.btn_Launch.setEnabled(True)

        for row, parent_folder in enumerate(self.selected_DQfiles, start=0):
            foldername = os.path.dirname(parent_folder)
            filename = os.path.basename(parent_folder)

            try: 
                data = np.loadtxt(parent_folder, delimiter=',')
                if data.shape[1] < 3:
                    QMessageBox.warning(self, "Invalid Data", f"The file {foldername} does not have at least 3 columns and will be removed from the table and file list.", QMessageBox.Ok)
                    table.removeRow(row)
                    del self.selected_DQfiles[row]

            except: 
                QMessageBox.warning(self, "Invalid Data", f"The file {foldername} does not have , as delimiter and will be removed from the table and file list.", QMessageBox.Ok)
                table.removeRow(row)
                del self.selected_DQfiles[row]

                
            item = QTableWidgetItem(filename)
            default_name = QTableWidgetItem(str(row+1))
            table.setItem(row, 0, item)
            table.setItem(row, 1, default_name)

    def launch(self):
        self.ui.radioButton_Log_2.setEnabled(True)
        self.ui.comboBox_5.setEnabled(True)
        self.ui.dq_min_2.setEnabled(True)
        self.ui.dq_max_2.setEnabled(True)

        self.dq_t2 = {}
        for row, parent_folder in enumerate(self.selected_DQfiles, start=0):
            # Read data from file
            data = np.loadtxt(parent_folder, delimiter=',')

            # Read the DQ filtering time, DQ amlitude and corresponding T2*
            dq_t2 = data[:, [0, 1, 3]]
            self.dq_t2[row] = dq_t2
        self.update_DQ_comparison_plot()

    def update_DQ_comparison_plot(self):
        cmap = pg.ColorMap([0, len(self.dq_t2)], [pg.mkColor('b'), pg.mkColor('r')])  # Blue to red

        #legend = self.ui.DQ_Widget_3.addLegend()
        legend1 = self.ui.DQ_Widget_4.addLegend()  # Get the legend object
        self.ui.DQ_Widget_3.clear()
        self.ui.DQ_Widget_4.clear()
        self.ui.DQ_Widget_5.clear()

        if legend1 is not None:
            #legend.clear()
            legend1.clear()
            #self.ui.DQ_Widget_3.addLegend()
            self.ui.DQ_Widget_4.addLegend()
            #legend.setPen((0, 0, 0))  
            legend1.setPen((0, 0, 0))

        # This is a mess :(

        center = []
        comparison_par = []
        for row, (key, data) in zip(range(self.ui.table_DQ_2.rowCount()), self.dq_t2.items()):
            file_name_item = self.ui.table_DQ_2.item(row, 1)
            if file_name_item is not None:
                file_name = file_name_item.text()
                if file_name != 'hide':
                    try:
                        _comp_par = float(file_name)
                    except:
                        _comp_par = 0
                        self.ui.table_DQ_2.setItem(row, 1, QTableWidgetItem('0'))
                    comparison_par.append(_comp_par)
                else:
                    continue

            # Initial arrays


            dq_time = data[:,0] #DQ filtering time
            dq = data[:,1] #DQ amlitude
            t2 = data[:,2] #T2*

            # Linear
            # Read boxes for ALL the graph just opne for the sake of consistency
            time_min = self.ui.dq_min_2.value()
            time_max = self.ui.dq_max_2.value()

            t2_x = dq_time[(dq_time >= time_min) & (dq_time <= time_max)]
            if len(t2_x) < 3:
                time_min = 0
                time_max = 20
                self.ui.dq_min_2.setValue(time_min)
                self.ui.dq_max_2.setValue(time_max)
                t2_x = dq_time[(dq_time >= time_min) & (dq_time <= time_max)]


            t2_y = t2[(dq_time >= time_min) & (dq_time <= time_max)]
            coeff = np.polyfit(t2_x, t2_y, 1)
            t2_lin_ = coeff[0] * dq_time + coeff[1] # Will it work?


            # Dq time on T2
            Integral = np.trapz(dq)
            dq_norm = dq/Integral

            if self.ui.radioButton_Log_2.isChecked():
                t2_lin = np.log10(t2_lin_)
                self.ui.DQ_Widget_4.getAxis('bottom').setLabel("log(T₂*)")
                p = [1, 1, 1, 0]
                b=([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1, np.inf])
            else:
                t2_lin = t2_lin_
                self.ui.DQ_Widget_4.getAxis('bottom').setLabel("T₂*")
                p = [1, 5, 5, 0]
                b=([0, 0, 0, 0, 0], [ np.inf, np.inf, np.inf, 1, np.inf])

            text = self.ui.comboBox_5.currentText()
            dq_fit = np.arange(0, np.max(t2_lin) + 0.001, 0.01)
            b1=([0, 0, 0, -10], [np.inf, np.inf, np.inf, np.inf])

            if text == 'Gauss':
                params, _ = curve_fit(self.gaussian, t2_lin, dq_norm, p0=p, bounds=b1)
                y_fit = self.gaussian(dq_fit, *params)
                y_r2 = self.gaussian(t2_lin, *params)
                cen = params[1]
                center.append(cen)
                fwhm = params[2]
                w = 0
            elif text == 'Lorenz':
                params, _ = curve_fit(self.lorenz, t2_lin, dq_norm, p0=p, bounds=b1)
                y_fit = self.lorenz(dq_fit, *params)
                y_r2 = self.lorenz(t2_lin, *params)
                cen = params[1]
                center.append(cen)
                fwhm = params[2]
                w = 1
            elif text == 'Pseudo Voigt':
                params, _ = curve_fit(self.voigt, t2_lin, dq_norm,  bounds = b)
                y_fit = self.voigt(dq_fit, *params)
                y_r2 = self.voigt(t2_lin, *params)
                cen = params[1]
                center.append(cen)
                fwhm = params[2]
                w = params[3]

            # Draw a graph
            color = tuple(cmap.map(key))
            self.ui.DQ_Widget_3.plot(dq_time, t2, pen=None, symbolPen=None, symbol='o', symbolBrush=color, symbolSize=5, name=file_name)
            self.ui.DQ_Widget_4.plot(t2_lin, dq_norm, pen=None, symbolPen=None, symbol='o', symbolBrush=color, symbolSize=5, name=file_name)
            self.ui.DQ_Widget_4.plot(dq_fit, y_fit, pen=color)

            if coeff is not None:
                x_line = np.arange(0, 105.1, 0.1)
                y_line = np.polyval(coeff, x_line)
                
                self.ui.DQ_Widget_3.plot(x_line, y_line, pen=color) 

        self.ui.DQ_Widget_5.plot(comparison_par, center, pen=None, symbolPen=None, symbol='o', symbolBrush='r')

    def nan_value(self, table, row, column_index):
        table.setItem(row, column_index, QTableWidgetItem('NaN'))

    def fill_table(self, table, c1, c2, c3, c4, i):
        j = i-1

        c2 = round(c2, 6)
        c3 = round(c3, 6)
        c4 = round(c4, 6)

        table.setItem(j, 0, QTableWidgetItem(c1))
        table.setItem(j, 1, QTableWidgetItem(str(c2)))
        table.setItem(j, 2, QTableWidgetItem(str(c3)))
        table.setItem(j, 3, QTableWidgetItem(str(c4)))
    
    def read_column_values(self, table, column_index):
        column_values = []
        for row in range(table.rowCount()):
            item = table.item(row, column_index)
            if item is not None and item.text() != '':
                column_values.append(float(item.text()))  # Assuming the values are numeric!!!!!
            else:
                self.nan_value(table, row, column_index)
        return column_values

    # Save and load data
    def save_data(self):
        parent_folder = os.path.dirname(self.selected_files[0])

        current_tab_index = self.ui.tabWidget.currentIndex()
     
        DQ_points = self.ui.DQ_Widget_1
        DQ_distribution = self.ui.DQ_Widget_2

        SE_table = self.ui.table_SE
        DQ_table = self.ui.table_DQ

        os.makedirs(parent_folder + '/Result/', exist_ok=True)
        basefolder = os.path.dirname(parent_folder)
        os.makedirs(basefolder + '/DQ_table/', exist_ok=True)
        if current_tab_index == 0:
            table_file_path = os.path.join(parent_folder, 'Result', f"SE_table.csv")
        elif current_tab_index == 1:
            graph_file_path = os.path.join(parent_folder, 'Result', f"DQ_points.png")
            graph_file_path_2 = os.path.join(parent_folder, 'Result', f"DQ_distribution.png")
            table_file_path = os.path.join(basefolder, 'DQ_table', f"DQ_table_{os.path.basename(parent_folder)}.csv")

        pg.QtGui.QGuiApplication.processEvents()  # Make sure all events are processed before exporting

            # Check if the file already exists, if so, append an index to the filename
        if os.path.exists(table_file_path):
            index = 1
            while os.path.exists(f"{table_file_path[:-4]} ({index}).csv"):
                index += 1
            table_file_path = f"{table_file_path[:-4]} ({index}).csv"

        if current_tab_index == 0:
            #Table
            self.save_table_to_csv(table_file_path, SE_table)

        elif current_tab_index == 1:
            #Image
            exporter_dq1 = pg.exporters.ImageExporter(DQ_points.plotItem)
            exporter_dq1.parameters()['width'] = 1000
            exporter_dq1.export(graph_file_path)

            exporter_dq2 = pg.exporters.ImageExporter(DQ_distribution.plotItem)
            exporter_dq2.parameters()['width'] = 1000
            exporter_dq2.export(graph_file_path_2)

            #Table
            self.save_table_to_csv(table_file_path, DQ_table)

    def save_table_to_csv(self, path, table):
        with open(path, 'w') as f:
            # Write data row by row
            for row in range(table.rowCount()):
                row_values = []
                for col in range(table.columnCount()):
                    item = table.item(row, col)
                    if item is not None:
                        row_values.append(item.text())
                    else:
                        row_values.append("")  # Handle empty cells
                f.write(','.join(row_values) + '\n')

    def load_data(self):
        dlg = OpenFilesDialog(self)
        self.ui.btn_Save.setEnabled(False)
        if dlg.exec():
            tableName = dlg.selectedFiles()
            self.selected_table = tableName
            self.load_table_from_csv(tableName)

    def load_table_from_csv(self, tableName):
        current_tab_index = self.ui.tabWidget.currentIndex()
        if current_tab_index == 0:
            table = self.ui.table_SE
        elif current_tab_index == 1:
            table = self.ui.table_DQ

        file_path = tableName[0]

        while self.ui.comboBox_4.count()>0:
            self.ui.comboBox_4.removeItem(0)
            self.ui.comboBox_4.setCurrentIndex(-1)

        with open(file_path, 'r') as f:
            lines = f.readlines()
            table.setRowCount(len(lines))
            for row, line in enumerate(lines):
                values = line.strip().split(',')
                self.ui.comboBox_4.addItem(f"File #{row+1}")
                for col, value in enumerate(values):
                    item = QTableWidgetItem(value)
                    table.setItem(row, col, item)

        if current_tab_index == 0:
            self.update_yaxis()
        elif current_tab_index == 1:
            self.selected_files = [
                "AmIanIdiot.txt",
                "Yes.txt"
            ]
            self.update_dq_graphs()  

        self.enable_buttons()
        self.ui.comboBox_4.setEnabled(False)
        #self.ui.checkBox.setEnabled(False)
        self.ui.btn_Phasing.setEnabled(False)

    def save_figures(self, file_path, variable):
        # Set names
        parent_folder = os.path.dirname(file_path)  

        graph_fft = self.ui.FFTWidget
        graph_fid = self.ui.FidWidget

        fft_file_path = (parent_folder + '/Result/' + f"FFT_{variable}.png")
        fid_file_path = (parent_folder + '/Result/' + f"FID_{variable}.png")

        pg.QtGui.QGuiApplication.processEvents()

        exporter_fft = pg.exporters.ImageExporter(graph_fft.plotItem)
        exporter_fft.parameters()['width'] = 1000
        exporter_fft.export(fft_file_path)

        exporter_fid = pg.exporters.ImageExporter(graph_fid.plotItem)
        exporter_fid.parameters()['width'] = 1000
        exporter_fid.export(fid_file_path)

    def load_data_and_check_validity(self, file_path, current_tab_index):
        try: 
            data = np.loadtxt(file_path)
            filename = os.path.basename(file_path)

            # Check that there are 3 columns
            if data.shape[1] != 3:
                dialog = AlertDialog()
                if dialog.exec() == QDialog.Rejected:
                    self.ui.btn_SelectFiles.setEnabled(True)
                    self.ui.radioButton.setEnabled(True)
                    self.selected_files.clear()
                    return False
                
            elif not ((filename.startswith("SE") or filename.startswith("XS")) and current_tab_index == 0) and \
                not (filename.startswith("DQ") and current_tab_index == 1) and not current_tab_index == 2:
                dialog = NotificationDialog()
                if dialog.exec() == QDialog.Rejected:
                    self.ui.btn_SelectFiles.setEnabled(True)
                    self.ui.radioButton.setEnabled(True)
                    return False
            return True
        except:
            QMessageBox.warning(self, "Invalid Data", f"I can't read the {file_path} file, deleting it.", QMessageBox.Ok)
            self.selected_files.remove(file_path)

            if self.selected_files == []:
                self.disable_buttons()
                return False
            else:
                return False

    # Math procedures
    def FFT_handmade(self, Fid, Time, Freq):
        N = len(Freq)
        Fur = np.zeros(N, dtype=complex)


        cos_values = np.cos(2 * np.pi * Time[:, None] * Freq)
        sin_values = np.sin(2 * np.pi * Time[:, None] * Freq)

        for i in range(N):
            c = i + 1
            progress = round((c / N) * 100, 3)
            self.ui.progressBar.setValue(progress)

            Fur[i] = np.sum(Fid * (cos_values[:, i] - 1j * sin_values[:, i]))
        QCoreApplication.processEvents()
        return Fur

class OpenFilesDialog(QFileDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFileMode(QFileDialog.ExistingFiles)  # Allow selecting multiple files
        self.setNameFilter(str("Data (*.dat *.txt *.csv)"))
        self.setDirectory(str("C:/Mega/NMR/003_Temperature"))
        
        self.selected_files = []  # Variable to store selected file paths

        

    def on_file_selected(self):
        options = QFileDialog.Options()
        files, _ = QFileDialog.getOpenFileNames(self, "Load Files", "", "Data Files (*.dat *.txt *.csv)", options=options)
        if files:
            self.selected_files.extend(files)

class NotificationDialog(QDialog, Ui_Note):
    stateChanged = Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.buttonBox.accepted.connect(self.accepted_handler)
        self.buttonBox.rejected.connect(self.rejected_handler)
        self.flag_ok = False

    def accepted_handler(self):
        self.flag_ok = True
        self.stateChanged.emit(self.flag_ok)
        self.accept()

    def rejected_handler(self):
        self.flag_ok = False
        self.stateChanged.emit(self.flag_ok)
        self.reject()

class AlertDialog(QDialog, Ui_Error):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.pushButton.clicked.connect(self.close_dialog) 

    def close_dialog(self):
        self.reject() 
        
class PhasingManual(QDialog):
    # TODO somehow create an araay of ORIGINAL Re_spectra and be able to restore it...
    closed = Signal()
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_PhasingManual()
        self.ui.setupUi(self)

        graph_phasing = self.ui.PhasingGraph
        graph_phasing.getAxis('bottom').setLabel("Frequency, MHz")
        graph_phasing.getAxis('left').setLabel("Amplitude, a.u.")
        graph_phasing.setTitle("Phasing")
        self.ui.PhasingGraph.addLegend()

        self.ui.pushButton_2.clicked.connect(self.zero)
        self.ui.verticalSlider_a.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_b.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_c.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_d.valueChanged.connect(self.value_changed)
        self.ui.dial.valueChanged.connect(self.smoothing_changed)
        self.ui.pushButton_3.clicked.connect(self.manual_read)
        self.ui.pushButton.clicked.connect(self.save_data)

        self.Real_freq_phased = None

    def read_data(self):
        global Frequency, Re_spectra, Im_spectra

        self.zero()

    def save_data(self):
        global Re_spectra

        if self.Real_freq_phased is not None:
            Re_spectra = self.Real_freq_phased
        else: 
            Re_spectra = Re_spectra
        

        self.close()

    def closeEvent(self, event):
        self.closed.emit()
        super().closeEvent(event)

    def zero(self):
        self.a = 0
        self.b = 0
        self.c = 0
        self.d = 0
        self.Smooth = 0

        self.ui.verticalSlider_a.setValue(0)
        self.ui.verticalSlider_b.setValue(0)
        self.ui.verticalSlider_c.setValue(0)
        self.ui.verticalSlider_d.setValue(0)
        self.ui.Box_a.setValue(0)
        self.ui.Box_b.setValue(0)
        self.ui.Box_c.setValue(0)
        self.ui.Box_d.setValue(0)
        self.ui.dial.setValue(0)
        self.ui.Box_smooth.setValue(0)

        if self.Real_freq_phased is not None:
            self.process_data()
    
    def process_data(self):
        self.Real_freq_phased = self.calculate_phase()
        self.update_plot()
        self.update_text()
    
    def calculate_phase(self):
        global Frequency, Re_spectra, Im_spectra
        phi = self.a + self.b * Frequency + self.c * Frequency ** 2 + self.d * Frequency ** 3
        self.Real_freq_phased = Re_spectra * np.cos(np.deg2rad(phi)) - Im_spectra * np.sin(np.deg2rad(phi))

        if self.Smooth > 1:
            self.Smooth = int(self.Smooth)
            self.Real_freq_phased = self.smooth(self.Real_freq_phased)

        return self.Real_freq_phased

    def smooth(self,y,):
        new_array = savgol_filter(y, window_length=self.Smooth, polyorder=1)

        return new_array

    def update_plot(self):
        global Frequency, Re_spectra, Im_spectra

        self.ui.PhasingGraph.clear()
        self.ui.PhasingGraph.plot(Frequency, Re_spectra, pen='r', name = 'Original') #Real origin
        self.ui.PhasingGraph.plot(Frequency, self.Real_freq_phased, pen='b', name = 'Phased') #Real phased

    def update_text(self):
        Integral = np.trapz(self.Real_freq_phased)
        left_mean = np.mean(self.Real_freq_phased[:100])
        right_mean = np.mean(self.Real_freq_phased[-100:])
        delta = left_mean - right_mean

        self.ui.Integral.setText(f"Integral: {round(Integral,3)}")
        self.ui.Delta.setText(f"Delta: {round(delta,7)}")

    def value_changed(self):
        global Frequency, Re_spectra, Im_spectra

        self.a = self.ui.verticalSlider_a.value()
        self.b = self.ui.verticalSlider_b.value()
        self.c = self.ui.verticalSlider_c.value()
        self.d = self.ui.verticalSlider_d.value()

        self.ui.Box_a.setValue(self.a)
        self.ui.Box_b.setValue(self.b)
        self.ui.Box_c.setValue(self.c)
        self.ui.Box_d.setValue(self.d)


        if Frequency is not None and Re_spectra is not None and Im_spectra is not None:
            self.process_data()
        else:
            return

    def smoothing_changed(self):
        self.Smooth = self.ui.dial.value()
        self.ui.Box_smooth.setValue(self.Smooth)
        self.process_data()

    def manual_read(self):
        self.a = self.ui.Box_a.value()
        self.b = self.ui.Box_b.value()
        self.c = self.ui.Box_c.value()
        self.d = self.ui.Box_d.value()
        self.Smooth = self.ui.Box_smooth.value()

        self.ui.verticalSlider_a.valueChanged.disconnect(self.value_changed)
        self.ui.verticalSlider_b.valueChanged.disconnect(self.value_changed)
        self.ui.verticalSlider_c.valueChanged.disconnect(self.value_changed)
        self.ui.verticalSlider_d.valueChanged.disconnect(self.value_changed)
        self.ui.dial.valueChanged.disconnect(self.smoothing_changed)
        

        self.ui.dial.setValue(self.Smooth)
        self.ui.verticalSlider_a.setValue(self.a)
        self.ui.verticalSlider_b.setValue(self.b)
        self.ui.verticalSlider_c.setValue(self.c)
        self.ui.verticalSlider_d.setValue(self.d)

        self.process_data()

        self.ui.verticalSlider_a.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_b.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_c.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_d.valueChanged.connect(self.value_changed)
        self.ui.dial.valueChanged.connect(self.smoothing_changed)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    sys.exit(app.exec())
