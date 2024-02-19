# This Python file uses the following encoding: utf-8
import sys, os, re

from PySide6.QtWidgets import QApplication, QMainWindow, QFileDialog, QTableWidgetItem, QDialog
from PySide6.QtCore import QCoreApplication, Signal
import numpy as np
from scipy.optimize import curve_fit
import pyqtgraph as pg
from pyqtgraph.exporters import ImageExporter
from ui_Form import Ui_Analyzer
from ui_ChooseFiles import Ui_ChooseFiles
from ui_Notification import Ui_Dialog
from ui_Error import Ui_Error

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_Analyzer()
        self.ui.setupUi(self)

        self.ui.btn_SelectFiles.clicked.connect(self.open_select_dialog)

        self.ui.btn_Start.clicked.connect(self.analysis)

        self.ui.btn_Save.clicked.connect(self.save_data)
        self.ui.btn_Load.clicked.connect(self.load_data)
        #Graphs
        # Seetings for FFT
        graph_fft = self.ui.FFTWidget
        graph_fft.getAxis('left').setLabel("Amplitude, a.u")
        graph_fft.getAxis('bottom').setLabel("Frequency, MHz")
        graph_fft.setBackground('w')
        graph_fft.setTitle("FFT") 

        # Settings for graph with FID
        graph_fid = self.ui.FidWidget
        graph_fid.getAxis('left').setLabel("Amplitude")
        graph_fid.getAxis('bottom').setLabel("Time, ms")
        graph_fid.setBackground('w')
        graph_fid.setTitle("FID")   

        # Settings for graph SE temperature
        graph_SE = self.ui.SEWidget
        graph_SE.getAxis('bottom').setLabel("Temperature, °C")
        graph_SE.setBackground('w')
        graph_SE.setTitle("Temperature dependence")   

        # Settings for graph DQ points
        graph_DQ_points = self.ui.DQ_Widget_1
        graph_DQ_points.getAxis('bottom').setLabel("DQ Filtering Time")
        graph_DQ_points.getAxis('left').setLabel("T2*")
        graph_DQ_points.setBackground('w')
        # Initialize coefficients
        self.coeff = None
        self.line_item = None
          

        # Settings for graph DQ distribution
        graph_DQ_distr = self.ui.DQ_Widget_2
        graph_DQ_distr.getAxis('left').setLabel("Norm. DQ Intensity")
        graph_DQ_distr.setBackground('w')   

        # Tables
        table_SE = self.ui.table_SE
        table_SE.setColumnWidth(0,70)
        table_SE.setColumnWidth(1,70)
        table_SE.setColumnWidth(2,70)
        table_SE.setColumnWidth(3,70)

        table_DQ = self.ui.table_DQ
        table_DQ.setColumnWidth(0,70)
        table_DQ.setColumnWidth(1,70)
        table_DQ.setColumnWidth(2,70)
        table_DQ.setColumnWidth(3,70)

        # Buttons
        self.ui.btn_Start.setEnabled(False)
        self.ui.btn_Save.setEnabled(False)
        
        self.ui.radioButton_Log.clicked.connect(self.t2_dq_graph)


        # Draw graph Temperature vs combobox option
        self.ui.comboBox.currentIndexChanged.connect(self.updateGraph)
        self.ui.comboBox_2.currentIndexChanged.connect(self.plot_fit)

        # Change event
        self.ui.dq_min.valueChanged.connect(self.linearization)
        self.ui.dq_max.valueChanged.connect(self.linearization)


    def open_select_dialog(self):
        dlg = OpenFilesDialog(self)
        if dlg.exec():
            fileNames = dlg.selectedFiles()
            self.selected_files = fileNames
        self.ui.btn_Start.setEnabled(True)


    def analysis(self):
        current_tab_index = self.ui.tabWidget.currentIndex()
        
        graph_fft = self.ui.FFTWidget
        graph_fid = self.ui.FidWidget

        i = 1

        self.ui.table_SE.setRowCount(len(self.selected_files))
        self.ui.table_DQ.setRowCount(len(self.selected_files))

        # Add legend
        legend = pg.LegendItem(offset=(300, 10))  # Adjust the offset as needed
        legend.setParentItem(graph_fid.graphicsItem()) 
        
        for file_path in self.selected_files:
            parent_folder = os.path.dirname(file_path)
            
            # Disable all Buttons during analysis
            self.ui.btn_SelectFiles.setEnabled(False)
            self.ui.btn_Start.setEnabled(False)
            self.ui.radioButton.setEnabled(False)


            self.ui.textEdit_4.clear()
            
            # Read data
            data = np.loadtxt(file_path)  
            if data.shape[1] != 3:
                # Open AlertDialog
                dialog = AlertDialog()
                if dialog.exec() == QDialog.Rejected:
                    # Clear data files and return
                    self.ui.btn_SelectFiles.setEnabled(True)
                    self.ui.radioButton.setEnabled(True)

                    self.selected_files.clear()
                    return
                
            x, y, z = data[:, 0], data[:, 1], data[:, 2]

            filename = re.search(r'[^/\\]+$', file_path).group()  # Extract filename using regex

            if filename.startswith("SE") is False and current_tab_index == 0:
                dialog = NotificationDialog()
                if dialog.exec() == QDialog.Rejected:
                    self.ui.btn_SelectFiles.setEnabled(True)
                    self.ui.btn_Start.setEnabled(False)
                    self.ui.radioButton.setEnabled(True)
                    self.ui.btn_Save.setEnabled(False)

                    return
            elif filename.startswith("DQ") is False and current_tab_index == 1:
                dialog = NotificationDialog()
                if dialog.exec() == QDialog.Rejected:
                    self.ui.btn_SelectFiles.setEnabled(True)
                    self.ui.btn_Start.setEnabled(False)
                    self.ui.btn_Save.setEnabled(False)
                    self.ui.radioButton.setEnabled(True)
                    return


            
            self.ui.textEdit_4.append(f"File: {filename}") # Write the filename
            self.ui.textEdit_6.setText(f"Analysing file {i} out of {len(self.selected_files)}")

            # Crop, phase, normalize
            Time_p, Amp, Re, Im = self.pre_processing(x, y, z)
            
            #plot FID
            graph_fid.clear()
            graph_fid.plot(Time_p, Amp, pen='k')
            graph_fid.plot(Time_p, Re, pen='r')
            graph_fid.plot(Time_p, Im, pen='b')

            # Apodization in time domain:
            Re_ap, Im_ap = self.apodization(Time_p, Amp, Re, Im)

            # Add zeros
            numberp = 16384
            Time, Fid = self.add_zeros(Time_p, Re_ap, Im_ap, numberp)

            # Recalculate time domain to frequency domain:
            Frequency = self.calculate_frequency_scale(Time)

            # Perform FFT
            QCoreApplication.processEvents()
            FFT = self.FFT_handmade(Fid, Time, Frequency)

            # Baseline:
            if len(Frequency) != len(FFT):
            # Adjust the length of Time to match Fid
                Frequency = np.linspace(Frequency[0], Frequency[-1], len(FFT))
                print('olalala')
            Amp_spectra, Re_spectra, Im_spectra = self.simple_baseline_correction(FFT)

            # Perform the apodization of the spectra:
            Real_apod = self.calculate_apodization(Re_spectra, Frequency)

            #plot FFT
            graph_fft.clear()
            graph_fft.plot(Frequency, Amp_spectra, pen='k')
            graph_fft.plot(Frequency, Re_spectra, pen='r')
            graph_fft.plot(Frequency, Im_spectra, pen='b')

            #Fill the table
            M2 = self.calculate_M2(Real_apod, Frequency)
            T2 = self.calculate_T2(M2)

            if current_tab_index == 0 and filename.startswith("SE"):
                Temperature =  self.extract_temperature(filename)
                SFC = self.calculate_SFC(Amp)
                self.fill_table_SE(SFC, M2, T2, Temperature, i, self.ui.table_SE)
            elif current_tab_index == 1 and filename.startswith("DQ"):
                dq_filtering_time = self.extract_dq_filtering_time(filename)
                Amplitude = self.calculate_amplitude(y, z)
                DQ = self.calculate_DQ_intensity(x, Amplitude)
                self.fill_table_DQ(DQ, M2, T2, dq_filtering_time, i, self.ui.table_DQ)
            
            i=i+1
        

            # Save figures
            if self.ui.radioButton.isChecked():
                if not os.path.exists(os.path.join(parent_folder, 'Result')):
                    os.makedirs(os.path.join(parent_folder, 'Result'))


                fft_widget = self.ui.FFTWidget
                fid_widget = self.ui.FidWidget
                if current_tab_index == 0 and filename.startswith("SE"):
                    fft_file_path = os.path.join(parent_folder, 'Result', f"FFT_{Temperature}.png")
                    fid_file_path = os.path.join(parent_folder, 'Result', f"FID_{Temperature}.png")
                elif current_tab_index == 1 and filename.startswith("DQ"):
                    fft_file_path = os.path.join(parent_folder, 'Result', f"FFT_{dq_filtering_time}.png")
                    fid_file_path = os.path.join(parent_folder, 'Result', f"FID_{dq_filtering_time}.png")

                if os.path.exists(fft_file_path):
                    os.remove(fft_file_path)
                if os.path.exists(fid_file_path):
                    os.remove(fid_file_path)

                pg.QtGui.QGuiApplication.processEvents()  # Make sure all events are processed before exporting

                exporter_fft = pg.exporters.ImageExporter(fft_widget.plotItem)
                exporter_fft.parameters()['width'] = 1000
                exporter_fft.export(fft_file_path)

                exporter_fid = pg.exporters.ImageExporter(fid_widget.plotItem)
                exporter_fid.parameters()['width'] = 1000
                exporter_fid.export(fid_file_path)

            
        self.ui.textEdit_6.setText(f"Finished")

        legend.removeItem('Amplitude')
        legend.removeItem('Re')
        legend.removeItem('Im')
        legend.addItem(graph_fid.plotItem.listDataItems()[0], name='Amp')
        legend.addItem(graph_fid.plotItem.listDataItems()[1], name='Re')
        legend.addItem(graph_fid.plotItem.listDataItems()[2], name='Im')   

        # Enable all Buttons after analysis
        self.ui.btn_SelectFiles.setEnabled(True)
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Save.setEnabled(True)
        self.ui.radioButton.setEnabled(True)

        if current_tab_index == 1 and filename.startswith("DQ"):
            self.dq_t2_graph()
            self.t2_dq_graph()  
            self.linearization() 

    def linearization(self):
        time_min = self.ui.dq_min.value()
        time_max = self.ui.dq_max.value()
        dq_time = np.array(self.read_column_values(self.ui.table_DQ, 0))
        t2 = np.array(self.read_column_values(self.ui.table_DQ, 3))
        dq = np.array(self.read_column_values(self.ui.table_DQ, 1))

        x = dq_time[(dq_time >= time_min) & (dq_time <= time_max)]
        y = t2[(dq_time >= time_min) & (dq_time <= time_max)]

        if len(x) <= 1 or len(y) <= 1:
            print("Not enough data points")
            return

        self.coeff = np.polyfit(x, y, 1)
        

        Integral = np.trapz(dq)
        DQ_norm = dq/Integral

        self.ui.table_DQ.setColumnCount(5)
        self.ui.table_DQ.setColumnWidth(4,70)
        self.ui.table_DQ.setHorizontalHeaderItem(4, QTableWidgetItem("T2* lin"))

        self.ui.table_DQ.setColumnCount(6)
        self.ui.table_DQ.setColumnWidth(5,70)
        self.ui.table_DQ.setHorizontalHeaderItem(5, QTableWidgetItem("DQ Norm"))

        for row in range(self.ui.table_DQ.rowCount()):
            T2_lin = self.coeff[0] * dq_time[row] + self.coeff[1]
            T2_lin = round(T2_lin, 4)
            item = QTableWidgetItem(str(T2_lin))
            self.ui.table_DQ.setItem(row, 4, item)
            item2 = QTableWidgetItem(str(round(DQ_norm[row], 4)))
            self.ui.table_DQ.setItem(row, 5, item2)
    
        self.t2_dq_graph()
        self.plot_line()

    def plot_line(self):
        self.graph_line = self.ui.DQ_Widget_1
        
        if self.coeff is not None:
            # Generate x values for the line
            x_line = np.arange(0, 105.1, 0.1)
            
            # Calculate y values using the coefficients
            y_line = np.polyval(self.coeff, x_line)
            
            # Plot the line
            self.dq_t2_graph()
            self.graph_line.plot(x_line, y_line, pen='r')  

    def plot_fit(self, index):
        if index == 0:
            self.plot_gauss()
        elif index == 1:
            self.plot_lorenz()
        elif index == 2:
            self.plot_voigt()
        else:
            return

    def plot_gauss(self):
        _x = self.read_column_values(self.ui.table_DQ, 4)
        y = self.read_column_values(self.ui.table_DQ, 5)

        button = self.ui.radioButton_Log
        if button.isChecked():
            x = np.log10(_x)
        else:
            x = _x

        initial_guess = [10**(-4), 10, 1]  # Initial guess for the parameters [amplitude, center, width]
        params, covariance = curve_fit(self.gaussian, x, y, p0=initial_guess)

        amp_fit, cen_fit, wid_fit = params

        x_fit = np.arange(0, 105.1, 0.1)
        y_fit = self.gaussian(x_fit, *params)

        # Plot the line
        self.t2_dq_graph()
        self.ui.DQ_Widget_2.plot(x_fit, y_fit, pen='r')

    def gaussian(self, x, amp, cen, wid):
        return amp * np.exp(-(x - cen)**2 / (2 * wid**2))

    def plot_lorenz(self):
        _x = self.read_column_values(self.ui.table_DQ, 4)
        y = self.read_column_values(self.ui.table_DQ, 5)

        button = self.ui.radioButton_Log
        if button.isChecked():
            x = np.log10(_x)
        else:
            x = _x

        initial_guess = [1, 0, 1]  # Initial guess for the parameters [amplitude, center, width]
        params, covariance = curve_fit(self.lorenz, x, y, p0=initial_guess)

        amp_fit, cen_fit, wid_fit = params

        x_fit = np.arange(0, 105.1, 0.1)
        y_fit = self.lorenz(x_fit, *params)

        # Plot the line
        self.t2_dq_graph()
        self.ui.DQ_Widget_2.plot(x_fit, y_fit, pen='r')

    def lorenz(self, x, amp, cen, wid):
        return (amp * (wid**2)) / ((x - cen)**2 + (wid**2))

    def plot_voigt(self):
        def voigt(x, amp_gauss, cen_gauss, wid_gauss, amp_lorenz, cen_lorenz, wid_lorenz, frac):
            gauss_component = self.gaussian(x, amp_gauss, cen_gauss, wid_gauss)
            lorenz_component = self.lorenz(x, amp_lorenz, cen_lorenz, wid_lorenz)
            return frac * lorenz_component + (1 - frac) * gauss_component

        _x = self.read_column_values(self.ui.table_DQ, 4)
        y = self.read_column_values(self.ui.table_DQ, 5)

        button = self.ui.radioButton_Log
        if button.isChecked():
            x = np.log10(_x)
        else:
            x = _x

        initial_guess = [1, 0, 1, 1, 0, 1, 0.5]  # Initial guess for the parameters [amplitude, center, width]
        params, covariance = curve_fit(voigt, x, y, p0=initial_guess)

        #amp_fit, cen_fit, wid_fit = params

        x_fit = np.arange(0, 105.1, 0.1)
        y_fit = voigt(x_fit, *params)

        # Plot the line
        self.t2_dq_graph()
        self.ui.DQ_Widget_2.plot(x_fit, y_fit, pen='r')
        
    def updateGraph(self, index):
        x = self.read_column_values(self.ui.table_SE, 0)

        if index == 0:  # SFC
            y =  self.read_column_values(self.ui.table_SE, 1)
            self.ui.SEWidget.getAxis('left').setLabel("SFC")
            
        elif index == 1:  # M2
            y =  self.read_column_values(self.ui.table_SE, 2)
            self.ui.SEWidget.getAxis('left').setLabel("M2")
            
        elif index == 2:  # T2*
            y =  self.read_column_values(self.ui.table_SE, 3)
            self.ui.SEWidget.getAxis('left').setLabel("T2*")

        elif index == 3:  # Set y
            y =  [0]
            x = [0]
            self.ui.SEWidget.getAxis('left').setLabel("NaN")
            
            
        self.ui.SEWidget.clear()
        self.ui.SEWidget.plot(x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=5)

    def dq_t2_graph(self):
        x = self.read_column_values(self.ui.table_DQ, 0)
        y= self.read_column_values(self.ui.table_DQ, 3)
        self.ui.DQ_Widget_1.clear()
        self.ui.DQ_Widget_1.plot(x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=5)

    def t2_dq_graph(self):
        x = self.read_column_values(self.ui.table_DQ, 4)
        y= self.read_column_values(self.ui.table_DQ, 5)

        graph_DQ_distr = self.ui.DQ_Widget_2
        button = self.ui.radioButton_Log
        if button.isChecked():
            new_x = np.log10(x)
            graph_DQ_distr.getAxis('bottom').setLabel("log(T2*)")
        else:
            new_x = x
            graph_DQ_distr.getAxis('bottom').setLabel("T2*")

        graph_DQ_distr.clear()
        graph_DQ_distr.plot(new_x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=5)

    def read_column_values(self, table, column_index):
        column_values = []
        for row in range(table.rowCount()):
            item = table.item(row, column_index)
            if item is not None:
                column_values.append(float(item.text()))  # Assuming the values are numeric!!!!!
        return column_values
    
    def extract_temperature(self, filename):
        match = re.search(r'_(\d+)_c\.dat', filename)
        if match:
            temperature = match.group(1)
        else:
            temperature = '0'  # Set temperature to '0' if no match is found
        return temperature

    def extract_dq_filtering_time(self, filename):
        match = re.search(r'_(\d+\.\d+)_', filename)
        if match:
            dq_time = match.group(1)
        else:
            dq_time = '0'
        return dq_time

    def calculate_DQ_intensity(self, Time, Amplitude):
        idx_time = np.argmin(np.abs(Time - 4))
        DQ = np.mean(Amplitude[:idx_time])
        return DQ

    def calculate_SFC(self, Amplitude):
        S = np.mean(Amplitude[1:4])
        L = np.mean(Amplitude[50:70])
        SFC = (S-L)/S
        return SFC
    
    def calculate_M2(self, FFT_real, Frequency):
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
        
        return M2
    
    def calculate_T2(self, M2):
        if M2 == 0:
            T2 = 0
        else:
            T2 = np.sqrt(2/M2)
        return T2

    def fill_table_SE(self, SFC_i, M2_i, T2_i, Tem_i, i, table):

        j = i-1

        SFC_i = round(SFC_i, 6)
        M2_i = round(M2_i, 6)
        T2_i = round(T2_i, 6)

        table.setItem(j, 0, QTableWidgetItem(Tem_i))
        table.setItem(j, 1, QTableWidgetItem(str(SFC_i)))
        table.setItem(j, 2, QTableWidgetItem(str(M2_i)))
        table.setItem(j, 3, QTableWidgetItem(str(T2_i)))

    def fill_table_DQ(self, DQ, M2_i, T2_i, dq_time, i, table):
        j = i-1

        DQ = round(DQ, 6)
        M2_i = round(M2_i, 6)
        T2_i = round(T2_i, 6)

        table.setItem(j, 0, QTableWidgetItem(dq_time))
        table.setItem(j, 1, QTableWidgetItem(str(DQ)))
        table.setItem(j, 2, QTableWidgetItem(str(M2_i)))
        table.setItem(j, 3, QTableWidgetItem(str(T2_i)))

    def pre_processing(self, Time_initial, Real_initial, Imaginary_initial):
        # Crop the data time<0
        Time_cropped, Real_cropped, Imaginary_cropped = self.crop_time_zero(Time_initial, Real_initial, Imaginary_initial)
        # Perform time domain phasing
        phase_angle, Re_phased, Im_phased = self.time_domain_phase(Real_cropped, Imaginary_cropped)
        # Calculate amplitude
        Amplitude_cropped = self.calculate_amplitude(Re_phased, Im_phased)
        # Normalize data to max of Amplitude
        Amp, Re, Im = self.normalize(Amplitude_cropped, Re_phased, Im_phased)

        return Time_cropped, Amp, Re, Im  
        
    def calculate_amplitude(self, Real, Imaginary):
        Amp = np.sqrt(Real ** 2 + Imaginary ** 2)
        return Amp
    
    def normalize(self, Amplitude, Real, Imaginary):
        Amplitude_max = np.max(Amplitude)
        Amp = Amplitude/Amplitude_max
        Re = Real/Amplitude_max
        Im = Imaginary/Amplitude_max
        return Amp, Re, Im
    
    def crop_time_zero(self, Time, Real, Imaginary):
        if Time[0] < 0:
            Time_start = 0
            Time_crop_idx = np.where(Time >= Time_start)[0][0]
            Time_cropped = Time[Time_crop_idx:]
            Real_cropped = Real[Time_crop_idx:]
            Imaginary_cropped = Imaginary[Time_crop_idx:]
            print('brit mila')
            return Time_cropped, Real_cropped, Imaginary_cropped
        else:
            return Time, Real, Imaginary

    def time_domain_phase(self, Real, Imaginary):
        delta = np.zeros(360)
        
        for phi in range(360):
            Re_phased = Real * np.cos(np.deg2rad(phi)) - Imaginary * np.sin(np.deg2rad(phi))
            Im_phased = Real * np.sin(np.deg2rad(phi)) + Imaginary * np.cos(np.deg2rad(phi))
            Magnitude_phased = self.calculate_amplitude(Re_phased, Im_phased)
            
            Re_cut = Re_phased[:10]
            Ma_cut = Magnitude_phased[:10]
            
            delta[phi] = np.mean(Ma_cut - Re_cut)
        
        idx = np.argmin(delta)

        Re = Real * np.cos(np.deg2rad(idx)) - Imaginary * np.sin(np.deg2rad(idx))
        Im = Real * np.sin(np.deg2rad(idx)) + Imaginary * np.cos(np.deg2rad(idx))

        
        return idx, Re, Im

    def apodization(self, Time, Amplitude, Real, Imaginary):
        coeffs = np.polyfit(Time, Amplitude, 1)  # Fit an exponential decay function
        c = np.polyval(coeffs, Time)
        d = np.argmin(np.abs(c - 1e-5))
        sigma = Time[d]
        apodization_function = np.exp(-(Time / sigma) ** 4)
        Re_ap = Real * apodization_function
        Im_ap = Imaginary * apodization_function
        return Re_ap, Im_ap

    def add_zeros(self, Time, Real, Imaginary, number_of_points):
        length_diff = number_of_points - len(Time)
        amount_to_add = np.zeros(length_diff)

        Re_zero = np.concatenate((Real, amount_to_add))
        Im_zero = np.concatenate((Imaginary, amount_to_add))

        dt = Time[1] - Time[0]
        Time_to_add = Time[-1] + np.arange(1, length_diff + 1) * dt

        Time = np.concatenate((Time, Time_to_add))
        Fid = np.array(Re_zero + 1j * Im_zero)

        return Time, Fid

    def calculate_frequency_scale(self, Time):
        numberp = len(Time)

        dt = Time[1] - Time[0]
        f_range = 1 / dt
        f_nyquist = f_range / 2
        df = 2 * (f_nyquist / numberp)
        Freq = np.arange(-f_nyquist, f_nyquist + df, df)

        return Freq

    def FFT_handmade(self, Fid, Time, Freq):
        M = len(Time)
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
        
    def simple_baseline_correction(self, FFT):
        twentyperc = int(round(len(FFT) * 0.02))
        Baseline = np.mean(np.real(FFT[:twentyperc]))
        FFT_corrected = FFT - Baseline
        Re = np.real(FFT_corrected)
        Im = np.imag(FFT_corrected)
        Amp = self.calculate_amplitude(Re, Im)
        return Amp, Re, Im
    
    def calculate_apodization(self, Real, Freq):
        # Find sigma at 2% from the max amplitude of the spectra
        Maximum = np.max(np.abs(Real))
        idx_max = np.argmax(np.abs(Real))
        ten_percent = Maximum * 0.02

        b = np.argmin(np.abs(Real[idx_max:] - ten_percent))
        Amplitudes = np.interp(ten_percent, Real, Real)
        sigma_ap = Freq[idx_max + b]

        apodization_function_s = np.exp(-(Freq / sigma_ap) ** 6)

        Real_apod = Real * apodization_function_s
        
        return Real_apod

    def save_data(self):
        parent_folder = os.path.dirname(self.selected_files[0])

        current_tab_index = self.ui.tabWidget.currentIndex()

        SE_temperature = self.ui.SEWidget        
        DQ_points = self.ui.DQ_Widget_1
        DQ_distribution = self.ui.DQ_Widget_2

        SE_table = self.ui.table_SE
        DQ_table = self.ui.table_DQ

        if current_tab_index == 0:
            graph_file_path = os.path.join(parent_folder, 'Result', f"SE_temperature.png")
            table_file_path = os.path.join(parent_folder, 'Result', f"SE_table.csv")
        elif current_tab_index == 1:
            graph_file_path = os.path.join(parent_folder, 'Result', f"DQ_points.png")
            graph_file_path_2 = os.path.join(parent_folder, 'Result', f"DQ_distribution.png")
            table_file_path = os.path.join(parent_folder, 'Result', f"DQ_table.csv")

        if os.path.exists(graph_file_path):
            os.remove(graph_file_path)
        if os.path.exists(graph_file_path_2):
            os.remove(graph_file_path_2)
        if os.path.exists(table_file_path):
            os.remove(table_file_path)

        pg.QtGui.QGuiApplication.processEvents()  # Make sure all events are processed before exporting

        if current_tab_index == 0:
            #Image
            exporter_se = pg.exporters.ImageExporter(SE_temperature.plotItem)
            exporter_se.parameters()['width'] = 1000
            exporter_se.export(graph_file_path)

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

        with open(file_path, 'r') as f:
            lines = f.readlines()
            table.setRowCount(len(lines))
            for row, line in enumerate(lines):
                values = line.strip().split(',')
                for col, value in enumerate(values):
                    item = QTableWidgetItem(value)
                    table.setItem(row, col, item)

        if current_tab_index == 0:
            self.updateGraph()
        elif current_tab_index == 1:
            self.dq_t2_graph()
            self.t2_dq_graph()  
            self.linearization() 



class OpenFilesDialog(QFileDialog, Ui_ChooseFiles):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFileMode(QFileDialog.ExistingFiles)  # Allow selecting multiple files
        self.setNameFilter(str("Data (*.dat *.txt *.csv)"))
        self.setDirectory(str("C:/Mega/NMR/003_Temperature"))
        #self.setDirectory(str("C:/Mega/NMR/003_Temperature/2023_12_21_SE_Temperature_PS35000"))
        
        self.selected_files = []  # Variable to store selected file paths
        self.setupUi(self)
        

    def on_file_selected(self, files):
        self.selected_files = files

class NotificationDialog(QDialog, Ui_Dialog):
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
        


if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = MainWindow()
    widget.show()
    sys.exit(app.exec())
