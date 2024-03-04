# This Python file uses the following encoding: utf-8
import sys, os, re

from PySide6.QtWidgets import QApplication, QMainWindow, QFileDialog, QTableWidgetItem, QDialog
from PySide6.QtCore import QCoreApplication, Signal
from PySide6.QtGui import QColor
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import pyqtgraph as pg
from pyqtgraph.exporters import ImageExporter
from ui_Form import Ui_Analyzer
from ui_ChooseFiles import Ui_ChooseFiles
from ui_Notification import Ui_Dialog
from ui_Error import Ui_Error
from ui_PhasingManual import Ui_Form as Ui_PhasingManual

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_Analyzer()
        self.ui.setupUi(self)

        self.selected_files = []

        # Connect buttons to their respective slots
        self.ui.btn_SelectFiles.clicked.connect(self.clear_list)
        self.ui.btn_SelectFiles.clicked.connect(self.open_select_dialog)
        self.ui.btn_Add.clicked.connect(self.open_select_dialog)
        self.ui.btn_Start.clicked.connect(self.analysis)
        self.ui.btn_Save.clicked.connect(self.save_data)
        self.ui.btn_Load.clicked.connect(self.load_data)
        self.ui.btn_Phasing.clicked.connect(self.open_phasing_manual)
        self.ui.radioButton_Log.clicked.connect(self.plot_fit)

        # Graph setup
        self.setup_graph(self.ui.FFTWidget, "Frequency, MHz", "Amplitude, a.u", "FFT")
        self.setup_graph(self.ui.FidWidget, "Time, ms", "Amplitude", "FID")
        self.setup_graph(self.ui.SEWidget, "Temperature, °C", "", "")
        self.setup_graph(self.ui.DQ_Widget_1, "DQ Filtering Time", "T2*", "")
        self.setup_graph(self.ui.DQ_Widget_2, "", "Norm. DQ Intensity", "")

        # Table setup
        self.setup_table(self.ui.table_SE)
        self.setup_table(self.ui.table_DQ)

        # Connect table signals to slots
        self.ui.table_DQ.currentItemChanged.connect(self.update_dq_graphs)
        self.ui.table_SE.currentItemChanged.connect(self.update_yaxis)

        # Connect combobox signals to slots
        self.ui.comboBox.currentIndexChanged.connect(self.update_yaxis)
        self.ui.comboBox_3.currentIndexChanged.connect(self.update_yaxis)
        self.ui.comboBox_2.currentIndexChanged.connect(self.plot_fit)

        # Connect change events
        self.ui.dq_min.valueChanged.connect(self.update_dq_graphs)
        self.ui.dq_max.valueChanged.connect(self.update_dq_graphs)

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
        
        self.process_file_data(file_path, current_tab_index, i)

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
            
    def highlight_row(self, table, row_selected):
        for col in range(table.columnCount()):
            for row in range(table.rowCount()):
                item = table.item(row, col)
                if item is not None:
                    item.setBackground(QColor(255, 255, 255))
            item_selected = table.item(row_selected-1, col)
            if item_selected is not None:
                item_selected.setBackground(QColor(255, 255, 0))
        
    def setup_graph(self, graph_widget, xlabel="", ylabel="", title=""):
        graph_widget.getAxis('left').setLabel(ylabel)
        graph_widget.getAxis('bottom').setLabel(xlabel)
        graph_widget.setBackground('w')
        graph_widget.setTitle(title)

    def setup_table(self, table_widget):
        column_widths = [70] * 4
        for i, width in enumerate(column_widths):
            table_widget.setColumnWidth(i, width)

    def open_select_dialog(self):
        dlg = OpenFilesDialog(self)
        if dlg.exec():
            fileNames = dlg.selectedFiles()
            self.selected_files.extend(fileNames)
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)

    def open_phasing_manual(self):
        self.phasing_manual_window = PhasingManual()
        self.phasing_manual_window.read_data(self.Frequency, self.Re_spectra, self.Im_spectra)
        self.phasing_manual_window.show()
        

    def disable_buttons(self):
        self.ui.btn_Start.setEnabled(False)
        self.ui.btn_Save.setEnabled(False)
        self.ui.btn_Phasing.setEnabled(False)
        self.ui.dq_min.setEnabled(False)
        self.ui.dq_max.setEnabled(False)
        self.ui.comboBox.setEnabled(False)
        self.ui.comboBox_2.setEnabled(False)
        self.ui.radioButton_Log.setEnabled(False)
        self.ui.checkBox.setEnabled(False)
        self.ui.btn_Add.setEnabled(False)

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
        self.ui.radioButton_Log.setEnabled(True)
        self.ui.checkBox.setEnabled(True)
        self.ui.comboBox_4.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)
    
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

        current_tab_index = self.ui.tabWidget.currentIndex()

        self.ui.table_SE.setRowCount(len(self.selected_files))
        self.ui.table_DQ.setRowCount(len(self.selected_files))

        legend = pg.LegendItem(offset=(300, 10), parent=self.ui.FidWidget.graphicsItem())

        if not self.load_data_and_check_validity((self.selected_files[0]), current_tab_index):
            return

        for i, file_path in enumerate(self.selected_files, start=1):
            filename = os.path.basename(file_path)

            self.disable_buttons()
            self.ui.btn_SelectFiles.setEnabled(False)
            self.ui.btn_Load.setEnabled(False)
            self.ui.radioButton.setEnabled(False)

            self.ui.comboBox_4.addItem(f"{filename}")
            self.ui.comboBox_4.setCurrentIndex(-1)
            
            self.ui.textEdit_6.setText(f"Analysing file {i} out of {len(self.selected_files)}")

            self.process_file_data(file_path, current_tab_index, i)

        self.finalize_analysis(legend, current_tab_index, filename)

    def finalize_analysis(self, legend, current_tab_index, filename):
        self.ui.textEdit_6.setText(f"Finished")
        legend.removeItem('Amplitude')
        legend.removeItem('Re')
        legend.removeItem('Im')
        legend.addItem(self.ui.FidWidget.plotItem.listDataItems()[0], name='Amp')
        legend.addItem(self.ui.FidWidget.plotItem.listDataItems()[1], name='Re')
        legend.addItem(self.ui.FidWidget.plotItem.listDataItems()[2], name='Im')
        self.enable_buttons()
        self.ui.btn_Phasing.setEnabled(True) # TODO: delete after setting phasing manual.
        if current_tab_index == 1:
            self.update_dq_graphs()

    def process_file_data(self, file_path, current_tab_index, i):
        # Read data
        data = np.loadtxt(file_path)
        x, y, z = data[:, 0], data[:, 1], data[:, 2]

        # Read name of filename
        filename = os.path.basename(file_path)
        
        # General mamth procedures and graph updates
        Frequency, Real_apod, Amp = self.general_analysis(x,y,z)

        M2, T2 = self.calculate_M2(Real_apod, Frequency)

        if current_tab_index == 0:
            match = re.search(r'.*_(-?\s*\d+\.?\d*).*.dat', filename)
            temperature = self.extract_info(match)

            SFC = self.calculate_SFC(Amp)
            self.fill_table(self.ui.table_SE, temperature, SFC, M2, T2, i)

            if self.ui.radioButton.isChecked():
                self.save_figures(file_path, temperature)


        elif current_tab_index == 1:
            match = re.search(r'_(\d+\.\d+)_', filename)
            dq_time = self.extract_info(match)

            Amplitude = self.calculate_amplitude(y, z)
            DQ = self.calculate_DQ_intensity(x, Amplitude)            
            self.fill_table(self.ui.table_DQ, dq_time, DQ, M2, T2, i)

            if self.ui.radioButton.isChecked():
                self.save_figures(file_path, dq_time)

    def general_analysis(self, x, y, z):
        # Simple preprocessing (math procedure)
        Time_p, Amp, Re, Im = self.pre_processing(x, y, z)

        # Update FID graph
        self.update_graphs(Time_p, Amp, Re, Im, self.ui.FidWidget)

        Re_ap, Im_ap = self.apodization(Time_p, Amp, Re, Im) #(math procedure)
        Time, Fid = self.add_zeros(Time_p, Re_ap, Im_ap, 16384)  #(math procedure)
        self.Frequency = self.calculate_frequency_scale(Time)  #(math procedure)

        if self.ui.checkBox.isChecked():
            FFT = self.FFT_handmade(Fid, Time, self.Frequency)  #(math procedure)
        else:
            FFT = np.fft.fftshift(np.fft.fft(Fid))
            
        # This condition is never met
        if len(self.Frequency) != len(FFT):
            self.Frequency = np.linspace(self.Frequency[0], self.Frequency[-1], len(FFT))

        self.Amp_spectra, self.Re_spectra, self.Im_spectra = self.simple_baseline_correction(FFT) #(math procedure)
        Real_apod = self.calculate_apodization(self.Re_spectra, self.Frequency) #(math procedure)

        # Update FFT graph
        self.update_graphs(self.Frequency, self.Amp_spectra, self.Re_spectra, self.Im_spectra, self.ui.FFTWidget)

        return self.Frequency, Real_apod, Amp
    
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
            self.ui.table_SE.setHorizontalHeaderLabels(["Temp.", "SFC", "M2", "T2*"])
        elif axis_x == "XS, %":
            self.ui.SEWidget.getAxis('bottom').setLabel("XS, %")
            self.ui.table_SE.setHorizontalHeaderLabels(["XS", "SFC", "M2", "T2*"])
        else:
            return
        
        x = self.read_column_values(self.ui.table_SE, 0)

        text = self.ui.comboBox.currentText()

        if text == "SFC":  
            y =  self.read_column_values(self.ui.table_SE, 1)
            self.ui.SEWidget.getAxis('left').setLabel("SFC")
            
        elif text == "M2":  
            y =  self.read_column_values(self.ui.table_SE, 2)
            self.ui.SEWidget.getAxis('left').setLabel("M2")
            
        elif text == "T2*":  
            y =  self.read_column_values(self.ui.table_SE, 3)
            self.ui.SEWidget.getAxis('left').setLabel("T2*")

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
        y = t2[(dq_time >= time_min) & (dq_time <= time_max)]

        if len(x) <= 1 or len(y) <= 1:
            return

        coeff = np.polyfit(x, y, 1)
        
        Integral = np.trapz(dq)
        DQ_norm = dq/Integral

        self.ui.table_DQ.setColumnCount(5)
        self.ui.table_DQ.setColumnWidth(4,70)
        self.ui.table_DQ.setHorizontalHeaderItem(4, QTableWidgetItem("T2* lin"))

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

    def plot_fit(self):        
        _x = self.read_column_values(self.ui.table_DQ, 4)
        y = self.read_column_values(self.ui.table_DQ, 5)

        button = self.ui.radioButton_Log
        if button.isChecked():
            x = np.log10(_x)
            p = [1, 1, 1, 0]
            b=([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1, np.inf])
        else:
            x = _x
            p = [1, 5, 5, 0]
            b=([0, 0, 0, 0, 0], [ np.inf, np.inf, np.inf, 1, np.inf])

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
        elif text == 'Lorenz':
            params, _ = curve_fit(self.lorenz, x, y, p0=p, bounds=b1)
            y_fit = self.lorenz(x_fit, *params)
            y_r2 = self.lorenz(x, *params)
            cen = params[1]
            fwhm = params[2]
            w = 1
        elif text == 'Pseudo Voigt':
            params, _ = curve_fit(self.voigt, x, y,  bounds = b)
            y_fit = self.voigt(x_fit, *params)
            y_r2 = self.voigt(x, *params)
            cen = params[1]
            fwhm = params[2]
            w = params[3]
        else:
            x_fit = []
            y_fit = []
            return

        #R2
        R = self.calculate_r_squared(y, y_r2)
        # Plot the line
        self.t2_dq_graph()
        self.ui.DQ_Widget_2.plot(x_fit, y_fit, pen='r')

        # Display R2, Xo and FWHM
        self.ui.textEdit_4.setText(f"R\u00B2: {round(R, 4)} \nX\u2080: {round(cen, 4)} \nFWHM: {round(fwhm, 4)} \nFraction (Lorenz): {round(w,2)}")

    def gaussian(self, x, amp, cen, wid, y0):
        return amp * np.exp(-(x - cen)**2 / (2 * wid**2)) + y0

    def lorenz(self, x, amp, cen, wid, y0):
        return (amp * (wid**2)) / ((x - cen)**2 + (wid**2)) + y0
    
    def voigt(self, x, amp, cen, wid, frac, y0):
        lorentzian = amp *(2 * wid) / (np.pi * (4 * (x - cen)**2 + wid**2))
        gaussian = amp *(np.exp((-4 * np.log(2) * (x - cen)**2) / wid**2)) / (wid * np.sqrt(np.pi / (4 * np.log(2))))
        return  (frac * lorentzian + (1 - frac) * gaussian) + y0

    def calculate_r_squared(self, y_true, y_pred):

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
    
    # Working with tables
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

        if current_tab_index == 0:
            table_file_path = os.path.join(parent_folder, 'Result', f"SE_table.csv")
        elif current_tab_index == 1:
            graph_file_path = os.path.join(parent_folder, 'Result', f"DQ_points.png")
            graph_file_path_2 = os.path.join(parent_folder, 'Result', f"DQ_distribution.png")
            table_file_path = os.path.join(parent_folder, 'Result', f"DQ_table.csv")

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
            self.update_dq_graphs()  

        self.enable_buttons()
        self.ui.comboBox_4.setEnabled(False)
        self.ui.checkBox.setEnabled(False)
        self.ui.btn_Phasing.setEnabled(False)

    def save_figures(self, file_path, variable):
        # Set names
        parent_folder = os.path.dirname(file_path)  

        graph_fft = self.ui.FFTWidget
        graph_fid = self.ui.FidWidget

        fft_file_path = (parent_folder + '/Result/' + f"FFT_{variable}.png")
        fid_file_path = (parent_folder + '/Result/' + f"FID_{variable}.png")
        os.makedirs(parent_folder + '/Result/', exist_ok=True)

        pg.QtGui.QGuiApplication.processEvents()

        exporter_fft = pg.exporters.ImageExporter(graph_fft.plotItem)
        exporter_fft.parameters()['width'] = 1000
        exporter_fft.export(fft_file_path)

        exporter_fid = pg.exporters.ImageExporter(graph_fid.plotItem)
        exporter_fid.parameters()['width'] = 1000
        exporter_fid.export(fid_file_path)

    def load_data_and_check_validity(self, file_path, current_tab_index):
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
            not (filename.startswith("DQ") and current_tab_index == 1):
            dialog = NotificationDialog()
            if dialog.exec() == QDialog.Rejected:
                self.ui.btn_SelectFiles.setEnabled(True)
                self.ui.radioButton.setEnabled(True)
                return False
            
        return True

    # Math procedures        
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
            T2 = 0
        else:
            T2 = np.sqrt(2/M2)
        
        return M2, T2

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
        self.selected_files.extend(files)

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
        
class PhasingManual(QDialog):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_PhasingManual()
        self.ui.setupUi(self)

        self.Frequency = None
        self.Re_spectra = None
        self.Im_spectra = None

        graph_phasing = self.ui.PhasingGraph
        graph_phasing.getAxis('bottom').setLabel("Frequency, MHz")
        graph_phasing.getAxis('left').setLabel("Amplitude, a.u.")
        graph_phasing.setBackground('w')
        graph_phasing.setTitle("Phasing")
        self.ui.PhasingGraph.addLegend()

        self.ui.pushButton_2.clicked.connect(self.zero)
        self.ui.verticalSlider_a.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_b.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_c.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_d.valueChanged.connect(self.value_changed)
        self.ui.dial.valueChanged.connect(self.smoothing_changed)
        self.ui.pushButton_3.clicked.connect(self.manual_read)

    def read_data(self, Frequency, Re_spectra, Im_spectra):
        self.Frequency = Frequency
        self.Re = Im_spectra
        self.Im = Re_spectra

        self.zero()

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

        self.process_data()

    #TODO: save and exit
    
    def process_data(self):
        Re_phased = self.calculate_phase()
        self.update_plot(Re_phased)
        self.update_text(Re_phased)
    
    def calculate_phase(self):
        phi = self.a + self.b * self.Frequency + self.c * self.Frequency ** 2 + self.d * self.Frequency ** 3
        Real_phased = self.Re * np.cos(np.deg2rad(phi)) - self.Im * np.sin(np.deg2rad(phi))

        if self.Smooth > 1:
            self.Smooth = int(self.Smooth)
            Real_phased = self.smooth(Real_phased)

        return Real_phased
    
    def smooth(self,y,): #hyinya redo
        new_array = savgol_filter(y, window_length=self.Smooth, polyorder=1)

        return new_array

    def update_plot(self, Re_phased):
        self.ui.PhasingGraph.clear()
        self.ui.PhasingGraph.plot(self.Frequency, self.Re, pen='r', name = 'Original') #Real origin
        self.ui.PhasingGraph.plot(self.Frequency, Re_phased, pen='b', name = 'Phased') #Real phased

    def update_text(self, Re_phased):
        Integral = np.trapz(Re_phased)
        left_mean = np.mean(Re_phased[:100])
        right_mean = np.mean(Re_phased[-100:])
        delta = left_mean - right_mean

        self.ui.Integral.setText(f"Integral: {round(Integral,3)}")
        self.ui.Delta.setText(f"Delta: {round(delta,7)}")


    def value_changed(self):
        self.a = self.ui.verticalSlider_a.value()
        self.b = self.ui.verticalSlider_b.value()
        self.c = self.ui.verticalSlider_c.value()
        self.d = self.ui.verticalSlider_d.value()

        self.ui.Box_a.setValue(self.a)
        self.ui.Box_b.setValue(self.b)
        self.ui.Box_c.setValue(self.c)
        self.ui.Box_d.setValue(self.d)


        if self.Frequency is not None and self.Re is not None and self.Im is not None:
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
