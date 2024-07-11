# This Python file uses the following encoding: utf-8
import sys, os, re
from PySide6.QtWidgets import QApplication, QMainWindow, QFileDialog, QTableWidgetItem, QDialog, QMessageBox, QScrollArea, QWidget
from PySide6.QtCore import QCoreApplication, Signal
from PySide6.QtGui import QColor, QIcon
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from webbrowser import open as open_application
from itertools import islice
import pyqtgraph as pg
from pyqtgraph.exporters import ImageExporter
from ui_Form import Ui_NMR
from ui_Notification import Ui_Note
from ui_PhasingManual import Ui_Phasing as Ui_PhasingManual
import Calculator as Cal # Mathematical procedures

pg.CONFIG_OPTIONS['background'] = 'w'
pg.CONFIG_OPTIONS['foreground'] = 'k'

# Global 
Frequency = []
Re_spectra = []
Im_spectra = []
State_multiple_files = None

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_NMR()
        self.ui.setupUi(self)
        self.setWindowIcon(QIcon('BrandIcon.ico'))

        # Set window geometry
        screen = QApplication.primaryScreen()

        if screen:
            self.showMaximized()

        self.scroll = QScrollArea(self)
        self.scroll.setWidget(self.ui.centralwidget)
        self.scroll.setWidgetResizable(True)
        self.setCentralWidget(self.scroll)
        # Maximize minimize
        self.setWindowFlags(self.windowFlags() | self.windowFlags().WindowMaximizeButtonHint)

        self.selected_files = []
        self.selected_folders = []
        self.selected_files_gly = []
        self.selected_DQfiles = []
        self.selected_T1files = []
        self.selected_FFCfiles = []
        self.selected_DQMQfile = []
        self.dq_t2 = {}
        self.tau_dictionary = {}
        self.ffc_dictionary = {}
        self.tab = None

        self.state()

        # Connect buttons to their respective slots
        self.ui.pushButton_DefaultFolder.clicked.connect(self.default_folder)
        self.ui.commandLinkButton.clicked.connect(self.open_url)

        self.ui.tabWidget.currentChanged.connect(self.state)
        
        self.ui.btn_Save.clicked.connect(self.save_data)
        self.ui.btn_Save_2.clicked.connect(self.save_data)
        self.ui.btn_Load.clicked.connect(self.load_data)
        self.ui.btn_Phasing.clicked.connect(self.open_phasing_manual)
        
        self.ui.btn_SelectFiles.clicked.connect(self.open_select_dialog)
        self.ui.btn_Add.clicked.connect(self.open_select_dialog)
        self.ui.btn_SelectFiles_T1.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.btn_SelectFilesDQMQ.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.btn_SelectFiles_FFC.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.btn_SelectFilesDQ.clicked.connect(self.open_select_comparison_files_dialog)
        
        #self.ui.btn_SelectFiles.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_2.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_3.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable.clicked.connect(self.clear_list)
        self.ui.btn_DeleteRow.clicked.connect(self.delete_row)
        self.ui.btn_DeleteRow_1.clicked.connect(self.delete_row)
        self.ui.btn_DeleteRow_2.clicked.connect(self.delete_row)

        self.ui.btn_Start.clicked.connect(self.analysis)
        self.ui.btn_Launch.clicked.connect(self.launch)
        
        self.ui.pushButton_DQMQ_1.clicked.connect(self.plot_original)
        self.ui.radioButton_Log.clicked.connect(self.plot_fit)
        self.ui.pushButton_DQMQ_4.clicked.connect(self.plot_norm)
        self.ui.pushButton_DQMQ_2.clicked.connect(self.plot_diff)
        self.ui.pushButton_DQMQ_3.clicked.connect(self.plot_nDQ)
        self.ui.btn_Plot1.clicked.connect(self.plot_relaxation_time)
        

        # Graph setup
        self.setup_graph(self.ui.FFTWidget, "Frequency, MHz", "Amplitude, a.u", "FFT")
        self.setup_graph(self.ui.FidWidget, "Time, μs", "Amplitude", "FID")
        self.setup_graph(self.ui.SEWidget, "Temperature, °C", "Choose", "")
        self.setup_graph(self.ui.DQ_Widget_1, "DQ Filtering Time", "T₂*", "")
        self.setup_graph(self.ui.DQ_Widget_2, "X axis", "Norm. DQ Intensity", "")
        self.setup_graph(self.ui.DQ_Widget_3, "DQ Filtering Time", "T₂*", "")
        self.setup_graph(self.ui.DQ_Widget_4, "", "Norm. DQ Intensity", "")
        self.setup_graph(self.ui.DQ_Widget_5, "", "Center", "")
        self.setup_graph(self.ui.DQ_Widget_6, "Name", "FWHM", "")
        self.setup_graph(self.ui.T1_Widget_1, "Time, μs", "Signal", "")
        self.setup_graph(self.ui.T1_Widget_2, "X axis", "τ, μs", "")
        self.setup_graph(self.ui.DQMQ_Widget, "Time", "NMR signal", "")
        self.setup_graph(self.ui.FFC_Widget_1,"Frequency, MHz", "1/T₁", "")
        self.setup_graph(self.ui.FFC_Widget_2, "X axis", "Y Axis", "")

        # Connect table signals to slots
        self.ui.table_DQ.currentItemChanged.connect(self.update_dq_graphs)
        self.ui.table_SE.currentItemChanged.connect(self.update_yaxis)
        self.ui.radioButton_4.clicked.connect(self.calculate_relaxation_time)
        self.ui.radioButton_5.clicked.connect(self.calculate_relaxation_time)
        self.ui.radioButton_6.clicked.connect(self.calculate_relaxation_time)
        self.ui.radioButton_16.clicked.connect(self.calculate_relaxation_time)
        self.ui.radioButton_17.clicked.connect(self.calculate_relaxation_time)
        self.ui.radioButton_3.clicked.connect(self.hide_FFT_progress)
        self.ui.radioButton_2.clicked.connect(self.hide_FFT_progress)
        # Connect combobox signals to slots
        self.ui.comboBox.activated.connect(self.update_yaxis)
        self.ui.comboBox_3.activated.connect(self.update_yaxis)
        self.ui.comboBox_2.activated.connect(self.plot_fit)
        self.ui.comboBox_6.activated.connect(self.calculate_relaxation_time)
        self.ui.comboBox_8.activated.connect(self.calculate_23_model)
        
        # Connect change events
        self.ui.dq_min.valueChanged.connect(self.update_dq_graphs)
        self.ui.dq_max.valueChanged.connect(self.update_dq_graphs)

        self.ui.dq_min_2.valueChanged.connect(self.update_DQ_comparison_plot)
        self.ui.dq_max_2.valueChanged.connect(self.update_DQ_comparison_plot)

        self.ui.radioButton_Log_2.clicked.connect(self.update_DQ_comparison_plot) # this is a bad coding
        self.ui.comboBox_5.activated.connect(self.update_DQ_comparison_plot)
        self.ui.comboBox_4.activated.connect(self.update_file)

        self.ui.dq_min_3.valueChanged.connect(self.plot_diff)
        self.ui.dq_max_3.valueChanged.connect(self.plot_diff)
        self.ui.power.valueChanged.connect(self.plot_diff)

        # Disable buttons initially
        self.disable_buttons()
        self.ui.progressBar.setHidden(True)
        self.ui.textEdit_5.setHidden(True)
        self.ui.textEdit_6.setHidden(True)

    def open_url(self):
        open_application('https://github.com/timesama/M2-Approach/releases')

    def update_file(self):
        i = self.ui.comboBox_4.currentIndex() + 1

        try:
            file_path = self.selected_files[i-1]
        except:
            return
        
        if self.ui.checkBox_2.isChecked():
            try:
                file_path_gly = self.selected_files_gly[i-1]
                self.process_file_data(file_path, file_path_gly, i)
            except:
                return
        else:
            self.process_file_data(file_path, [], i)

            
        # Update general figures
        if self.tab == 'DQ':
            self.highlight_row(self.ui.table_DQ, i) 
            self.update_dq_graphs()
        elif self.tab == 'SE':
            self.highlight_row(self.ui.table_SE, i)
            self.update_yaxis()

        self.ui.btn_Phasing.setEnabled(True)

        #TODO sometime I should add the highlight of the certain point on graph, but I am too lazy
            
    def clear_list(self):
        self.selected_files = []
        self.selected_folders = []
        self.selected_files_gly = []
        self.selected_DQfiles = []
        self.selected_T1files = []
        self.selected_FFCfiles = []
        self.selected_DQMQfile = []
        self.dq_t2 = {}
        self.tau_dictionary = {}
        self.ffc_dictionary = {}

        self.ui.table_SE.setRowCount(0)
        self.ui.table_DQ.setRowCount(0)
        self.ui.table_DQ_2.setRowCount(0)
        self.ui.table_T1.setRowCount(0)
        self.ui.table_FFC_1.setColumnCount(0)

        self.ui.T1_Widget_1.clear()
        self.ui.T1_Widget_2.clear()
        #self.ui.comboBox_6.currentIndexChanged.disconnect(self.calculate_relaxation_time)

        if self.tab == 'T1T2':
            combobox = self.ui.comboBox_6
        elif self.tab == '23Model':
            combobox = self.ui.comboBox_8

        while combobox.count()>0:          
            combobox.removeItem(0)
        
        #self.ui.comboBox_6.currentIndexChanged.connect(self.calculate_relaxation_time)

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
        self.ui.DQ_Widget_6.clear()
        self.ui.SEWidget.clear()
        self.ui.FFTWidget.clear()
        self.ui.FidWidget.clear()

    def delete_row(self):

        if self.tab == 'SE':
            table = self.ui.table_SE
        elif self.tab == 'DQ':
            table = self.ui.table_DQ
        elif self.tab =='T1T2':
            table = self.ui.table_T1
        elif self.tab == '23Model':
            table = self.ui.table_FFC_1
        else:
            return


        row = table.currentRow()
        table.removeRow(row)
            
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
          
    def setup_graph(self, graph_widget, xlabel="", ylabel="", title=""):
        graph_widget.getAxis('left').setLabel(ylabel)
        graph_widget.getAxis('bottom').setLabel(xlabel)
        graph_widget.setTitle(title)

    def setup_table(self, table_widget):
        column_widths = [70] * 4
        for i, width in enumerate(column_widths):
            table_widget.setColumnWidth(i, width)

    def open_select_comparison_files_dialog(self):
        global State_multiple_files

        if self.tab == 'DQMQ':
            State_multiple_files = False
        else:
            State_multiple_files = True
            
        dlg = OpenFilesDialog(self)
        if dlg.exec():
            if self.tab == 'DQ_Temp':
                DQfileNames = dlg.selectedFiles()
                self.selected_DQfiles.extend(DQfileNames)
                self.update_DQ_comparison()
            elif self.tab == 'T1T2':
                while self.ui.comboBox_6.count()>0:
                    self.ui.comboBox_6.removeItem(0)
                T1fileNames = dlg.selectedFiles()
                self.selected_T1files.extend(T1fileNames)
                self.update_T12_table()
            elif self.tab == 'DQMQ':
                self.selected_DQMQfile = dlg.selectedFiles()
                self.dq_mq_analysis()
            elif self.tab == '23Model':
                while self.ui.comboBox_8.count()>0:
                    self.ui.comboBox_8.removeItem(0)
                FFCfileNames = dlg.selectedFiles()
                self.selected_FFCfiles.extend(FFCfileNames)
                self.update_FFC_table()

    def open_select_dialog(self):
        global State_multiple_files
        State_multiple_files = True
        dlg = OpenFilesDialog(self)
        if dlg.exec():
            self.selected_files = []
            fileNames = dlg.selectedFiles()
            self.selected_files.extend(fileNames)
            self.ui.btn_Start.setEnabled(True)
            self.ui.btn_Start.setStyleSheet("background-color: green")
            self.ui.btn_Add.setEnabled(True)
    
    def open_select_dialog_glycerol(self):
        dlg = OpenFilesDialog(self)
        dlg.setWindowTitle("Select Reference Files")
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

    def state(self):
        current_tab_index =  self.ui.tabWidget.currentIndex()

        if current_tab_index == 0:
            self.tab = 'SE'
        elif current_tab_index == 1:
            self.tab = 'DQ'
        elif current_tab_index == 2:
            self.tab = 'DQ_Temp'
        elif current_tab_index == 3:
            self.tab = 'T1T2'
        elif current_tab_index == 4:
            self.tab = 'DQMQ'
        elif current_tab_index == 5:
            self.tab = '23Model'
        elif current_tab_index == 6:
            self.tab = 'Extra'

        print(f"state: {self.tab}")

        if not (self.tab == 'SE' or self.tab == 'DQ'):
            self.ui.BOX_up.setHidden(True)
        else:
            self.ui.BOX_up.setHidden(False)

    def disable_buttons(self):
        self.ui.btn_Start.setEnabled(False)
        self.ui.btn_Save.setEnabled(False)
        self.ui.btn_Phasing.setEnabled(False)
        self.ui.dq_min.setEnabled(False)
        self.ui.dq_max.setEnabled(False)
        self.ui.comboBox.setEnabled(False)
        self.ui.comboBox_2.setEnabled(False)
        self.ui.radioButton_Log.setEnabled(False)
        self.ui.btn_Add.setEnabled(False)
        self.ui.radioButton_Log_2.setEnabled(False)
        self.ui.comboBox_5.setEnabled(False)
        self.ui.dq_min_2.setEnabled(False)
        self.ui.dq_max_2.setEnabled(False)
        self.ui.btn_Launch.setEnabled(False)
        self.ui.btn_Plot1.setEnabled(False)
        self.ui.pushButton_DQMQ_1.setEnabled(False)
        self.ui.pushButton_DQMQ_2.setEnabled(False)
        self.ui.pushButton_DQMQ_3.setEnabled(False)
        self.ui.pushButton_DQMQ_4.setEnabled(False)

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
        self.ui.comboBox_4.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)
    
    def hide_FFT_progress(self):
        if self.ui.radioButton_3.isChecked():
            self.ui.progressBar.setHidden(True)
            self.ui.textEdit_5.setHidden(True)
            self.ui.textEdit_6.setHidden(True)
        else:
            self.ui.progressBar.setHidden(False)
            self.ui.textEdit_5.setHidden(False)
            self.ui.textEdit_6.setHidden(False)
    
    def default_folder(self):
        folder_path = str(QFileDialog.getExistingDirectory(self, "Select Default Directory"))

        if os.path.exists("selected_folder.txt"):
            os.remove("selected_folder.txt")

        if folder_path:
            # Write the selected folder's full path to a text file
            with open("selected_folder.txt", "w") as file:
                file.write(folder_path)

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


        legend = pg.LegendItem(offset=(300, 10), parent=self.ui.FidWidget.graphicsItem())

        if not self.load_data_and_check_validity((self.selected_files[0])):
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
                if len(self.selected_files_gly) < len(self.selected_files):
                    QMessageBox.warning(self, "Invalid Data", f"The amount of reference files is not the same as sample files. Adding glycerol files automatically.", QMessageBox.Ok)
                    last_file = self.selected_files_gly[-1]
                    num_to_add = len(self.selected_files) - len(self.selected_files_gly)
                    self.selected_files_gly.extend([last_file] * num_to_add)
                
                file_path_gly = self.selected_files_gly[0]
                try:
                    self.process_file_data(file_path, file_path_gly, i)
                except:
                    self.analysis_error(file_path)
                    return
            else:
                try:
                    self.process_file_data(file_path, [], i)
                except:
                    self.analysis_error(file_path)
                    return
                    
        self.finalize_analysis(legend)
        self.ui.btn_Start.setStyleSheet("background-color: none")

    def analysis_error(self, file_path):
        QMessageBox.warning(self, "Invalid Data", f"Error in {file_path}. Restarting analysis", QMessageBox.Ok)
        self.selected_files.remove(file_path)
        self.ui.btn_SelectFiles.setEnabled(True)
        self.analysis()

    def finalize_analysis(self, legend):
        self.ui.textEdit_6.setText(f"Finished")
        legend.removeItem('Amplitude')
        legend.removeItem('Re')
        legend.removeItem('Im')
        legend.addItem(self.ui.FidWidget.plotItem.listDataItems()[0], name='Amp')
        legend.addItem(self.ui.FidWidget.plotItem.listDataItems()[1], name='Re')
        legend.addItem(self.ui.FidWidget.plotItem.listDataItems()[2], name='Im')
        self.enable_buttons()

        if self.tab == 'DQ':
            self.update_dq_graphs()

    def process_file_data(self, file_path, file_path_gly, i):
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
            Time_r, Re_r, Im_r = Cal.analysis_time_domain(file_path_gly)
            Time_s, Re_s, Im_s = Cal.analysis_time_domain(file_path)
            Time, Re, Im = Cal.long_component(Time_s, Time_r, Re_s, Re_r, Im_s, Im_r) 
        else:
            Time, Re, Im = Cal.analysis_time_domain(file_path)
        
        Amp = Cal.calculate_amplitude(Re, Im)
        self.update_graphs(Time, Amp, Re, Im, self.ui.FidWidget)

        Time_fid, Fid =  Cal.final_analysis_time_domain(Time, Re, Im)

        Frequency = Cal.calculate_frequency_scale(Time_fid)
        if self.ui.radioButton_2.isChecked():
            FFT = self.FFT_handmade(Fid, Time_fid, Frequency)  #(math procedure)
        else:
            FFT = np.fft.fftshift(np.fft.fft(Fid))

        # 8. Simple baseline
        Amp_spectra, Re_spectra, Im_spectra = Cal.simple_baseline_correction(FFT)
        # 9. Cal.apodization
        Real_apod = Cal.calculate_apodization(Re_spectra, Frequency)

        # Update FFT graphs
        self.update_graphs(Frequency, Amp_spectra, Re_spectra, Im_spectra, self.ui.FFTWidget)

        if self.ui.comboBox_4.currentIndex() == -1:
            M2, T2 = Cal.calculate_M2(Real_apod, Frequency)

            if self.tab == 'SE':
                match = re.search(r'.*_(-?\s*\d+\.?\d*).*.dat', filename)
                temperature = self.extract_info(match)

                SFC = Cal.calculate_SFC(Amp)
                self.ui.table_SE.setRowCount(len(self.selected_files))
                self.fill_table(self.ui.table_SE, temperature, SFC, M2, T2, i)

                self.ui.table_SE.setItem(i-1, 4, QTableWidgetItem(filename))

                if self.ui.radioButton.isChecked():
                    self.save_figures(file_path, filename)

            elif self.tab == 'DQ':
                match = re.search(r'_(\d+\.\d+)_', filename)
                dq_time = self.extract_info(match)

                Amplitude = Cal.calculate_amplitude(y, z)
                DQ = Cal.calculate_DQ_intensity(x, Amplitude)     
                self.ui.table_DQ.setRowCount(len(self.selected_files))       
                self.fill_table(self.ui.table_DQ, dq_time, DQ, M2, T2, i)

                if self.ui.radioButton.isChecked():
                    self.save_figures(file_path, dq_time)
        else:
            pass
     
    def after_phasing(self):
        global Frequency, Re_spectra, Im_spectra

        i = self.ui.comboBox_4.currentIndex()

        Real_apod   = Cal.calculate_apodization(Re_spectra, Frequency) #(math procedure)
        Amp_spectra = Cal.calculate_amplitude(Re_spectra, Im_spectra)

        # Update FFT graph
        self.update_graphs(Frequency, Amp_spectra, Re_spectra, Im_spectra, self.ui.FFTWidget)

        M2, T2 = Cal.calculate_M2(Real_apod, Frequency)
        M2_r = round(M2, 6)
        T2_r = round(T2, 6)

        if self.tab == 'SE':
            table = self.ui.table_SE
        elif self.tab == 'DQ':
            table = self.ui.table_DQ


        table.setItem(i, 2, QTableWidgetItem(str(M2_r)))
        table.setItem(i, 3, QTableWidgetItem(str(T2_r)))

        if self.tab == 'SE':
            self.update_yaxis()
        elif self.tab == 'DQ':
            self.update_dq_graphs()

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
            self.ui.table_SE.setHorizontalHeaderLabels(["Temp.", "SC", "M₂", "T₂*"])
        elif axis_x == "XS, %":
            self.ui.SEWidget.getAxis('bottom').setLabel("XS, %")
            self.ui.table_SE.setHorizontalHeaderLabels(["XS", "SC", "M₂", "T₂*"])
        else:
            return
        
        x = self.read_column_values(self.ui.table_SE, 0)

        text = self.ui.comboBox.currentText()

        if text == "SC":  
            y =  self.read_column_values(self.ui.table_SE, 1)
            self.ui.SEWidget.getAxis('left').setLabel("SC")
            
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
        
        try:
            x_fit = np.arange(0, np.max(x) + 0.001, 0.01)
        except:
            x_fit = np.arange(0, 100 + 0.001, 0.01)

        b1=([0, 0, 0, -10], [np.inf, np.inf, np.inf, np.inf])

        if text == 'Gauss':
            params, _ = curve_fit(Cal.gaussian, x, y, p0=p, bounds=b1)
            y_fit = Cal.gaussian(x_fit, *params)
            y_r2 = Cal.gaussian(x, *params)
            cen = params[1]
            fwhm = params[2]
            w = 0
            button.setEnabled(True)
        elif text == 'Lorenz':
            params, _ = curve_fit(Cal.lorenz, x, y, p0=p, bounds=b1)
            y_fit = Cal.lorenz(x_fit, *params)
            y_r2 = Cal.lorenz(x, *params)
            cen = params[1]
            fwhm = params[2]
            w = 1
            button.setEnabled(True)
        elif text == 'Pseudo Voigt':
            params, _ = curve_fit(Cal.voigt, x, y,  bounds = b)
            y_fit = Cal.voigt(x_fit, *params)
            y_r2 = Cal.voigt(x, *params)
            cen = params[1]
            fwhm = params[2]
            w = params[3]
            button.setEnabled(True)
        else:
            button.setEnabled(False)
            return

        #R2
        R = Cal.calculate_r_squared(y, y_r2)
        # Plot the line
        self.t2_dq_graph()
        self.ui.DQ_Widget_2.plot(x_fit, y_fit, pen='r')

        # Display R2, Xo and FWHM
        self.ui.textEdit_4.setText(f"R\u00B2: {round(R, 4)} \nX\u2080: {round(cen, 4)} \nFWHM: {round(fwhm, 4)} \nFraction (Lorenz): {round(w,2)}")

    # DQ MQ section
    def dq_mq_analysis(self):
        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        table = self.ui.table_DQMQ
        file_path = self.selected_DQMQfile[0]

        Time, DQ, Ref = Cal.read_data(file_path, 1)
        num_rows  = len(Time)
        table.setRowCount(num_rows)

        for row in range(num_rows):
            table.setItem(row, 0, QTableWidgetItem(str(Time[row])))
            table.setItem(row, 1, QTableWidgetItem(str(DQ[row])))
            table.setItem(row, 2, QTableWidgetItem(str(Ref[row])))

        table.resizeColumnsToContents() # TODO ADD THIS EVERYWHERE!!!!
        self.ui.pushButton_DQMQ_1.setEnabled(True)
        self.ui.pushButton_DQMQ_2.setEnabled(False)
        self.ui.pushButton_DQMQ_3.setEnabled(False)
        self.ui.pushButton_DQMQ_4.setEnabled(True)

    def plot_original(self):
        file_path = self.selected_DQMQfile[0]
        # fit_from = self.ui.dq_min_3.value()
        # fit_to = self.ui.dq_max_3.value()
        # p = self.ui.power.value()

        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        legend = figure.addLegend()

        Time, DQ, Ref = Cal.read_data(file_path, 1)

        figure.plot(Time, DQ, pen='r', name = 'DQ')
        figure.plot(Time, Ref, pen='b', name = 'Ref')

    def plot_norm(self):
        file_path = self.selected_DQMQfile[0]
        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        legend = figure.addLegend()

        Time, DQ_norm, Ref_norm, _, _, _, _, _, _ = Cal.dqmq(file_path, 40, 100, 1)

        figure.plot(Time, DQ_norm, pen='r', name = 'DQ')
        figure.plot(Time, Ref_norm, pen='b', name = 'Ref')

        self.ui.pushButton_DQMQ_2.setEnabled(True)
        self.ui.pushButton_DQMQ_3.setEnabled(True)     
        self.ui.dq_min_3.setEnabled(True)
        self.ui.dq_max_3.setEnabled(True)
        self.ui.power.setEnabled(True)
    
    def plot_diff(self):
        file_path = self.selected_DQMQfile[0]
        fit_from = self.ui.dq_min_3.value()
        fit_to = self.ui.dq_max_3.value()
        p = self.ui.power.value()

        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        legend = figure.addLegend()

        Time, DQ_norm, Ref_norm, Diff, _, _, _, _, fitted_curve = Cal.dqmq(file_path, fit_from, fit_to, p)

        figure.plot(Time, DQ_norm, pen='r', name = 'DQ')
        figure.plot(Time, Ref_norm, pen='b', name = 'Ref')
        figure.plot(Time, Diff, pen='k', name = 'Diff')
        figure.plot(Time, fitted_curve, pen='m', name = 'fitting')

        self.ui.pushButton_DQMQ_2.setEnabled(True)
        self.ui.pushButton_DQMQ_3.setEnabled(True)    

    def plot_nDQ(self):
        file_path = self.selected_DQMQfile[0]
        fit_from = self.ui.dq_min_3.value()
        fit_to = self.ui.dq_max_3.value()
        p = self.ui.power.value()
        table = self.ui.table_DQMQ

        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        legend = figure.addLegend()

        Time, _, _, _, DQ_normal, MQ_normal, Time0, nDQ, _ = Cal.dqmq(file_path, fit_from, fit_to, p)

        figure.plot(Time, DQ_normal, pen='r', name = 'DQ')
        figure.plot(Time, MQ_normal, pen='b', name = 'Ref')
        figure.plot(Time0, nDQ, pen='k', symbol='o', symbolPen='k', symbolSize=5, name='nDQ')

        num_rows  = len(Time)
        for row in range(num_rows):
            table.setItem(row, 3, QTableWidgetItem(str(round(nDQ[row+1],4))))
        table.resizeColumnsToContents()
        
    # Relaxation time section
    def update_T12_table(self):
        selected_files = self.selected_T1files
        table = self.ui.table_T1
        combobox = self.ui.comboBox_6
        pattern = r'(T1_.*\.dat|T2_.*\.dat)'
        dictionary = self.tau_dictionary

        file_name = []
        x_axis = []
        
        for file in selected_files:
            current_file = [os.path.basename(file)]
            try:
                x_axis += [re.search(pattern,file).group(1)]
            except:
                x_axis += ['0']

            Time = []
            Signal = []

            try:
                # Read files as I create them from excel in spintrack
                with open(file, "r") as data:
                    lines = [line.replace('\t\t\t', '').rstrip('\n') for line in data if not (line.rstrip('\n').endswith('\t\t\t\t'))]

                for line in lines[1:]:  # Skip the first line !!!
                    parts = line.split('\t')
                    time_value = float(parts[0])
                    signal_value = float(parts[1])
                    Time.append(time_value)
                    Signal.append(signal_value)
            except:
                try:
                    # read files regularly
                    data = np.loadtxt(file)
                    Time, Signal = data[:, 0], data[:, 1]
                except Exception:
                    QMessageBox.warning(self, "Invalid Data", f"I couldn't read {current_file}, removing file from the table and file list.", QMessageBox.Ok)
                    for file_to_delete in selected_files:
                        if file_to_delete == file:
                            selected_files.remove(file)

            dictionary[file] = {"X Axis": [], "Time": [], "Signal": []}
            dictionary[file]["X Axis"].append(x_axis)
            dictionary[file]["Time"].extend(Time)
            dictionary[file]["Signal"].extend(Signal)

            
            table.setRowCount(len(selected_files))
            self.ui.btn_Plot1.setEnabled(True)            
        
        for row, file in zip(range(table.rowCount()), selected_files):
            Folder = QTableWidgetItem(file)
            file_name = os.path.basename(file)
            try:
                x_axis = re.search(pattern,file).group(1)
            except:
                x_axis = row

            Filename = QTableWidgetItem(file_name)
            Temp = QTableWidgetItem(x_axis)
            table.setItem(row, 0, Folder)
            table.setItem(row, 1, Filename)
            table.setItem(row, 2, Temp)

            combobox.addItem(f"{file_name}")
            combobox.setCurrentIndex(-1)

    def calculate_relaxation_time(self):
        table = self.ui.table_T1
        figure = self.ui.T1_Widget_1
        selected_file_idx = self.ui.comboBox_6.currentIndex()
        dictionary = self.tau_dictionary

        if self.ui.radioButton_16.isChecked():
        # milisec units
            denominator = 1
        else:
        # microsec units
            denominator = 1000

        if selected_file_idx == -1:
            return
        
        value_from_row = table.item(selected_file_idx, 0).text()
        Time = np.array(dictionary[value_from_row]['Time'])/denominator
        Signal = np.array(dictionary[value_from_row]['Signal'])

        Time_fit = np.arange(min(Time), max(Time) + 1, 1)

        try: 
            if self.ui.radioButton_4.isChecked():
                order = 1
                Time_fit, fitted_curve, tau1, _, _, R2, decrease_order = Cal.fit_exponent(Time, Signal, order)
                tau2 = 0
                tau3 = 0

            elif self.ui.radioButton_5.isChecked():
                order = 2
                Time_fit, fitted_curve, tau1, tau2, _, R2, decrease_order = Cal.fit_exponent(Time, Signal, order)
                tau3 = 0

            
            else: 
                order = 3
                Time_fit, fitted_curve, tau1, tau2, tau3, R2, decrease_order = Cal.fit_exponent(Time, Signal, order)

        except:
            QMessageBox.warning(self, "No covariance", f"I am sorry, I couldn't fit with {order} exponents. Fitting with one.", QMessageBox.Ok)
            order = 1
            Time_fit, fitted_curve, tau1, _, _, R2, decrease_order = Cal.fit_exponent(Time, Signal, order)
            tau2 = 0
            tau3 = 0
        
        self.ui.textEdit_error.setText(f"R² {R2}") 
        
        
        tau_str = str(tau1)
        tau_str2 = str(tau2)
        tau_str3 = str(tau3)
        item = QTableWidgetItem(tau_str)
        item2 = QTableWidgetItem(tau_str2)
        item3 = QTableWidgetItem(tau_str3)
        table.setItem(selected_file_idx,3,item)
        table.setItem(selected_file_idx,4,item2)
        table.setItem(selected_file_idx,5,item3)
        
        
        figure.clear()
        figure.plot(Time, Signal, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=5)
        figure.plot(Time_fit, fitted_curve, pen='b')

    def plot_relaxation_time(self):

        table = self.ui.table_T1
        graph = self.ui.T1_Widget_2

        if self.ui.radioButton_10.isChecked():
            column = 3
        elif self.ui.radioButton_11.isChecked():
            column = 4
        elif self.ui.radioButton_12.isChecked():
            column = 5

        graph.clear()
        if table.rowCount() < 1:
            return

        x_axis = []
        relaxation_time =[]
        number = 1
        for row in range(table.rowCount()):
            try:
                C = float(table.item(row, 2).text())
                x_axis.append(C)
            except:
                x_axis.append(number)
            
            try:
                T = float(table.item(row, column).text())
                relaxation_time.append(T)
            except:
                relaxation_time.append(0)
            number = number + 1
                    

        graph.plot(x_axis, relaxation_time, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=5)

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

            pattern = r'DQ_table_(\d+)(?!\([^)]*\))(?![A-Za-z])'
            try:
                default_name = QTableWidgetItem(str(re.search(pattern, filename).group(1)))
            except:
                default_name = QTableWidgetItem(str(row+1))

            table.setItem(row, 0, item)
            table.setItem(row, 1, default_name)
        
        table.resizeColumnsToContents()

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
        self.ui.DQ_Widget_6.clear()
        self.ui.DQ_Widget_5.clear()
        self.ui.DQ_Widget_4.clear()
        self.ui.DQ_Widget_3.clear()

        if legend1 is not None:
            #legend.clear()
            legend1.clear()
            #self.ui.DQ_Widget_3.addLegend()
            self.ui.DQ_Widget_4.addLegend()
            #legend.setPen((0, 0, 0))  
            legend1.setPen((0, 0, 0))

        # This is a mess :(

        center = []
        weight = []
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
                params, _ = curve_fit(Cal.gaussian, t2_lin, dq_norm, p0=p, bounds=b1)
                y_fit = Cal.gaussian(dq_fit, *params)
                y_r2 = Cal.gaussian(t2_lin, *params)
                cen = params[1]
                center.append(cen)
                fwhm = params[2]
                weight.append(fwhm)
                w = 0
            elif text == 'Lorenz':
                params, _ = curve_fit(Cal.lorenz, t2_lin, dq_norm, p0=p, bounds=b1)
                y_fit = Cal.lorenz(dq_fit, *params)
                y_r2 = Cal.lorenz(t2_lin, *params)
                cen = params[1]
                center.append(cen)
                fwhm = params[2]
                weight.append(fwhm)
                w = 1
            elif text == 'Pseudo Voigt':
                params, _ = curve_fit(Cal.voigt, t2_lin, dq_norm,  bounds = b)
                y_fit = Cal.voigt(dq_fit, *params)
                y_r2 = Cal.voigt(t2_lin, *params)
                cen = params[1]
                center.append(cen)
                fwhm = params[2]
                weight.append(fwhm)
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


        self.ui.DQ_Widget_5.plot(comparison_par, center, pen='r', symbolPen=None, symbol='o', symbolBrush='r')
        self.ui.DQ_Widget_6.plot(comparison_par, weight, pen='b', symbolPen=None, symbol='o', symbolBrush='b')

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

        table.resizeColumnsToContents()
    
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
        try:
            parent_folder = os.path.dirname(self.selected_files[0])
        except:
            parent_folder = os.path.dirname(self.selected_DQMQfile[0])
     
        DQ_points = self.ui.DQ_Widget_1
        DQ_distribution = self.ui.DQ_Widget_2

        SE_table = self.ui.table_SE
        DQ_table = self.ui.table_DQ
        DQMQ_table = self.ui.table_DQMQ

        os.makedirs(parent_folder + '/Result/', exist_ok=True)
        basefolder = os.path.dirname(parent_folder)
        os.makedirs(basefolder + '/DQ_table/', exist_ok=True)
        os.makedirs(basefolder + '/DQMQ_Result/', exist_ok=True)
        if self.tab == 'SE':
            table_file_path = os.path.join(parent_folder, 'Result', f"SE_table.csv")
        elif self.tab == 'DQ':
            graph_file_path = os.path.join(parent_folder, 'Result', f"DQ_points.png")
            graph_file_path_2 = os.path.join(parent_folder, 'Result', f"DQ_distribution.png")
            table_file_path = os.path.join(basefolder, 'DQ_table', f"DQ_table_{os.path.basename(parent_folder)}.csv")
        elif self.tab == 'DQMQ':
            table_file_path = os.path.join(basefolder,'DQMQ_Result', f"DQMQ_Result_{os.path.basename(parent_folder)}.csv")
        pg.QtGui.QGuiApplication.processEvents()  # Make sure all events are processed before exporting

            # Check if the file already exists, if so, append an index to the filename
        if os.path.exists(table_file_path):
            index = 1
            while os.path.exists(f"{table_file_path[:-4]} ({index}).csv"):
                index += 1
            table_file_path = f"{table_file_path[:-4]} ({index}).csv"

        if self.tab == 'SE':
            #Table
            self.save_table_to_csv(table_file_path, SE_table)

        elif self.tab == 'DQ':
            #Image
            exporter_dq1 = pg.exporters.ImageExporter(DQ_points.plotItem)
            exporter_dq1.parameters()['width'] = 1000
            exporter_dq1.export(graph_file_path)

            exporter_dq2 = pg.exporters.ImageExporter(DQ_distribution.plotItem)
            exporter_dq2.parameters()['width'] = 1000
            exporter_dq2.export(graph_file_path_2)

            #Table
            self.save_table_to_csv(table_file_path, DQ_table)
        elif self.tab == 'DQMQ':
            self.save_table_to_csv(table_file_path, DQMQ_table)

        QMessageBox.information(self, "Data Saved", f"The data have been saved to {parent_folder}", QMessageBox.Ok)

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
        if self.tab == 'SE':
            table = self.ui.table_SE
        elif self.tab == 'DQ':
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

        if self.tab == 'SE':
            self.update_yaxis()
        elif self.tab == 'DQ':
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

    def load_data_and_check_validity(self, file_path):
        try: 
            data = np.loadtxt(file_path)
            filename = os.path.basename(file_path)
            # Check that there are 3 columns
            if data.shape[1] != 3:
                QMessageBox.warning(self, "Invalid Data Format", f"I can't read the {filename} file, it should have 3 columns exactly, deleting it.", QMessageBox.Ok)
                self.ui.btn_SelectFiles.setEnabled(True)
                self.ui.radioButton.setEnabled(True)
                self.selected_files.clear()
                return False
            return True

        except:
            QMessageBox.warning(self, "Invalid Data", f"I can't read the {filename} file, deleting it.", QMessageBox.Ok)
            self.selected_files.remove(file_path)

            if self.selected_files == []:
                self.disable_buttons()
                return False
            else:
                return False

    # FFC analysis
    def update_FFC_table(self):
        selected_files = self.selected_FFCfiles
        table = self.ui.table_FFC_1
        combobox = self.ui.comboBox_8

        dictionary = self.ffc_dictionary

        file_name = []
        
        for file in selected_files:
            current_file = [os.path.basename(file)]

            Omega = []
            Rate = []

            try:
                with open(file, "r") as data:
                    for line in islice(data, 4, None):
                        columns = line.split()
                        Omega.append(columns[0])
                        Rate.append(columns[2])
            except Exception:
                    QMessageBox.warning(self, "Invalid Data", f"I couldn't read {current_file}, removing file from the table and file list.", QMessageBox.Ok)
                    for file_to_delete in selected_files:
                        if file_to_delete == file:
                            selected_files.remove(file)

            dictionary[file] = {"File": [], "Freq": [], "Rate": []}
            dictionary[file]["File"].append(current_file)
            dictionary[file]["Freq"].extend(Omega)
            dictionary[file]["Rate"].extend(Rate)

            
            table.setColumnCount(len(selected_files))
            self.ui.btn_Plot1_2.setEnabled(True)            
        
        for column, file in zip(range(table.columnCount()), selected_files):
            Folder = QTableWidgetItem(file)
            file_name = os.path.basename(file)

            Filename = QTableWidgetItem(file_name)

            table.setItem(0, column, Folder)
            table.setItem(1, column, Filename)

            combobox.addItem(f"{file_name}")
            combobox.setCurrentIndex(-1)

    def calculate_23_model(self):
        table = self.ui.table_FFC_1
        figure = self.ui.FFC_Widget_1
        selected_file_idx = self.ui.comboBox_8.currentIndex()
        dictionary = self.ffc_dictionary

        if selected_file_idx == -1:
            return
        
        value_from_row = table.item(0, selected_file_idx).text()
        Omega = np.array(dictionary[value_from_row]['Freq'], dtype=float)
        #Omega = Omeg * 10**6

        Rate = np.array(dictionary[value_from_row]['Rate'], dtype=float)

        if self.ui.checkBox_3.isChecked():
            fixed_CDD = self.ui.doubleSpinBox.value()
        else:
            fixed_CDD = None

        Omega_fit, fitted_curve, popt, R2= Cal.fit_model(Omega, Rate, fixed_CDD)
        
        self.ui.textEdit_error_2.setText(f"R² {R2}") 
        print(f'filename: {value_from_row}\nCDD: {popt[0]}\ntauc: {popt[1]}\nA: {popt[2]}\nCtrans: {popt[3]}\ntautrans: {popt[4]}\ntaures: {popt[5]}\n')


        CDD         = QTableWidgetItem(str("{:.5e}".format(popt[0])))
        tau_c       = QTableWidgetItem(str("{:.5e}".format(popt[1])))
        A           = QTableWidgetItem(str("{:.5e}".format(popt[2])))
        C_trans     = QTableWidgetItem(str("{:.5e}".format(popt[3])))
        tau_trans   = QTableWidgetItem(str("{:.5e}".format(popt[4])))
        tau_res     = QTableWidgetItem(str("{:.5e}".format(popt[5])))

        table.setItem(3, selected_file_idx,CDD)
        table.setItem(5, selected_file_idx,tau_c)
        table.setItem(2, selected_file_idx,A)
        table.setItem(4, selected_file_idx,C_trans)
        table.setItem(6, selected_file_idx,tau_trans)
        table.setItem(7, selected_file_idx,tau_res)

        figure.clear()
        figure.plot(Omega, Rate, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=5)
        figure.plot(Omega_fit, fitted_curve, pen='b')

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
        global State_multiple_files
        super().__init__(parent)

        if State_multiple_files:
            self.setFileMode(QFileDialog.ExistingFiles)  # Allow selecting multiple files
        else:
            pass

        self.setNameFilter(str("Data (*.dat *.txt *.csv *.sef)"))

        initial_directory_file = "selected_folder.txt"
        if os.path.exists(initial_directory_file):
            with open(initial_directory_file, 'r') as file:
                initial_directory = file.read().strip()
            if os.path.isdir(initial_directory):
                self.setDirectory(initial_directory)
            else:
                # Fallback if the content of the file is not a valid directory
                exe_dir = os.path.dirname(sys.argv[0])
                self.setDirectory(exe_dir)
        else:
            exe_dir = os.path.dirname(sys.argv[0])
            self.setDirectory(exe_dir)

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

        self.set_zero(self.ui.verticalSlider_a, self.ui.Box_a)
        self.set_zero(self.ui.verticalSlider_b, self.ui.Box_b)
        self.set_zero(self.ui.verticalSlider_c, self.ui.Box_c)
        self.set_zero(self.ui.verticalSlider_d, self.ui.Box_d)

        self.ui.dial.setValue(0)
        self.ui.Box_smooth.setValue(0)

        if self.Real_freq_phased is not None:
            self.process_data()

    def set_zero(self,slider,box):
        slider.setValue(0)
        box.setValue(0)
    
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

        self.check_borders(self.ui.verticalSlider_a)
        self.check_borders(self.ui.verticalSlider_b)
        self.check_borders(self.ui.verticalSlider_c)
        self.check_borders(self.ui.verticalSlider_d)


        if Frequency is not None and Re_spectra is not None and Im_spectra is not None:
            self.process_data()
        else:
            return


    def check_borders(self, slider):
        max = slider.maximum()
        min = slider.minimum()

        val = slider.value()

        if val + 10 > max:
            slider.setMaximum(val+10)
        elif val - 10 < min:
            slider.setMinimum(val-10)


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
