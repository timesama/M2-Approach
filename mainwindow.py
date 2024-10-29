# This Python file uses the following encoding: utf-8
import sys, os, re, csv, requests
from PySide6.QtWidgets import QApplication, QMainWindow, QFileDialog, QTableWidgetItem, QInputDialog, QDialog, QMessageBox, QScrollArea
from PySide6.QtCore import QCoreApplication, Signal
from PySide6.QtGui import QColor, QIcon
import numpy as np
import json
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from webbrowser import open as open_application
from itertools import islice
import pyqtgraph as pg
import pyqtgraph.exporters
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
        self.selected_files_gly = []
        self.selected_DQfiles = []
        self.selected_T1files = []
        self.selected_FFCfiles = []
        self.selected_DQMQfile = []
        self.dq_t2 = {}
        self.dq_comparison_linear = {}
        self.dq_comparison_distribution = {}
        self.tau_dictionary = {}
        self.ffc_dictionary = {}
        self.tab = None

        self.state()

        # Connect buttons to their respective slots
        self.ui.pushButton_DefaultFolder.clicked.connect(self.default_folder)
        self.ui.commandLinkButton.clicked.connect(self.open_url)
        self.check_for_updates()

        self.ui.tabWidget.currentChanged.connect(self.state)
    
        self.ui.btn_Save.clicked.connect(self.save_data)
        self.ui.btn_Save_2.clicked.connect(self.save_data)
        self.ui.btn_Save_6.clicked.connect(self.save_data)
        self.ui.btn_Save_3.clicked.connect(self.save_data)
        self.ui.btn_Save_7.clicked.connect(self.save_data)
        self.ui.btn_Load.clicked.connect(self.load_data)
        self.ui.btn_Load_2.clicked.connect(self.load_data)
        self.ui.btn_Load_3.clicked.connect(self.load_data)
        self.ui.btn_Load_4.clicked.connect(self.load_data)
        self.ui.btn_Load_5.clicked.connect(self.load_data)
        self.ui.btn_Phasing.clicked.connect(self.open_phasing_manual)
        
        self.ui.btn_SelectFiles.clicked.connect(self.open_select_dialog)
        self.ui.btn_Add.clicked.connect(self.add_select_dialog)
        self.ui.btn_SelectFiles_T1.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.btn_SelectFilesDQMQ.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.btn_SelectFiles_FFC.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.btn_SelectFilesDQ.clicked.connect(self.open_select_comparison_files_dialog)
        
        #self.ui.btn_SelectFiles.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_2.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_3.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_4.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_5.clicked.connect(self.clear_list)
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
        self.ui.btn_Plot1_2.clicked.connect(self.simulation)

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

        # Table Headers
        self.ui.table_SE.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=0: self.renameSection(self.ui.table_SE, index=0)
)
        self.ui.table_T1.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=2: self.renameSection(self.ui.table_T1, index=2)
)
        self.ui.table_DQ_2.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=2: self.renameSection(self.ui.table_DQ_2, index=1)
)
        # Connect table signals to slots
        self.ui.table_DQ.currentItemChanged.connect(self.update_dq_graphs)
        #self.ui.table_SE.currentItemChanged.connect(self.update_xaxis)
        self.ui.radioButton_4.clicked.connect(self.calculate_relaxation_time)
        self.ui.radioButton_5.clicked.connect(self.calculate_relaxation_time)
        self.ui.radioButton_6.clicked.connect(self.calculate_relaxation_time)
        self.ui.radioButton_16.clicked.connect(self.calculate_relaxation_time)
        self.ui.radioButton_17.clicked.connect(self.calculate_relaxation_time)
        self.ui.T1T2_fit_from.valueChanged.connect(self.calculate_relaxation_time)
        self.ui.radioButton_3.clicked.connect(self.hide_FFT_progress)
        self.ui.radioButton_2.clicked.connect(self.hide_FFT_progress)
        # Connect combobox signals to slots
        self.ui.comboBox.activated.connect(self.update_se_graphs)
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

        self.ui.checkBox_tau_1.clicked.connect(self.simulation)
        self.ui.checkBox_tau_2.clicked.connect(self.simulation)
        self.ui.checkBox_tau_3.clicked.connect(self.simulation)

        # Disable buttons initially
        self.disable_buttons()
        self.ui.progressBar.setHidden(True)
        self.ui.textEdit_5.setHidden(True)
        self.ui.textEdit_6.setHidden(True)

    def check_for_updates(self):
        current_version = '0.0.4'
        url = 'https://api.github.com/repos/timesama/M2-Approach/releases/latest'
        try:
            # Make a GET request to fetch the latest release data
            response = requests.get(url)
            response.raise_for_status()  # Raise an error for bad status codes
            latest_release = response.json()
            latest_version = latest_release['tag_name']

            # Compare the versions
            if latest_version != current_version:
                result = QMessageBox.information(self, "New Relaxyzer Available", f"A new version (Relaxyzer {latest_version}) is available.\nWould you like to update?", QMessageBox.Yes | QMessageBox.No)
                if result == QMessageBox.Yes:
                    self.open_url()
                
        except requests.RequestException as e:
            print(f"Failed to check for updates: {e}")

    def open_url(self):
        open_application('https://github.com/timesama/M2-Approach/releases')
    
    def update_file(self):
        i = self.ui.comboBox_4.currentIndex() + 1
        self.ui.btn_Phasing.setEnabled(True)

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
            self.update_se_graphs()

        #TODO sometime I should add the highlight of the certain point on graph, but I am too lazy
            
    def clear_list(self):
        if self.tab == 'SE' or self.tab == 'DQ':
            self.selected_files = []
            self.selected_files_gly = []
            self.selected_files = []  
            self.selected_files_gly = []
            self.ui.table_SE.setRowCount(0)
            self.ui.table_DQ.setRowCount(0)
            self.ui.SEWidget.clear()
            self.ui.DQ_Widget_1.clear()
            self.ui.DQ_Widget_2.clear()
            self.ui.DQ_Widget_3.clear()
            self.ui.DQ_Widget_4.clear()   
            self.ui.FFTWidget.clear()
            self.ui.FidWidget.clear()
        elif self.tab == 'DQ_Temp':
            self.selected_DQfiles = []
            self.dq_t2 = {}
            self.ui.table_DQ_2.setRowCount(0)
            self.ui.DQ_Widget_3.clear()
            self.ui.DQ_Widget_4.clear()
            self.ui.DQ_Widget_5.clear()
            self.ui.DQ_Widget_6.clear()
        elif self.tab == 'T1T2':
            self.selected_T1files = []
            self.tau_dictionary = {}
            self.ui.table_T1.setRowCount(0)
            self.ui.T1_Widget_1.clear()
            self.ui.T1_Widget_2.clear()
        elif self.tab == 'DQMQ':
            self.selected_DQMQfile = []
        elif self.tab == '23Model':
            self.selected_FFCfiles = []
            self.ffc_dictionary = {}
            self.ui.table_FFC_1.setRowCount(0)
            self.ui.btn_Plot1_2.setEnabled(False)
            self.ui.groupBox_7.setEnabled(False)
            self.ui.groupBox_6.setEnabled(False)
            self.ui.checkBox_3.setEnabled(False)
            self.ui.checkBox_3.setChecked(False)
            self.ui.table_FFC_1.clear()
        elif self.tab == 'Extra':
            pass


        if self.tab == 'T1T2':
            combobox = self.ui.comboBox_6
        elif self.tab == 'SE' or 'DQ':
            combobox = self.ui.comboBox_4
        elif self.tab == '23Model':
            combobox = self.ui.comboBox_8

        while combobox.count()>0:          
            combobox.removeItem(0)

    def delete_row(self):

        if self.tab == 'SE':
            table = self.ui.table_SE
            combobox = self.ui.comboBox_4
        elif self.tab == 'DQ':
            table = self.ui.table_DQ
            combobox = self.ui.comboBox_4
        elif self.tab =='T1T2':
            table = self.ui.table_T1
            combobox = self.ui.comboBox_6
        elif self.tab == '23Model':
            table = self.ui.table_FFC_1
            combobox = self.ui.comboBox_8
        else:
            return

        row = table.currentRow()
        table.removeRow(row)
        combobox.removeItem(row)
            
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

    def add_select_dialog(self):
        global State_multiple_files
        State_multiple_files = True
        dlg = OpenFilesDialog(self)
        if dlg.exec():
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

    def renameSection(self, table, index):
        current_header = table.horizontalHeaderItem(index).text()
        new_header, ok = QInputDialog.getText(self, "Rename Column", 
                                            f"Enter new name for the column '{current_header}':")
        if ok and new_header:
            table.horizontalHeaderItem(index).setText(new_header)
            self.update_xaxis(table, index)
            
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

        if self.ui.radioButton.isChecked():
            QMessageBox.information(self, "Data Saved", f"The figures have been saved to {os.path.dirname(file_path) + '/Result'}", QMessageBox.Ok)

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

        # Read name of filename
        filename = os.path.basename(file_path)

        # Read data
        try: 
            data = np.loadtxt(file_path)
            x, y, z = data[:, 0], data[:, 1], data[:, 2]
    
        except Exception as e:
            QMessageBox.warning(self, "Corrupt File", f"Couldn't read the {filename} beacuse {e} Only table data is available", QMessageBox.Ok)
            self.ui.FidWidget.clear()
            self.ui.FFTWidget.clear()
            self.ui.btn_Phasing.setEnabled(False)
            return
    
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

        if self.ui.radioButton_2.isChecked():
            number_of_points = 2**14
        else:
            number_of_points = 2**16

        Time_fid, Fid =  Cal.final_analysis_time_domain(Time, Re, Im, number_of_points)
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
            self.update_se_graphs()
        elif self.tab == 'DQ':
            self.update_dq_graphs()

    def extract_info(self, pattern):
        if pattern:
            info = pattern.group(1)
        else:
            info = '0'
        return info

    # Changed table name -> changed x axis #TODO add figure
    def update_xaxis(self, table, index):

        if self.tab == 'SE':
            figure = self.ui.SEWidget
        elif self.tab == 'T1T2':
            figure = self.ui.T1_Widget_2
        elif self.tab == 'DQ_Temp':
            figure = self.ui.DQ_Widget_6
        else:
            return

        name = table.horizontalHeaderItem(index).text()
        figure.getAxis('bottom').setLabel(name)

    # Working with graphs
    def update_graphs(self, x, y1, y2, y3, graph):
        graph.clear()
        graph.plot(x, y1, pen='k')
        graph.plot(x, y2, pen='r')
        graph.plot(x, y3, pen='b')

    # Working with SE graphs
    def update_se_graphs(self):
        
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
        self.ui.SEWidget.plot(x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=10)

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
        self.ui.DQ_Widget_1.plot(x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=10)

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
        graph_DQ_distr.plot(new_x, y, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=10)

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

        table = self.ui.table_DQMQ
        table.clear()

        try:
            Time, DQ, Ref = self.plot_original()
        except Exception as e:
            QMessageBox.warning(self, "Corrupt File", f"Couldn't read the file beacuse {e}", QMessageBox.Ok)
            self.clear_list()
            return
        num_rows  = len(Time)
        table.setRowCount(num_rows)

        for row in range(num_rows):
            table.setItem(row, 0, QTableWidgetItem(str(Time[row])))
            table.setItem(row, 1, QTableWidgetItem(str(DQ[row])))
            table.setItem(row, 2, QTableWidgetItem(str(Ref[row])))

        table.resizeColumnsToContents() # TODO ADD THIS EVERYWHERE!!!!
        table.setHorizontalHeaderLabels(['Time','DQ','Ref','nDQ']) # I don't know why i need this, but otherwise the columns get 1234 headres.
        self.ui.pushButton_DQMQ_1.setEnabled(True)
        self.ui.pushButton_DQMQ_2.setEnabled(False)
        self.ui.pushButton_DQMQ_3.setEnabled(False)
        self.ui.pushButton_DQMQ_4.setEnabled(True)
        
    def plot_original(self):
        file_path = self.selected_DQMQfile[0]

        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        legend = figure.addLegend()

        try:
            Time, DQ, Ref = Cal.read_data(file_path, 1)
        except:
            Time, DQ, Ref = Cal.read_data(file_path, 0)
    

        figure.plot(Time, DQ, pen='r', name = 'DQ')
        figure.plot(Time, Ref, pen='b', name = 'Ref')

        return Time, DQ, Ref

    def plot_norm(self):
        file_path = self.selected_DQMQfile[0]
        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        legend = figure.addLegend()
        noise_level = self.ui.noise.value()

        Time, DQ_norm, Ref_norm, _, _, _, _, _, _ = Cal.dqmq(file_path, 40, 100, 1, noise_level)

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
        noise_level = self.ui.noise.value()

        figure = self.ui.DQMQ_Widget
        figure.clear()
        legend = figure.addLegend()
        legend.clear()
        legend = figure.addLegend()

        Time, DQ_norm, Ref_norm, Diff, _, _, _, _, fitted_curve = Cal.dqmq(file_path, fit_from, fit_to, p, noise_level)

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
        noise_level = self.ui.noise.value()

        Time, _, _, _, DQ_normal, MQ_normal, Time0, nDQ, _ = Cal.dqmq(file_path, fit_from, fit_to, p, noise_level)

        figure.plot(Time, DQ_normal, pen='r', name = 'DQ')
        figure.plot(Time, MQ_normal, pen='b', name = 'Ref')
        figure.plot(Time0, nDQ, pen='k', symbol='o', symbolPen='k', symbolSize=10, name='nDQ')

        num_rows  = len(Time)
        for row in range(num_rows):
            table.setItem(row, 3, QTableWidgetItem(str(round(nDQ[row+1],4))))
        table.resizeColumnsToContents()
        
    # Relaxation time section
    def update_T12_table(self):
        def clean_line(line):
            while '\t\t' in line:
                line = line.replace('\t\t', '\t')
            return line.strip()
        
        selected_files = self.selected_T1files
        table = self.ui.table_T1
        combobox = self.ui.comboBox_6
        pattern = r'T1_(.*)\.dat|T2_(.*)\.dat'
        dictionary = self.tau_dictionary
        table.setRowCount(len(selected_files))

        x_axis = []

        for row, file in zip(range(table.rowCount()), selected_files):
            Folder = QTableWidgetItem(file)

            current_file = os.path.basename(file)

            try:
                x_axis = re.search(pattern,file).group(1)
            except:
                x_axis = row

            Filename = QTableWidgetItem(current_file)
            Temp = QTableWidgetItem(x_axis)
            table.setItem(row, 0, Folder)
            table.setItem(row, 1, Filename)
            table.setItem(row, 2, Temp)
        
            Time = []
            Signal = []

            try:
                # Read files as I create them from excel in spintrack
                with open(file, "r") as data:
                    #lines = [line.replace('\t\t\t', '').rstrip('\n') for line in data if not (line.rstrip('\n').endswith('\t\t\t\t'))]
                    lines = [clean_line(line.rstrip('\n')) for line in data if line.strip()]
                for line in lines[1:]:  # Skip the first line !!!
                    parts = line.split('\t')
                    time_value = float(parts[0])
                    signal_value = float(parts[1])
                    Time.append(time_value)
                    Signal.append(signal_value)

                combobox.addItem(f"{current_file}")
            except:
                try:
                    # read files regularly
                    data = np.loadtxt(file)
                    Time, Signal = data[:, 0], data[:, 1]

                    combobox.addItem(f"{current_file}")
                except FileNotFoundError as fnf_error:
                    QMessageBox.warning(self, "File Not Found", f"The file {file} was not found. Only table data is available.", QMessageBox.Ok)

                except Exception as e:
                    QMessageBox.warning(self, "Invalid Data", f"I couldn't read {current_file} because {e} Removing file from the table and file list.", QMessageBox.Ok)
                    for file_to_delete in selected_files:
                        if file_to_delete == file:
                            selected_files.remove(file)

            dictionary[file] = {"X Axis": [], "Time": [], "Signal": []}
            dictionary[file]["X Axis"].append(x_axis)
            dictionary[file]["Time"].extend(Time)
            dictionary[file]["Signal"].extend(Signal)
        
        self.ui.btn_Plot1.setEnabled(True)            
        combobox.setCurrentIndex(-1)

    def calculate_relaxation_time(self):
        table = self.ui.table_T1
        figure = self.ui.T1_Widget_1
        selected_file_idx = self.ui.comboBox_6.currentIndex()
        dictionary = self.tau_dictionary
        starting_point = int(self.ui.T1T2_fit_from.value())

        if self.ui.radioButton_16.isChecked():
        # milisec units
            denominator = 1
        else:
        # microsec units
            denominator = 1000

        if selected_file_idx == -1:
            return
        
        value_from_row = table.item(selected_file_idx, 0).text()
        Time_original = np.array(dictionary[value_from_row]['Time'])/denominator
        Signal_original = np.array(dictionary[value_from_row]['Signal'])

        Time = Time_original[starting_point:]
        Signal = Signal_original[starting_point:]

        Time_fit = np.arange(min(Time), max(Time) + 1, 1)

        try: 
            if self.ui.radioButton_4.isChecked():
                order = 1
            elif self.ui.radioButton_5.isChecked():
                order = 2
            else: 
                order = 3
        except:
            QMessageBox.warning(self, "No covariance", f"I am sorry, I couldn't fit with {order} exponents. Fitting with one.", QMessageBox.Ok)
            order = 1

        
        Time_fit, fitted_curve, tau1, tau2, tau3, R2, A1, A2, A3, decrease_order = Cal.fit_exponent(Time, Signal, order)
        
        self.ui.textEdit_error.setText(f"R² {R2}") 
        
        
        item1 = QTableWidgetItem(str(tau1))
        item2 = QTableWidgetItem(str(tau2))
        item3 = QTableWidgetItem(str(tau3))
        item11 = QTableWidgetItem(str(A1))
        item22 = QTableWidgetItem(str(A2))
        item33 = QTableWidgetItem(str(A3))

        table.setItem(selected_file_idx,3,item1)
        table.setItem(selected_file_idx,5,item2)
        table.setItem(selected_file_idx,7,item3)
        table.setItem(selected_file_idx,4,item11)
        table.setItem(selected_file_idx,6,item22)
        table.setItem(selected_file_idx,8,item33)
        
        
        figure.clear()
        figure.plot(Time_original, Signal_original, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)
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
                    

        graph.plot(x_axis, relaxation_time, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)

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

            except Exception as e: 
                QMessageBox.warning(self, "Invalid Data", f"The file {foldername} {e} Removed from the table and file list.", QMessageBox.Ok)
                table.removeRow(row)
                del self.selected_DQfiles[row]

                
            item = QTableWidgetItem(filename)

            pattern = r'DQ_table_(-?[0-9]+).csv'
            try:
                default_name = QTableWidgetItem(str(re.search(pattern, filename).group(1)))
            except:
                default_name = QTableWidgetItem(str(row+1))

            table.setItem(row, 0, item)
            table.setItem(row, 1, default_name)
        
        table.resizeColumnsToContents()

    def launch(self):
        try:
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
        except Exception as e:
            QMessageBox.warning(self, "Corrupted file", f"Couldn't analyse the {os.path.dirname(parent_folder)} because {e}", QMessageBox.Ok)

    def update_DQ_comparison_plot(self):
        cmap = pg.ColorMap([0, len(self.dq_t2)], [pg.mkColor('b'), pg.mkColor('r')])  # Blue to red
    
        self.dq_comparison_distribution = {'File name': [], 
                'X axis': [], 'Center': [], 'FWHM': [], 'Lorentz ratio': [], 
                'Fitting type': [], 'T2 limit': []}

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
            file_item = self.ui.table_DQ_2.item(row, 0)
            if file_name_item is not None:
                file_name = file_name_item.text()
                file = file_item.text()
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
            self.ui.DQ_Widget_3.plot(dq_time, t2, pen=None, symbolPen=None, symbol='o', symbolBrush=color, symbolSize=10, name=file_name)
            self.ui.DQ_Widget_4.plot(t2_lin, dq_norm, pen=None, symbolPen=None, symbol='o', symbolBrush=color, symbolSize=10, name=file_name)
            self.ui.DQ_Widget_4.plot(dq_fit, y_fit, pen=color)

            if coeff is not None:
                x_line = np.arange(0, 105.1, 0.1)
                y_line = np.polyval(coeff, x_line)
                
                self.ui.DQ_Widget_3.plot(x_line, y_line, pen=color) 

            self.dq_comparison_distribution['File name'].append(file)
            self.dq_comparison_distribution['X axis'].append(file_name)
            self.dq_comparison_distribution['Center'].append(cen)
            self.dq_comparison_distribution['FWHM'].append(fwhm)
            self.dq_comparison_distribution['Lorentz ratio'].append(w)
            self.dq_comparison_distribution['Fitting type'].append(text)
            self.dq_comparison_distribution['T2 limit'].append(time_max)


        self.ui.DQ_Widget_5.plot(comparison_par, center, pen='r', symbolPen=None, symbol='o', symbolBrush='r')
        self.ui.DQ_Widget_6.plot(comparison_par, weight, pen='b', symbolPen=None, symbol='o', symbolBrush='b')

    def write_collective_dictionary(self, dictionary, save_path_name):
        # file = name of the file 
        data = dictionary
        headers = data.keys()

        max_len = max(len(v) for v in data.values() if isinstance(v, list))
        rows = zip(*[v + v[-1:]*(max_len-len(v)) if len(v) < max_len else v for v in data.values()])

        # Write to CSV
        with open(save_path_name, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(headers)
            writer.writerows(rows)
        
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
        if self.tab == 'SE':
            table = self.ui.table_SE
            files = self.selected_files
            default_name = 'SE_'
        elif self.tab == 'DQ':
            table = self.ui.table_DQ
            files = self.selected_files
            default_name = 'Table_DQ_' + os.path.split(os.path.dirname(files[0]))[1]
        elif self.tab == 'DQ_Temp':
            table = self.ui.table_DQ_2
            files = self.selected_DQfiles
            path = os.path.dirname(files[0]) + '/Table_DQ_comparison_parametrs'
            self.write_collective_dictionary(self.dq_comparison_distribution, path)
            default_name = 'Table_DQ_comparison'

        elif self.tab == 'T1T2':
            table = self.ui.table_T1
            files = self.selected_T1files
            default_name = 'T'

        elif self.tab == 'DQMQ':
            table = self.ui.table_DQMQ
            files = self.selected_DQMQfile
            pattern = r'.*_(.*)'
            default_name = 'DQMQ_data_' + re.search(pattern, os.path.split(os.path.dirname(files[0]))[1] ).group(1)

        elif self.tab == '23Model':
            table = self.ui.table_FFC_1
            files = self.selected_FFCfiles
            default_name = 'FFC_'

        dialog = SaveFilesDialog(self)
        dialog.save_data_as_csv(self, table, files, default_name)

    def load_data(self):
        dlg = OpenFilesDialog(self)
        self.ui.btn_Save.setEnabled(False)
        if dlg.exec():
            tableName = dlg.selectedFiles()
            self.load_table_from_csv(tableName)

    def load_table_from_csv(self, tableName):

        self.clear_list()
        self.enable_buttons()

        file_path = tableName[0]
        try:
            files_list = file_path.strip().split('.')[0] + '_files.json'
            with open(files_list, 'r') as file:
                files = json.load(file)
        except Exception as e:  
            QMessageBox.warning(self, "File missing", f"Didn't find file list, only the tabular result is available", QMessageBox.Ok)
            files = None
            self.ui.comboBox_4.setEnabled(False)
            self.ui.btn_Phasing.setEnabled(False)
        
        if self.tab == 'SE':
            table = self.ui.table_SE
            self.selected_files = files
        elif self.tab == 'DQ':
            table = self.ui.table_DQ
            self.selected_files = files
        elif self.tab == 'DQ_Temp':
            table = self.ui.table_DQ_2
            self.selected_DQfiles = files
        elif self.tab == 'T1T2':
            table = self.ui.table_T1
            self.selected_T1files = files
        elif self.tab == 'DQMQ':
            table = self.ui.table_DQMQ
            self.selected_DQMQfile = files
        elif self.tab == '23Model':
            table = self.ui.table_FFC_1
            self.selected_FFCfiles = files

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
            self.update_se_graphs()
        elif self.tab == 'DQ':
            self.update_dq_graphs()
        elif self.tab == 'DQ_Temp':
            self.update_DQ_comparison()
        elif self.tab == 'T1T2':
            self.update_T12_table()
        elif self.tab == 'DQMQ':
            self.dq_mq_analysis()
            self.ui.pushButton_DQMQ_2.setEnabled(True)
            self.ui.pushButton_DQMQ_3.setEnabled(True)     
            self.ui.dq_min_3.setEnabled(True)
            self.ui.dq_max_3.setEnabled(True)
            self.ui.power.setEnabled(True)
        elif self.tab == '23Model':
            self.update_FFC_table()

    def save_figures(self, file_path, variable):
        # Set names
        parent_folder = os.path.dirname(file_path)
        os.makedirs(parent_folder + '/Result/', exist_ok=True)  

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

            dictionary[file] = {"File": [], "Freq": [], "Rate": [], "popt": None}
            dictionary[file]["File"].append(current_file)
            dictionary[file]["Freq"].extend(Omega)
            dictionary[file]["Rate"].extend(Rate)
            
            table.setRowCount(len(selected_files))
            self.ui.btn_Plot1_2.setEnabled(True)            
        
        for row, file in zip(range(table.rowCount()), selected_files):
            Folder = QTableWidgetItem(file)
            file_name = os.path.basename(file)

            Filename = QTableWidgetItem(file_name)

            table.setItem(row, 0, Folder)
            table.setItem(row, 1, Filename)

            combobox.addItem(f"{file_name}")
            combobox.setCurrentIndex(-1)

    def calculate_23_model(self):
        table = self.ui.table_FFC_1
        figure = self.ui.FFC_Widget_1
        selected_file_idx = self.ui.comboBox_8.currentIndex()
        dictionary = self.ffc_dictionary
        legend = figure.addLegend()  

        if legend is not None:
            legend.clear()

        if selected_file_idx == -1:
            print('does it ever happen?')
            return
        
        value_from_row = table.item(selected_file_idx, 0).text()
        Omega = np.array(dictionary[value_from_row]['Freq'], dtype=float)
        #Omega = Omeg * 10**6

        Rate = np.array(dictionary[value_from_row]['Rate'], dtype=float)
        Initial_coefficients = dictionary[value_from_row]['popt']

        if self.ui.checkBox_3.isChecked():
            fixed_CDD = self.ui.doubleSpinBox.value()
        else:
            fixed_CDD = None

        difference = 1

        try: 
            while difference > 1e-10:
                Omega_fit, fitted_curve, popt, R2 = Cal.fit_model(Omega, Rate, fixed_CDD, Initial_coefficients)
                difference = max(np.array(Initial_coefficients) - np.array(popt))
                Initial_coefficients = popt
                print('I use python as calculator, nahnahnah')
        except:
            Omega_fit, fitted_curve, popt, R2 = Cal.fit_model(Omega, Rate, fixed_CDD, Initial_coefficients)
            Initial_coefficients = popt


        Omega_fit, fitted_curve, popt, R2 = Cal.fit_model(Omega, Rate, fixed_CDD, Initial_coefficients)
        Initial_coefficients = popt    
        dictionary[value_from_row]['popt'] = popt

        CDD         = popt[0]
        tauc       = popt[1]
        A           = popt[2]
        C_trans     = popt[3]
        tau_trans   = popt[4]
        tau_res     = popt[5]

        tau_c_rate      = Cal.twod_model(Omega_fit, CDD, tauc, A, 0, tau_trans, tau_res)
        tau_trans_rate  = Cal.twod_model(Omega_fit, 0, 0, A, C_trans, tau_trans, tau_res)


        CDD_t         = QTableWidgetItem(str(popt[0]))
        tau_c_t       = QTableWidgetItem(str("{:.5e}".format(popt[1])))
        A_t           = QTableWidgetItem(str("{:.5e}".format(popt[2])))
        C_trans_t     = QTableWidgetItem(str("{:.5e}".format(popt[3])))
        tau_trans_t   = QTableWidgetItem(str("{:.5e}".format(popt[4])))
        tau_res_t     = QTableWidgetItem(str("{:.5e}".format(popt[5])))

        table.setItem(selected_file_idx, 3, CDD_t)
        table.setItem(selected_file_idx, 5, tau_c_t)
        table.setItem(selected_file_idx, 2, A_t)
        table.setItem(selected_file_idx, 4, C_trans_t)
        table.setItem(selected_file_idx, 6, tau_trans_t)
        table.setItem(selected_file_idx, 7, tau_res_t)

        for col in range(table.columnCount()):
            table.setColumnWidth(col, 120)

        figure.clear()
        figure.plot(Omega, Rate, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)
        figure.plot(Omega_fit, fitted_curve, pen='b', name = 'Original data')
        figure.plot(Omega_fit, tau_c_rate, pen = 'c',name='tau c')
        figure.plot(Omega_fit, tau_trans_rate, pen = 'm', name = 'tau trans')
        legend = figure.addLegend()  

        self.ui.checkBox_3.setEnabled(True)
        self.ui.btn_Plot1_2.setEnabled(True)
        self.ui.groupBox_7.setEnabled(True)
        self.ui.groupBox_6.setEnabled(True)

        self.ui.textEdit_error_2.setText(f"R² {R2}") 

    def simulation(self):
        table = self.ui.table_FFC_1
        figure = self.ui.FFC_Widget_2
        selected_file_idx = self.ui.comboBox_8.currentIndex()
        dictionary = self.ffc_dictionary

        value_from_row = table.item(selected_file_idx, 0).text()

        Omega = np.array(dictionary[value_from_row]['Freq'], dtype=float)
        Rate = np.array(dictionary[value_from_row]['Rate'], dtype=float)

        state = self.get_checkbox_state()

        if state == 'None':
            figure.clear()
            self.ui.textEdit_error_3.setText(f"No Data")     
        else:
            fitting, short, middle, long, popt, R2 = Cal.simulation(Omega, Rate, state)
            figure.clear()
            figure.plot(Omega, Rate, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)
            figure.plot(Omega, fitting, pen='b')
            if state == 'Two':
                figure.plot(Omega, short, pen='m')
                figure.plot(Omega, middle, pen='g')
                C1, C2, A, t1, t2= np.round(popt, decimals=5)
                C3=0
                t3 = 0
            elif state == 'Three':
                figure.plot(Omega, short, pen='m')
                figure.plot(Omega, middle, pen='g')
                figure.plot(Omega, long, pen='c')
                C1, C2, C3, A, t1, t2, t3 = np.round(popt, decimals=5)
            elif state == 'One':
                self.ui.textEdit_error_3.setText(f"R² {R2}")
                C1,A,t1 = np.round(popt, decimals=5)
                C2 = 0
                C3 = 0
                t2 = 0 
                t3 = 0
            else:
                return

            self.ui.textEdit_error_3.setText(f"R² {R2}\nA = {A}\nCDD = {C1}\tτ = {t1}\nCDD = {C2}\tτ = {t2}\nCDD = {C3}\tτ = {t3}") 

    def get_checkbox_state(self):
        checked_count = sum([self.ui.checkBox_tau_1.isChecked(),
                            self.ui.checkBox_tau_2.isChecked(),
                            self.ui.checkBox_tau_3.isChecked()])
        
        if checked_count == 0:
            state = 'None'
        elif checked_count == 1:
            state = 'One'
        elif checked_count == 2:
            state = 'Two'
        elif checked_count == 3:
            state = 'Three'
        else:
            state = 'Unknown'  # This case should technically never occur
        
        return state

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
    
class SaveFilesDialog(QFileDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        #self.setFileMode(QFileDialog.AnyFile)  # Allow selecting any file for saving
        self.setAcceptMode(QFileDialog.AcceptSave)  # Set the dialog to save mode

    def save_data_as_csv(self, directory, table, files, default_filename):
        # I have no fucjing idea, why the hell this function takes these arguments, but it doesn't work otherwise.
        initial_directory_file = "selected_folder.txt"
        try:
            with open(initial_directory_file, 'r') as file:
                directory = file.read().strip()
        except Exception as e:
            print(f"Couldn't read the initial directory: {e}")
            directory = os.path.dirname(sys.argv[0])


        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File As", directory + '/' + default_filename, "CSV files (*.csv)", options=options)

        if file_path:
            try:
                with open(file_path, 'w') as f:
                    for row in range(table.rowCount()):
                        row_values = []
                        for col in range(table.columnCount()):
                            item = table.item(row, col)
                            if item is not None:
                                row_values.append(item.text())
                            else:
                                row_values.append("")  # Handle empty cells
                        f.write(','.join(row_values) + '\n')

                files_list_path = os.path.splitext(file_path)[0] + '_files.json'

                with open(files_list_path, 'w') as file_list:
                    json.dump(files, file_list)

            except Exception as e:
                print(f"Failed to save file as CSV: {e}")

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
        global Frequency
        Integral = np.trapz(self.Real_freq_phased)
        left_mean = np.mean(self.Real_freq_phased[:100])
        right_mean = np.mean(self.Real_freq_phased[-100:])
        delta = left_mean - right_mean

        self.ui.Integral.setText(f"Integral: {round(Integral,3)}")
        self.ui.Delta.setText(f"Delta: {round(delta,7)}")

        Real_apod   = Cal.calculate_apodization(self.Real_freq_phased, Frequency)
        M2, T2 = Cal.calculate_M2(Real_apod, Frequency)

        self.ui.M2.setText(f"M₂: {round(M2,5)}")
        self.ui.T2.setText(f"T₂*: {round(T2,3)}")

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
