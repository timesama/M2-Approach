# This Python file uses the following encoding: utf-8
import sys, os, re, csv, requests, winreg
from PySide6.QtWidgets import QWidget,QTableWidgetItem, QTableWidget, QApplication, QMainWindow, QFileDialog, QTableWidgetItem, QInputDialog, QDialog, QMessageBox, QScrollArea
from PySide6.QtCore import QCoreApplication, Signal, QEvent, Qt
from PySide6.QtGui import QColor, QIcon, QKeySequence
import numpy as np
import json
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from webbrowser import open as open_application
from itertools import islice
import pyqtgraph as pg
from pyqtgraph import mkPen, ColorMap, mkColor, InfiniteLine
import pyqtgraph.exporters
from ui_Form import Ui_NMR
from ui_PhasingManual import Ui_Phasing as Ui_PhasingManual
import Calculator as Cal # Mathematical procedures
from controllers import (
    SETabController, DQTabController, DQTempTabController,
    T1T2TabController, DQMQTabController, GSTabController, ExtraTabController
)
from dialogs.open_files_dialog import OpenFilesDialog
import dialogs.open_files_dialog as open_files_dialog_module
from dialogs.notification_dialog import NotificationDialog
from dialogs.group_window import GroupWindow
from widgets.table_copy_enabler import TableCopyEnabler
from app_state import AppState

pg.CONFIG_OPTIONS['background'] = 'w'
pg.CONFIG_OPTIONS['foreground'] = 'k'

# IF you have ever wondered how bad code in Python looks like:
# here is the generous example of the masterpiece in bad coding.
# but it works and I don't care

# Global
Frequency = []
Re_spectra = []
Im_spectra = []

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

        self.app_state = AppState()

        self.selected_files = []
        self.selected_files_DQ_single = []
        self.app_state.dq_files = self.selected_files_DQ_single
        self.selected_files_gly = []
        self.selected_files_empty = []
        self.selected_DQfiles = []
        self.selected_T1files = []
        self.selected_DQMQfile = []
        self.app_state.dqmq_files = self.selected_DQMQfile
        self.selected_GSfiles = []
        self.window_array = np.array([])
        self.dq_t2 = {}
        self.dq_comparison_linear = {}
        self.dq_comparison_distribution = {}
        self.tau_dictionary = {}
        self.GS_dictionary = {}
        self.group_data_SE = {}
        self.group_data_T1T2 = {}
        self.group_data_SD = {}
        self.tab = None
        self.state_bad_code = False
        self.se_controller = SETabController(self)
        self.dq_controller = DQTabController(
            ui=self.ui,
            state=self.app_state,
            parent=self,
        )
        self.dq_temp_controller = DQTempTabController(self)
        self.t1t2_controller = T1T2TabController(ui=self.ui, state=self.app_state, parent=self)
        self.dqmq_controller = DQMQTabController(ui=self.ui, state=self.app_state, parent=self)
        self.gs_controller = GSTabController(self)
        self.extra_controller = ExtraTabController(self)

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
        self.ui.btn_Save_4.clicked.connect(self.save_data)

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
        self.ui.btn_SelectFilesDQ.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.btn_SelectFiles_GS.clicked.connect(self.open_select_comparison_files_dialog)

        self.ui.btn_ClearTable.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_2.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_3.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_4.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_5.clicked.connect(self.clear_list)

        self.ui.btn_DeleteRow.clicked.connect(self.delete_row)
        self.ui.btn_DeleteRow_1.clicked.connect(self.delete_row)
        self.ui.btn_DeleteRow_2.clicked.connect(self.delete_row)
        self.ui.btn_DeleteRow_3.clicked.connect(self.delete_row)

        self.ui.btn_Start.clicked.connect(self.analysis)

        self.ui.pushButton_DQMQ_1.clicked.connect(self.dqmq_controller.plot_original)
        self.ui.radioButton_Log.clicked.connect(self.dq_controller.plot_fit)
        self.ui.pushButton_DQMQ_4.clicked.connect(self.dqmq_controller.plot_norm)
        self.ui.pushButton_DQMQ_2.clicked.connect(self.dqmq_controller.plot_diff)
        self.ui.pushButton_DQMQ_3.clicked.connect(self.dqmq_controller.plot_nDQ)
        self.ui.btn_Plot1.clicked.connect(self.t1t2_controller.plot_relaxation_time)
        self.ui.btn_Plot_GS.clicked.connect(self.plot_sqrt_time)

        self.ui.pushButton_Eact.clicked.connect(self.plot_Arr)

        self.ui.pushButton_GroupSE.clicked.connect(self.open_group_window)
        self.ui.pushButton_GroupT1T2.clicked.connect(self.open_group_window)
        self.ui.pushButton_GroupSD.clicked.connect(self.open_group_window)

        # Graph setup
        self.setup_graph(self.ui.FFTWidget, "Frequency, MHz", "Amplitude, a.u", "FFT")
        self.setup_graph(self.ui.FidWidget, "Time, μs", "Amplitude", "NMR Signal")
        self.setup_graph(self.ui.SEWidget, "Temperature, °C", "Choose", "")
        self.setup_graph(self.ui.DQ_Widget_1, "DQ Filtering Time", "T₂*", "")
        self.setup_graph(self.ui.DQ_Widget_2, "X axis", "Norm. DQ Intensity", "")
        self.setup_graph(self.ui.DQ_Widget_4, "T₂*", "Norm. DQ Intensity", "FunctionFit")
        self.setup_graph(self.ui.DQ_Widget_5, "X axis", "Center", "")
        self.setup_graph(self.ui.T1_Widget_1, "Time, ms", "Signal", "")
        self.setup_graph(self.ui.T1_Widget_2, "X axis", "τ, ms", "")
        self.setup_graph(self.ui.DQMQ_Widget, "Time", "NMR signal", "")
        self.setup_graph(self.ui.DQMQ_Widget_DRes, "Dres/2pi, KHz", "P(Dres)", "")
        self.setup_graph(self.ui.GS_Widget_1, "√Time, √us", "Signal", "")
        self.setup_graph(self.ui.GS_Widget_2, "X axis", "√Time, √us", "")
        self.setup_graph(self.ui.DQ_Widget_polyFit, "T₂*", "Norm. DQ Intensity", "PolyFit")

        # Table Headers
        self.copy_enabler = TableCopyEnabler(self)

        self.ui.table_SE.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=0: self.renameSection(self.ui.table_SE, index=0)
)
        self.ui.table_T1.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=2: self.renameSection(self.ui.table_T1, index=2)
)
        self.ui.table_GS.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=2: self.renameSection(self.ui.table_GS, index=2)
)
        self.ui.table_DQ_2.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=2: self.renameSection(self.ui.table_DQ_2, index=1)
)
        # Connect table signals to slots
        self.ui.table_DQ.currentItemChanged.connect(self.dq_controller.update_graphs)
        self.ui.T1T2_FitWith1ExpButton.clicked.connect(self.t1t2_controller.change_exponential_order)
        self.ui.T1T2_FitWith2ExpButton.clicked.connect(self.t1t2_controller.change_exponential_order)
        self.ui.T1T2_FitWith3ExpButton.clicked.connect(self.t1t2_controller.change_exponential_order)

        self.ui.radioButton_16.clicked.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.radioButton_17.clicked.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.T1T2_fit_from.valueChanged.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.T1T2_fit_to.valueChanged.connect(self.t1t2_controller.calculate_relaxation_time)

        self.ui.checkBox_3.clicked.connect(self.calculate_sqrt_time)
        self.ui.radioButton_short.clicked.connect(self.calculate_sqrt_time)
        self.ui.radioButton_medium.clicked.connect(self.calculate_sqrt_time)
        self.ui.radioButton_long.clicked.connect(self.calculate_sqrt_time)
        self.ui.GS_fit_from_1.valueChanged.connect(self.calculate_sqrt_time)
        self.ui.GS_fit_to_1.valueChanged.connect(self.calculate_sqrt_time)

        self.ui.GS_beta.valueChanged.connect(self.calculate_sqrt_time)
        self.ui.GS_r2.valueChanged.connect(self.calculate_sqrt_time)
        self.ui.GS_m2.valueChanged.connect(self.calculate_sqrt_time)

        # Connect combobox signals to slots
        self.ui.comboBox_SE_chooseY.activated.connect(self.update_se_graphs)
        self.ui.comboBox_FunctionDQ.activated.connect(self.dq_controller.plot_fit)
        self.ui.T1T2_ChooseFileComboBox.activated.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.comboBox_7.activated.connect(self.calculate_sqrt_time)

        # Eact
        self.ui.radioButton_7.clicked.connect(self.plot_Arr)
        self.ui.radioButton_8.clicked.connect(self.plot_Arr)
        self.ui.checkBox_5.clicked.connect(self.plot_Arr)
        self.ui.Eact_start.valueChanged.connect(self.plot_Arr)
        self.ui.Eact_end.valueChanged.connect(self.plot_Arr)
        self.ui.pushButton_Done.clicked.connect(self.hide_Eact)


        # Connect change events
        self.ui.dq_min.valueChanged.connect(self.dq_controller.update_graphs)
        self.ui.dq_max.valueChanged.connect(self.dq_controller.update_graphs)
        self.ui.comboBox_4.activated.connect(self.update_file)
        self.ui.dq_min_3.valueChanged.connect(self.dqmq_controller.plot_diff)
        self.ui.dq_max_3.valueChanged.connect(self.dqmq_controller.plot_diff)
        self.ui.power.valueChanged.connect(self.dqmq_controller.plot_diff)

        # Disable buttons initially
        self.disable_buttons()
        self.ui.groupBox_EAct.setHidden(True)

        self.tab10_colors = [
            mkColor('#1f77b4'),  # blue
            mkColor('#ff7f0e'),  # orange
            mkColor('#2ca02c'),  # green
            mkColor('#d62728'),  # red
            mkColor('#9467bd'),  # purple
            mkColor('#8c564b'),  # brown
            mkColor('#e377c2'),  # pink
            mkColor('#7f7f7f'),  # gray
            mkColor('#bcbd22'),  # olive
            mkColor('#17becf')   # cyan
        ]

    def check_for_updates(self):
        current_version = '0.2.2'
        url = 'https://api.github.com/repos/timesama/M2-Approach/releases/latest'
        try:
            # Make a GET request to fetch the latest release data
            response = requests.get(url)
            response.raise_for_status()  # Raise an error for bad status codes
            latest_release = response.json()
            latest_version = latest_release['tag_name']

            # Compare the versions
            if latest_version != current_version:
                msg = QMessageBox(self)
                msg.setIcon(QMessageBox.Information)
                msg.setWindowTitle("New Relaxyzer Available")
                msg.setText(f"A new version (Relaxyzer {latest_version}) is available.\nWould you like to update?")
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
                msg.setWindowFlags(msg.windowFlags() | Qt.WindowStaysOnTopHint)
                result = msg.exec_()
                if result == QMessageBox.Yes:
                    self.open_url()

        except requests.RequestException as e:
            print(f"Failed to check for updates: {e}")

    def open_url(self):
        open_application('https://github.com/timesama/M2-Approach/releases')

    def update_file(self): # It detects the file ny its position in the table - should not be so. Take the filename from the dictionary
        self.ui.btn_Phasing.setEnabled(True)

        filename = self.ui.comboBox_4.currentText()
        if self.tab == 'SE':
            files = self.selected_files
        else:
            files = self.selected_files_DQ_single

        try:
            file_path = next((f for f in files if os.path.basename(f) == filename), None)
            i = next((i for i, f in enumerate(files) if os.path.basename(f) == filename), None)+1
        except:
            return

        if self.ui.checkBox_glycerol.isChecked() and self.ui.checkBox_baseline.isChecked():
            file_path_gly = self.selected_files_gly[i-1]
            file_path_empty = self.selected_files_empty[i-1]
        elif self.ui.checkBox_glycerol.isChecked():
            file_path_gly = self.selected_files_gly[i-1]
            file_path_empty=[]
        elif self.ui.checkBox_baseline.isChecked():
            file_path_gly = []
            file_path_empty = self.selected_files_empty[i-1]
        else:
            file_path_gly = []
            file_path_empty = []

        if self.ui.checkBox_Smooth.isChecked():
            window_first_value = self.ui.SmoothWindowFrom.value()
            window_last_value = self.ui.SmoothWindowTo.value()
            self.window_array = np.linspace(window_first_value, window_last_value, len(files), dtype=np.int32)
        else:
            self.window_array = np.array([])

        try:
            self.process_file_data(file_path, file_path_gly, file_path_empty, i)
        except:
            return

        # Update general figures
        if self.tab == 'DQ':
            self.highlight_row(self.ui.table_DQ, i)
            self.dq_controller.update_graphs()
        elif self.tab == 'SE':
            self.highlight_row(self.ui.table_SE, i)
            self.update_se_graphs()

        #TODO sometime I should add the highlight of the certain point on graph, but I am too lazy

    def clear_list(self):
        if self.tab == 'SE':
            self.app_state = AppState()
            self.selected_files = []
            self.selected_files_gly = []
            self.selected_files_empty = []
            self.ui.table_SE.setRowCount(0)
            self.ui.SEWidget.clear()
            self.ui.FFTWidget.clear()
            self.ui.FidWidget.clear()
            self.ui.btn_Start.setStyleSheet("background-color: none")
            self.group_data_SE = {}
        elif self.tab == 'DQ':
            self.selected_files_DQ_single = []
            self.app_state.dq_files = []
            self.selected_files_gly = []
            self.selected_files_empty = []
            self.ui.table_DQ.setRowCount(0)
            self.ui.DQ_Widget_1.clear()
            self.ui.DQ_Widget_2.clear()
            self.ui.DQ_Widget_4.clear()
            self.ui.FFTWidget.clear()
            self.ui.FidWidget.clear()
            self.ui.btn_Start.setStyleSheet("background-color: none")
        elif self.tab == 'DQ_Temp':
            self.selected_DQfiles = []
            self.dq_t2 = {}
            self.ui.table_DQ_2.setRowCount(0)
            self.ui.DQ_Widget_4.clear()
            self.ui.DQ_Widget_5.clear()
            self.ui.DQ_Widget_polyFit.clear()
        elif self.tab == 'T1T2':
            self.selected_T1files = []
            self.tau_dictionary = {}
            self.ui.table_T1.setRowCount(0)
            self.ui.T1_Widget_1.clear()
            self.ui.T1_Widget_2.clear()
            self.group_data_T1T2 = {}
        elif self.tab == 'DQMQ':
            self.selected_DQMQfile = []
            self.app_state.dqmq_files = []
        elif self.tab == 'GS':
            self.selected_GSfiles = []
            self.GS_dictionary = {}
            self.ui.table_GS.setRowCount(0)
            self.ui.GS_Widget_1.clear()
            self.ui.GS_Widget_2.clear()
            self.group_data_SD = {}
        elif self.tab == 'Extra':
            pass


        if self.tab == 'T1T2':
            combobox = self.ui.T1T2_ChooseFileComboBox
        elif self.tab in ('SE', 'DQ'):
            combobox = self.ui.comboBox_4
        elif self.tab == 'GS':
            combobox = self.ui.comboBox_7

        if self.tab != 'DQMQ' and  self.tab != 'DQ_Temp':
            while combobox.count()>0:
                combobox.removeItem(0)

        self.window_array = np.array([])

    def delete_row(self):

        if self.tab == 'SE':
            table = self.ui.table_SE
            combobox = self.ui.comboBox_4
        elif self.tab == 'DQ':
            table = self.ui.table_DQ
            combobox = self.ui.comboBox_4
        elif self.tab =='T1T2':
            table = self.ui.table_T1
            combobox = self.ui.T1T2_ChooseFileComboBox
            files = self.selected_T1files
        elif self.tab =='GS':
            table = self.ui.table_GS
            combobox = self.ui.comboBox_7
            files = self.selected_GSfiles
        else:
            return

        row = table.currentRow()
        if row == -1:
            QMessageBox.warning(self, "Cricket sounds", f"Select the row.", QMessageBox.Ok)
            return
        item = table.item(row,0).text()
        table.removeRow(row)
        combobox.removeItem(row)

        try:
            for file_to_delete in files:
                if file_to_delete == item:
                    files.remove(item)
        except:
            QMessageBox.warning(self, "Hidden", f"The row is hidden, but the file is not deleted.", QMessageBox.Ok)

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
        if self.tab == 'DQMQ':
            open_files_dialog_module.State_multiple_files = False
        else:
            open_files_dialog_module.State_multiple_files = True

        dlg = OpenFilesDialog(self)
        if dlg.exec():
            if self.tab == 'DQ_Temp':
                DQfileNames = dlg.selectedFiles()
                self.selected_DQfiles.extend(DQfileNames)
                self.update_DQ_comparison()
            elif self.tab == 'T1T2':
                while self.ui.T1T2_ChooseFileComboBox.count()>0:
                    self.ui.T1T2_ChooseFileComboBox.removeItem(0)
                T1fileNames = dlg.selectedFiles()
                self.selected_T1files.extend(T1fileNames)
                self.t1t2_controller.update_T12_table()
            elif self.tab == 'GS':
                while self.ui.comboBox_7.count()>0:
                    self.ui.comboBox_7.removeItem(0)
                GSfileNames = dlg.selectedFiles()
                self.selected_GSfiles.extend(GSfileNames)
                self.update_GS_table()
            elif self.tab == 'DQMQ':
                self.selected_DQMQfile = dlg.selectedFiles()
                self.dqmq_controller.dq_mq_analysis()

    def open_select_dialog(self):
        open_files_dialog_module.State_multiple_files = True
        dlg = OpenFilesDialog(self)

        if dlg.exec():
            self.clear_list()
            files = []
            fileNames = dlg.selectedFiles()
            files.extend(fileNames)
            self.ui.btn_Start.setEnabled(True)
            self.ui.btn_Start.setStyleSheet("background-color: green")
            self.ui.btn_Add.setEnabled(True)

        if self.tab == 'SE':
            self.selected_files = files
        else:
            self.selected_files_DQ_single = files
            self.app_state.dq_files = files

    def add_select_dialog(self):
        open_files_dialog_module.State_multiple_files = True
        dlg = OpenFilesDialog(self)

        if self.tab == 'SE':
            files = self.selected_files
        else:
            files = self.selected_files_DQ_single

        if dlg.exec():
            fileNames = dlg.selectedFiles()
            files.extend(fileNames)
            self.ui.btn_Start.setEnabled(True)
            self.ui.btn_Start.setStyleSheet("background-color: green")
            self.ui.btn_Add.setEnabled(True)

        if self.tab == 'SE':
            self.selected_files = files
        else:
            self.selected_files_DQ_single = files
            self.app_state.dq_files = files

    def open_select_dialog_glycerol(self):
        dlg = OpenFilesDialog(self)
        dlg.setWindowTitle("Select Glycerol Files")
        if dlg.exec():
            fileNames_gly = dlg.selectedFiles()
            self.selected_files_gly.extend(fileNames_gly)
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)

    def open_select_dialog_baseline(self):
        dlg = OpenFilesDialog(self)
        dlg.setWindowTitle("Select Baseline Files")
        if dlg.exec():
            fileNames_empty = dlg.selectedFiles()
            self.selected_files_empty.extend(fileNames_empty)
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)

    def open_phasing_manual(self):
        self.phasing_manual_window = PhasingManual()
        self.phasing_manual_window.read_data()
        self.phasing_manual_window.show()

        self.phasing_manual_window.closed.connect(self.after_phasing)

    def open_group_window(self):
        self.group_window = GroupWindow()

        if self.tab == 'SE':
            self.group_window.copy_table_data(self.ui.table_SE)
            self.group_data_SE = {}
        elif self.tab == 'T1T2':
            self.group_window.copy_table_data(self.ui.table_T1)
            self.group_data_T1T2 = {}
        elif self.tab == 'GS':
            self.group_window.copy_table_data(self.ui.table_GS)
            self.group_data_SD = {}

        if self.group_window.exec_() == QDialog.Accepted:
            data = self.group_window.group_dict
            # print("Group data received:", data)
        else:
            print("Group window was cancelled.")

        if self.tab == 'SE':
            self.group_data_SE = data
            self.update_se_graphs()

        elif self.tab == 'T1T2':
            self.group_data_T1T2 = data
            self.t1t2_controller.plot_relaxation_time()

        elif self.tab == 'GS':
            self.group_data_SD = data
            self.plot_sqrt_time()

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
            self.tab = 'GS'
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
        self.ui.comboBox_SE_chooseY.setEnabled(False)
        self.ui.comboBox_FunctionDQ.setEnabled(False)
        self.ui.radioButton_Log.setEnabled(False)
        self.ui.btn_Add.setEnabled(False)
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
        self.ui.comboBox_SE_chooseY.setEnabled(True)
        self.ui.comboBox_FunctionDQ.setEnabled(True)
        self.ui.comboBox_4.setEnabled(True)
        self.ui.btn_Add.setEnabled(True)

    def default_folder(self):
        folder_path = QFileDialog.getExistingDirectory(self, "Select Default Directory")

        if folder_path:
            try:
                key = winreg.CreateKey(winreg.HKEY_CURRENT_USER, r"Software\MyApp")
                winreg.SetValueEx(key, "SelectedFolder", 0, winreg.REG_SZ, folder_path)
                winreg.CloseKey(key)
                print(f"Default folder saved to registry: {folder_path}")
            except Exception as e:
                print(f"Failed to save the default folder to the registry: {e}")

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

        if self.tab == 'SE':
            files = self.selected_files
            self.ui.SEWidget.clear()
            self.ui.comboBox_SE_chooseY.setCurrentIndex(-1)
        else:
            files = self.selected_files_DQ_single
            self.ui.DQ_Widget_1.clear()
            self.ui.DQ_Widget_2.clear()
            self.ui.DQ_Widget_4.clear()
            self.ui.textEdit_4.setText("")
            self.ui.comboBox_FunctionDQ.setCurrentIndex(-1)

        if (len(files)==0):
            return

        if self.ui.checkBox_glycerol.isChecked():
            if self.selected_files_gly == []:
                self.open_select_dialog_glycerol()
                if len(self.selected_files_gly) < len(files):
                    QMessageBox.warning(self, "Invalid Data", f"The amount of Glycerol files is not the same as sample files. Adding glycerol files automatically.", QMessageBox.Ok)
                    last_file = self.selected_files_gly[-1]
                    num_to_add = len(files) - len(self.selected_files_gly)
                    self.selected_files_gly.extend([last_file] * num_to_add)
        else:
            self.selected_files_gly = []

        if self.ui.checkBox_baseline.isChecked():
            if self.selected_files_empty == []:
                self.open_select_dialog_baseline()
                if len(self.selected_files_empty) < len(files):
                    QMessageBox.warning(self, "Invalid Data", f"The amount of Empty files is not the same as sample files. Adding baseline files automatically.", QMessageBox.Ok)
                    last_file = self.selected_files_empty[-1]
                    num_to_add = len(files) - len(self.selected_files_empty)
                    self.selected_files_empty.extend([last_file] * num_to_add)
        else:
            self.selected_files_empty = []

        self.disable_buttons()
        self.ui.btn_SelectFiles.setEnabled(False)
        self.ui.btn_Load.setEnabled(False)
        self.ui.radioButton.setEnabled(False)
        self.ui.comboBox_4.setCurrentIndex(-1)

        if self.ui.checkBox_Smooth.isChecked():
            window_first_value = self.ui.SmoothWindowFrom.value()
            window_last_value = self.ui.SmoothWindowTo.value()
            self.window_array = np.linspace(window_first_value, window_last_value, len(files), dtype=np.int32)
        else:
            self.window_array = np.array([])

        for i, file_path in enumerate(files, start=1):

            if self.ui.checkBox_glycerol.isChecked():
                file_path_gly = self.selected_files_gly[i-1]
            else:
                file_path_gly = []
            if self.ui.checkBox_baseline.isChecked():
                file_path_empty = self.selected_files_empty[i-1]
            else:
                file_path_empty = []

            try:
                self.process_file_data(file_path, file_path_gly, file_path_empty, i)
            except:
                self.analysis_error(file_path, files)
                # return

        self.update_legends_and_dq_graphs()
        self.ui.btn_Start.setStyleSheet("background-color: none")

        if self.ui.radioButton.isChecked():
            QMessageBox.information(self, "Data Saved", f"The figures have been saved to {os.path.dirname(file_path) + '/Result'}", QMessageBox.Ok)

    def analysis_error(self, file_path, files):
        QMessageBox.warning(self, "Invalid Data", f"Couldn't read the {os.path.basename(file_path)}. Deleting the file.", QMessageBox.Ok)
        files.remove(file_path)
        self.ui.btn_SelectFiles.setEnabled(True)
        self.ui.btn_Start.setStyleSheet("background-color: none")

        self.ui.FidWidget.clear()
        self.ui.FFTWidget.clear()
        self.ui.btn_Phasing.setEnabled(False)
        self.enable_buttons()
        # self.analysis()

    def update_legends_and_dq_graphs(self):
        self.enable_buttons()

        if self.tab == 'DQ':
            self.dq_controller.update_graphs()

    def process_file_data(self, file_path, file_path_gly, file_path_empty, i):
        global Frequency, Re_spectra, Im_spectra

        # Read name of filename
        filename = os.path.basename(file_path)

        if self.ui.checkBox_baseline.isChecked():
            subtract = True
        else:
            subtract = False

        # Read data
        data = np.loadtxt(file_path)
        x, y, z = data[:, 0], data[:, 1], data[:, 2]

        Time, Re, Im = Cal.analysis_time_domain(file_path, file_path_empty, subtract)

        if self.ui.checkBox_glycerol.isChecked():
            # Glycerol
            Time_r, Re_r, Im_r = Cal.analysis_time_domain(file_path_gly, [], False)
            Time, Re, Im = Cal.magnet_inhomogenity_correction(Time, Time_r, Re, Re_r, Im, Im_r)

        if self.ui.checkBox_long_component.isChecked():
            # Longcomponent
            Time, Re, Im = Cal.subtract_long_component(Time, Re, Im)

        Amp = Cal._calculate_amplitude(Re, Im)
        self.update_graphs(Time, Amp, Re, Im, self.ui.FidWidget)

        number_of_points = 2**16

        Time_fid, Fid =  Cal.final_analysis_time_domain(Time, Re, Im, number_of_points)
        Frequency = Cal._calculate_frequency_scale(Time_fid)

        if len(Frequency) < len(Fid):
            diff = len(Fid) - len(Frequency)
            df = Frequency[1] - Frequency[0]
            extra_freq = Frequency[-1] + np.arange(1, diff+1) * df
            Frequency = np.concatenate([Frequency, extra_freq])
        elif len(Frequency) > len(Fid):
            Frequency = Frequency[:len(Fid)]


        FFT = np.fft.fftshift(np.fft.fft(Fid))

        if not self.window_array.size == 0:
            try:
                window = self.window_array[i-1]
                RealPart = savgol_filter(np.real(FFT), window, polyorder=1)
                ImaginaryPart = savgol_filter(np.imag(FFT), window, polyorder=1)
                FFT = np.array(RealPart + 1j * ImaginaryPart)
            except Exception as e:
                print(f'Couldnt smooth because {e}')

        # 8. Simple baseline
        Amp_spectra, Re_spectra, Im_spectra = Cal._simple_baseline_correction(FFT)
        # 9. Cal.apodization
        Real_apod = Cal._calculate_apodization(Re_spectra, Frequency)


        # Update FFT graphs
        self.update_graphs(Frequency, Amp_spectra, Re_spectra, Im_spectra, self.ui.FFTWidget)
        if self.ui.comboBox_4.findText(filename) == -1:
            self.ui.comboBox_4.addItem(filename)

        if self.ui.comboBox_4.currentIndex() == -1:
            M2, T2 = Cal._calculate_M2(Real_apod, Frequency)

            if self.tab == 'SE':
                match = re.search(r'.*_(-?\s*\d+\.?\d*).*.dat', filename)
                temperature = self.extract_info(match)

                short_from = self.ui.SC_short.value() *2
                short_to = self.ui.SC_short_2.value() *2

                long_from = self.ui.SC_long.value() *2
                long_to = self.ui.SC_long_2.value() *2

                absolute = self.ui.radioButton_absolute.isChecked()
                times = [int(short_from), int(short_to), int(long_from), int(long_to)]

                SFC = Cal.calculate_SC(Amp, times, absolute)
                self.ui.table_SE.setRowCount(i) # HERE to set the right amount of rows
                self.fill_table(self.ui.table_SE, temperature, SFC, M2, T2, i)

                self.ui.table_SE.setItem(i-1, 4, QTableWidgetItem(filename))

                if self.ui.radioButton.isChecked():
                    self.save_figures(file_path, filename)

            elif self.tab == 'DQ':
                match = re.search(r'_(\d+\.\d+)_', filename)
                dq_time = self.extract_info(match)

                Amplitude = Cal._calculate_amplitude(y, z)
                DQ = Cal.calculate_DQ_intensity(x, Amplitude)
                self.ui.table_DQ.setRowCount(i)
                self.fill_table(self.ui.table_DQ, dq_time, DQ, M2, T2, i)

                if self.ui.radioButton.isChecked():
                    self.save_figures(file_path, dq_time)
        else:
            pass

    def after_phasing(self):
        global Frequency, Re_spectra, Im_spectra

        i = self.ui.comboBox_4.currentIndex()

        Real_apod   = Cal._calculate_apodization(Re_spectra, Frequency) #(math procedure)
        Amp_spectra = Cal._calculate_amplitude(Re_spectra, Im_spectra)

        # Update FFT graph
        self.update_graphs(Frequency, Amp_spectra, Re_spectra, Im_spectra, self.ui.FFTWidget)

        M2, T2 = Cal._calculate_M2(Real_apod, Frequency)
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
            self.dq_controller.update_graphs()

    def extract_info(self, pattern):
        if pattern:
            info = pattern.group(1)
        else:
            info = '0'
        return info

    def update_xaxis(self, table, index):

        if self.tab == 'SE':
            figure = self.ui.SEWidget
        elif self.tab == 'T1T2':
            figure = self.ui.T1_Widget_2
        elif self.tab == 'GS':
            figure = self.ui.GS_Widget_2
        elif self.tab == 'DQ_Temp':
            figure = self.ui.DQ_Widget_5
        else:
            return

        name = table.horizontalHeaderItem(index).text()
        figure.getAxis('bottom').setLabel(name)

    # Working with graphs
    def update_graphs(self, x, y1, y2, y3, graph):
        graph.clear()
        graph.plot(x, y1, pen=mkPen('k', width=3))
        graph.plot(x, y2, pen=mkPen('r', width=2))
        graph.plot(x, y3, pen=mkPen('b', width=2))

    # Working with SE graphs
    def update_se_graphs(self):
        self.se_controller.update_graphs()


    # Working with tables
    def update_DQ_comparison(self):
        table = self.ui.table_DQ_2
        table.setRowCount(len(self.selected_DQfiles))

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
            item_name = QTableWidgetItem()

            pattern = r'Table_DQ_(-?[0-9]+).*.csv'
            try:
                item_xaxis = float(re.search(pattern, filename).group(1))
                item_name.setData(Qt.EditRole, item_xaxis)
            except:
                item_name.setData(Qt.EditRole, float(row+1))

            table.setItem(row, 0, item)
            table.setItem(row, 1, item_name)

        table.resizeColumnsToContents()

        self.launch()

    def launch(self):
        try:
            self.dq_t2 = {}
            for row, parent_folder in enumerate(self.selected_DQfiles, start=0):
                # Read data from file
                data = np.loadtxt(parent_folder, delimiter=',')

                # Read the DQ filtering time, DQ amlitude and corresponding T2*
                dq_t2 = data[:, [4, 5]]
                self.dq_t2[row] = dq_t2
            self.update_DQ_comparison_plot()
        except Exception as e:
            QMessageBox.warning(self, "Corrupted file", f"Couldn't analyse the {os.path.dirname(parent_folder)} because {e}", QMessageBox.Ok)

    def update_DQ_comparison_plot(self):
        cmap = pg.ColorMap([0, len(self.dq_t2)], [pg.mkColor('b'), pg.mkColor('r')])  # Blue to red

        self.dq_comparison_distribution = {'File name': [], 
                'X axis': [], 'Center': [], 'FWHM': [], 'Lorentz ratio': [], 
                'Fitting type': [], 'T2 limit': []}

        legend = self.ui.DQ_Widget_4.addLegend(offset=(0,0))
        legend_functionsCenters = self.ui.DQ_Widget_5.addLegend()
        self.ui.DQ_Widget_5.clear()
        self.ui.DQ_Widget_4.clear()
        self.ui.DQ_Widget_polyFit.clear()

        if legend is not None:
            legend.clear()
            legend.setPen((0, 0, 0)) 

        if legend_functionsCenters is not None:
            legend_functionsCenters.clear()
            legend_functionsCenters.setPen((0, 0, 0)) 


        legend.anchor(itemPos=(1, 0), parentPos=(1, 0))

        # This is a mess :(
        center_g = []
        center_l = []
        center_v = []
        center_d = []

        bold_g = []
        bold_l = []
        bold_v = []

        comparison_par = []

        for row, (key, data) in zip(range(self.ui.table_DQ_2.rowCount()), self.dq_t2.items()):
            file_name_item = self.ui.table_DQ_2.item(row, 1)
            file_item = self.ui.table_DQ_2.item(row, 0)
            if file_name_item is not None:
                file_name = file_name_item.text()
                file = file_item.text()
                if file_name != 'hide': #TODO OHMYGOD SERIOUSLY?????????????
                    try:
                        _comp_par = float(file_name)
                    except:
                        _comp_par = 0
                        self.ui.table_DQ_2.setItem(row, 1, QTableWidgetItem('0'))
                    comparison_par.append(_comp_par)
                else:
                    continue

            # Initial arrays
            t2_lin = data[:,0] #T2*
            dq_norm = data[:,1]

            # Initial fitting parameters with statistical functions
            p = [1, 5, 5, 0]
            b=([0, 0, 0, 0, 0], [ np.inf, np.inf, np.inf, 1, np.inf]) # Voigt
            b1=([0, 0, 0, -10], [np.inf, np.inf, np.inf, np.inf]) # Gauss and Lorenz

            # create an array of T2* for plotting fitted function
            t2_fit = np.arange(0, np.max(t2_lin) + 0.001, 0.01)

            # Fit with 3 different functions and get centers
            # if text == 'Gauss':
            params, _ = curve_fit(Cal.gaussian, t2_lin, dq_norm, p0=p, bounds=b1)
            cen_g = params[1]
            center_g.append(cen_g)
            fwhm_g = params[2]
            bold_g.append(fwhm_g)

            # elif text == 'Lorenz':
            params, _ = curve_fit(Cal.lorenz, t2_lin, dq_norm, p0=p, bounds=b1)
            cen_l = params[1]
            center_l.append(cen_l)
            fwhm_l = params[2]
            bold_l.append(fwhm_l)

            # elif text == 'Pseudo Voigt':
            params, _ = curve_fit(Cal.voigt, t2_lin, dq_norm,  bounds = b)
            cen_v = params[1]
            center_v.append(cen_v)
            fwhm_v = params[2]
            bold_v.append(fwhm_v)

            y_fit = Cal.voigt(t2_fit, *params)

            # Fitting with polynomial and derivative search:
            t2_detivative_plot, nDQ_derivative_plot, center_derivative = Cal.derivative_peak_find(t2_lin, dq_norm)
            center_d.append(center_derivative)

            # Draw a graph T2* versus nDQ data and fitted functions
            color = tuple(cmap.map(key))
            self.ui.DQ_Widget_4.plot(t2_lin, dq_norm, pen=None, symbolPen=None, symbol='o', symbolBrush=color, symbolSize=10, name=file_name)
            self.ui.DQ_Widget_4.plot(t2_fit, y_fit, pen=color)

            self.ui.DQ_Widget_polyFit.plot(t2_lin, dq_norm, pen=None, symbolPen=None, symbol='o', symbolBrush=color, symbolSize=10)
            self.ui.DQ_Widget_polyFit.plot(t2_detivative_plot, nDQ_derivative_plot, pen=color)

            # Add centers values to the table:
            def set_numeric_item(table, row, col, value):
                item = QTableWidgetItem()
                item.setData(Qt.EditRole, round(float(value), 2))
                table.setItem(row, col, item)

            # Add centers values to the table
            set_numeric_item(self.ui.table_DQ_2, row, 2, cen_g)
            set_numeric_item(self.ui.table_DQ_2, row, 3, cen_l)
            set_numeric_item(self.ui.table_DQ_2, row, 4, cen_v)
            set_numeric_item(self.ui.table_DQ_2, row, 5, center_derivative)

            # Add FWHM values to the table
            set_numeric_item(self.ui.table_DQ_2, row, 6, fwhm_g)
            set_numeric_item(self.ui.table_DQ_2, row, 7, fwhm_l)
            set_numeric_item(self.ui.table_DQ_2, row, 8, fwhm_v)

            #TODO: make disctionary self.dq_comparison_distribution['File name'].append(file) updates

        self.ui.DQ_Widget_5.plot(comparison_par, center_g, pen='r', symbolPen=None, symbol='o', symbolBrush='r', name='Gaus')
        self.ui.DQ_Widget_5.plot(comparison_par, center_l, pen='b', symbolPen=None, symbol='o', symbolBrush='b', name='Lorenz')
        self.ui.DQ_Widget_5.plot(comparison_par, center_v, pen='k', symbolPen=None, symbol='o', symbolBrush='k', name='Voigt')
        self.ui.DQ_Widget_5.plot(comparison_par, center_d, pen='g', symbolPen=None, symbol='o', symbolBrush='g', name='Derivative')

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
            files = self.selected_files_DQ_single
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

        elif self.tab == 'GS':
            table = self.ui.table_GS
            files = self.selected_GSfiles
            default_name = 'SpinDiffusion_'

        elif self.tab == 'DQMQ':
            table = self.ui.table_DQMQ
            files = self.selected_DQMQfile
            pattern = r'.*_(.*)'
            try:
                default_name = 'DQMQ_data_' + re.search(pattern, os.path.split(os.path.dirname(files[0]))[1] ).group(1)
            except:
                default_name = 'DQMQ_data_'
        dialog = SaveFilesDialog(self)
        dialog.save_data_as_csv(self, table, files, default_name)

        if self.state_bad_code == True:
            self.bad_code_makes_more_bad_code()

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
            self.selected_files_DQ_single = files
            self.app_state.dq_files = files
        elif self.tab == 'DQ_Temp':
            table = self.ui.table_DQ_2
            self.selected_DQfiles = files
        elif self.tab == 'T1T2':
            table = self.ui.table_T1
            self.selected_T1files = files
        elif self.tab == 'DQMQ':
            table = self.ui.table_DQMQ
            self.selected_DQMQfile = files
            self.app_state.dqmq_files = files
        elif self.tab == 'GS':
            table = self.ui.table_GS
            self.selected_GSfiles = files

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
            self.dq_controller.update_graphs()
        elif self.tab == 'DQ_Temp':
            self.update_DQ_comparison()
        elif self.tab == 'T1T2':
            self.t1t2_controller.update_T12_table()
        elif self.tab == 'DQMQ':
            self.ui.pushButton_DQMQ_1.setEnabled(True)
            self.ui.pushButton_DQMQ_2.setEnabled(True)
            self.ui.pushButton_DQMQ_3.setEnabled(True)
            self.ui.pushButton_DQMQ_4.setEnabled(True)
            self.ui.dq_min_3.setEnabled(True)
            self.ui.dq_max_3.setEnabled(True)
            self.ui.power.setEnabled(True)
            self.dqmq_controller.plot_nDQ_on_Load()

        elif self.tab == 'GS':
            self.update_GS_table()

    def save_figures(self, file_path, variable):
        parent_folder = os.path.dirname(file_path)
        result_folder = os.path.join(parent_folder, 'Result')
        os.makedirs(result_folder, exist_ok=True)

        graph_fft = self.ui.FFTWidget
        graph_fid = self.ui.FidWidget

        fft_file_path = os.path.join(result_folder, f"FFT_{variable}.png")
        fid_file_path = os.path.join(result_folder, f"NMR_{variable}.png")
        fft_csv_path  = os.path.join(result_folder, f"FFT_{variable}.csv")
        fid_csv_path  = os.path.join(result_folder, f"NMR_{variable}.csv")

        pg.QtGui.QGuiApplication.processEvents()

        # Export PNGs
        exporter_fft = pg.exporters.ImageExporter(graph_fft.plotItem)
        exporter_fft.parameters()['width'] = 1000
        exporter_fft.export(fft_file_path)

        exporter_fid = pg.exporters.ImageExporter(graph_fid.plotItem)
        exporter_fid.parameters()['width'] = 1000
        exporter_fid.export(fid_file_path)

        # Save FFT FID CSV
        curves = graph_fft.plotItem.listDataItems()

        _, y2 = curves[1].getData()
        x, y3 = curves[2].getData()

        with open(fft_csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for xi, y2i, y3i in zip(x, y2, y3):
                writer.writerow([xi, y2i, y3i])

        curves = graph_fid.plotItem.listDataItems()

        _, y2 = curves[1].getData()
        x, y3 = curves[2].getData()

        with open(fid_csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for xi, y2i, y3i in zip(x, y2, y3):
                writer.writerow([xi, y2i, y3i])

    def load_data_and_check_validity(self, file_path):
        try:
            data = np.loadtxt(file_path)
            filename = os.path.basename(file_path)
            # Check that there are 3 columns
            if data.shape[1] != 3:
                QMessageBox.warning(self, "Invalid Data Format", f"I can't read the {filename} file, it should have 3 columns exactly, deleting it.", QMessageBox.Ok)
                self.ui.btn_SelectFiles.setEnabled(True)
                self.ui.radioButton.setEnabled(True)

                if self.tab =='SE':
                    files = self.selected_files
                else:
                    files = self.selected_files_DQ_single
                files.clear()
                return False
            return True

        except:
            QMessageBox.warning(self, "Invalid Data", f"I can't read the {filename} file, deleting it.", QMessageBox.Ok)
            self.selected_files.remove(file_path)
            self.ui.btn_Start.setStyleSheet("background-color: none")

            if self.selected_files == []:
                self.disable_buttons()
                return False
            else:
                return False

    # FFC analysis
    def bad_code_makes_more_bad_code(self):
        dictionary = self.tau_dictionary
        dialog = SaveFilesDialog(self)
        basename = os.path.basename(self.selected_T1files[0])

        save = not all(data.get('T1 1', 0) == 0 for data in dictionary.values())
        dialog.save_file_in_sef(self, dictionary, 'T1 1', 1, basename, save)

        save = not all(data.get('T1 2', 0) ==0 for data in dictionary.values())
        dialog.save_file_in_sef(self, dictionary, 'T1 2', 2, basename, save)

        save = not all(data.get('T1 3', 0) == 0 for data in dictionary.values())
        dialog.save_file_in_sef(self, dictionary, 'T1 3', 3, basename, save)

    # GS Goldman Shen
    def update_GS_table(self):
        def clean_line(line):
            while '\t\t' in line:
                line = line.replace('\t\t', '\t')
            return line.strip()

        def create_dictionary(dictionary, file, addition, x_axis, sqrtTime, short, medium, long):
            dictionary[file + addition]["X Axis"].append(x_axis)
            dictionary[file + addition]["sqrtTime"].extend(sqrtTime)
            dictionary[file + addition]["short"].extend(short)
            dictionary[file + addition]["medium"].extend(medium)
            dictionary[file + addition]["long"].extend(long)
            return dictionary

        selected_files = self.selected_GSfiles
        table = self.ui.table_GS
        combobox = self.ui.comboBox_7

        pattern = r'_([0-9]+)\.dat$'
        dictionary = self.GS_dictionary

        try:
            table.setRowCount(len(selected_files))
            x_axis = []
            for row, file in zip(range(table.rowCount()), selected_files):
                dictionary[file] = {"X Axis": [], "sqrtTime": [], "short": [], "medium": [], "long": []}

                sqrtTime = []
                short = []
                medium = []
                long = []
                Folder = QTableWidgetItem(file)

                current_file = os.path.basename(file)

                try:
                    x_axis = re.search(pattern,current_file).group(1)
                except:
                    x_axis = row

                try:
                    with open(file, "r") as data:
                        lines = [clean_line(line.rstrip('\n')) for line in data if line.strip()]
                    for line in lines[1:]:  # Skip the first line !!!
                        parts = line.split('\t')
                        time_value = float(parts[0])
                        signal_value1 = float(parts[4])
                        signal_value2 = float(parts[5])
                        signal_value3 = float(parts[6])

                        sqrtTime.append(time_value)
                        short.append(signal_value1)
                        medium.append(signal_value2)
                        long.append(signal_value3)

                    combobox.addItem(f"{current_file}")
                except Exception as e:
                    QMessageBox.warning(self, "Invalid Data", f"I couldn't read {current_file} because {e} Removing file from the table and file list.", QMessageBox.Ok)
                    for file_to_delete in selected_files:
                        if file_to_delete == file:
                            selected_files.remove(file)
                    return

                Filename = QTableWidgetItem(current_file)
                XValue = QTableWidgetItem(x_axis)
                table.setItem(row, 0, Folder)
                table.setItem(row, 1, Filename)
                table.setItem(row, 2, XValue)
                dictionary[file]["X Axis"].append(x_axis)
                dictionary[file]["sqrtTime"].extend(sqrtTime)
                dictionary[file]["short"].extend(short)
                dictionary[file]["medium"].append(medium)
                dictionary[file]["long"].extend(long)

        except Exception as e:
            QMessageBox.warning(self, "Error", f"Something {e} went wrong. Try again.", QMessageBox.Ok)
            self.clear_list()

    def calculate_sqrt_time(self):
        selected_file_idx = self.ui.comboBox_7.currentIndex()
        if selected_file_idx == -1:
            return
        table = self.ui.table_GS
        figure = self.ui.GS_Widget_1

        dictionary = self.GS_dictionary

        value_from_row = table.item(selected_file_idx, 0).text()

        if self.ui.checkBox_3.isChecked():
        # transform to sqrt
            Time_original = np.sqrt(np.array(dictionary[value_from_row]['sqrtTime']).flatten())
        else: # time is already in sqrt
            Time_original = np.array(dictionary[value_from_row]['sqrtTime']).flatten()

        short_original = np.array(dictionary[value_from_row]['short']).flatten()
        medium_original = np.array(dictionary[value_from_row]['medium']).flatten()
        long_original = np.array(dictionary[value_from_row]['long']).flatten()

        self.ui.GS_fit_from_1.setMinimum(Time_original[0])
        self.ui.GS_fit_from_1.setMaximum(Time_original[-1])

        self.ui.GS_fit_to_1.setMinimum(Time_original[15])
        self.ui.GS_fit_to_1.setMaximum(Time_original[-1])

        from_val = self.ui.GS_fit_from_1.value()
        to_val = self.ui.GS_fit_to_1.value()

        # Find the closest indices to the selected Time range
        starting_point = (np.abs(Time_original - from_val)).argmin()
        ending_point = (np.abs(Time_original - to_val)).argmin() + 1  # +1 to include the endpoint

        Time = Time_original[starting_point:ending_point]
        short = short_original[starting_point:ending_point]
        medium = medium_original[starting_point:ending_point]
        long = long_original[starting_point:ending_point]


        if self.ui.radioButton_short.isChecked():
            Signal = short
            Signal_original = short_original
        elif self.ui.radioButton_medium.isChecked():
            Signal = medium
            Signal_original = medium_original
        elif self.ui.radioButton_long.isChecked():
            Signal = long
            Signal_original = long_original

        beta    =   self.ui.GS_beta.value()
        r2    =   self.ui.GS_r2.value()
        M2    =   self.ui.GS_m2.value()
        try:
            Time_fit, fitted_curve, sqrtT, R2 = Cal.linear_fit_GS(Time, Signal)
            d = Cal.calculate_domain_size(sqrtT, beta, r2, M2)

            self.ui.textEdit_error_2.setText(f"R² {R2}")

            item = QTableWidgetItem(str(sqrtT))

            table.setItem(selected_file_idx,3,item)

            item2 = QTableWidgetItem(str(d))

            table.setItem(selected_file_idx,4,item2)

            figure.clear()
            figure.plot(Time_original, Signal_original, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)
            figure.plot(Time_fit, fitted_curve, pen='b')

            dictionary[value_from_row]['sqrtT'] = sqrtT
            dictionary[value_from_row]['d'] = d

            self.ui.btn_Plot1.setEnabled(True)
        except Exception as e:
            figure.clear()
            QMessageBox.warning(self, "Error", f"Something {e} went wrong. Try again.", QMessageBox.Ok)

    def plot_sqrt_time(self):

        table = self.ui.table_GS
        graph = self.ui.GS_Widget_2
        column = 3 #time

        graph.clear()
        if table.rowCount() < 1:
            return

        x_axis = []
        sqrtT =[]
        number = 1
        for row in range(table.rowCount()):
            try:
                C = float(table.item(row, 2).text())
                x_axis.append(C)
            except:
                x_axis.append(number)

            try:
                T = float(table.item(row, column).text())
                sqrtT.append(T)
            except:
                sqrtT.append(0)
            number = number + 1

        graph.plot(x_axis, sqrtT, pen=None, symbolPen=None, symbol='o', symbolBrush='r', symbolSize=10)

        if self.group_data_SD:
            for i, (group_number, group_rows) in enumerate(self.group_data_SD.items()):
                group_x = []
                group_y = []

                for row_data in group_rows:
                    try:
                        # Adjust these indexes as per your group data structure
                        x_val = float(row_data[2])  # e.g., x value
                        y_val = float(row_data[column])
                        group_x.append(x_val)
                        group_y.append(y_val)
                    except (ValueError, IndexError):
                        continue

                sorted_points = sorted(zip(group_x, group_y), key=lambda p: p[0])
                if len(sorted_points) > 1:
                    xs, ys = zip(*sorted_points)
                    color = self.tab10_colors[i % len(self.tab10_colors)]

                    graph.plot(
                        xs, ys,
                        pen={'color': color, 'width': 2},
                        symbol='o',
                        symbolBrush=color,
                        symbolPen=None,
                        symbolSize=8
                    )

    # Eact
    def plot_Arr(self):
        self.ui.groupBox_EAct.setHidden(False)
        table = self.ui.table_SE
        graph = self.ui.SEWidget

        Temperature = np.array(self.read_column_values(table, 0))
        T2 = np.array(self.read_column_values(table, 3))

        sorted_indices = np.argsort(Temperature)
        Temperature = Temperature[sorted_indices]
        T2 = T2[sorted_indices]

        if self.ui.checkBox_5.isChecked():
            Temperature = Temperature + 273.15 #Transfer to K

        starting_point = int(self.ui.Eact_start.value())
        ending_point = -(int(self.ui.Eact_end.value()))
        if ending_point == 0:
            ending_point = None
        try:
            Temperature = Temperature[starting_point:ending_point]
            T2 = T2[starting_point:ending_point]
            reciprocal_temperature, lnT2 = Cal.calculate_Arrhenius_ax(Temperature, T2)
            # fit linear
            Temp_fit, fitted_curve, Eact, R2 = Cal.calculate_Eact(reciprocal_temperature, lnT2, self.ui.radioButton_8.isChecked())

            graph.clear()
            graph.plot(reciprocal_temperature, lnT2, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=10)
            graph.plot(Temp_fit, fitted_curve, pen='b')
            self.setup_graph(graph, "1000/T, 𝐾⁻¹", "ln(τ)", "")

            self.ui.textEdit_EAct.setText(f"Eact = {Eact}\nR² {R2}")

        except Exception as e:
            QMessageBox.warning(self, "Error", f"Something {e} went wrong. Try again.", QMessageBox.Ok)

    def hide_Eact(self):
        self.ui.groupBox_EAct.setHidden(True)
        self.ui.SEWidget.clear()
        self.setup_graph(self.ui.SEWidget, "", "", "")
        self.update_xaxis(self.ui.table_SE, 0)

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
        self.setAcceptMode(QFileDialog.AcceptSave)  # Set the dialog to save mode

    def save_data_as_csv(self, directory, table, files, default_filename):
        # I have no fucjing idea, why the hell this function takes these arguments, but it doesn't work otherwise.
        try:
            key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, r"Software\MyApp", 0, winreg.KEY_READ)
            directory, _ = winreg.QueryValueEx(key, "SelectedFolder")
            winreg.CloseKey(key)
        except FileNotFoundError:
            return os.path.dirname(sys.argv[0])
        except Exception as e:
            print(f"Couldn't read the initial directory: {e}")
            return os.path.dirname(sys.argv[0])


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

    def save_file_in_sef(self, wtf, dictionary, tau, n, begin, save):
        if save == False:
            return

        try:
            key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, r"Software\MyApp", 0, winreg.KEY_READ)
            directory, _ = winreg.QueryValueEx(key, "SelectedFolder")
            winreg.CloseKey(key)
        except FileNotFoundError:
            return os.path.dirname(sys.argv[0])
        except Exception as e:
            print(f"Couldn't read the initial directory: {e}")
            return os.path.dirname(sys.argv[0])

        try:
            options = QFileDialog.Options()
            default_add = '_recalculated_tau_' + str(n)
            file_path, _ = QFileDialog.getSaveFileName(self, "Save File As", directory + '/'+ begin + default_add, "SEF files (*.sef)", options=options)
            with open(file_path, 'w') as f:
                f.write('STELAR Export File\n')  
                f.write('\n')  
                f.write('_BRLX______\t_T1________\t_R1________\t____%err___\t+-err(R1)__\tZone\tFile\n')  
                f.write('\n')
                for key in dictionary:
                    Frequency = dictionary[key]['X Axis']
                    T1 = dictionary[key][tau]
                    if T1 == 0:
                        Omega_1 = 0
                    else:
                        Omega_1 = 1000/T1
                    f.write(f'{Frequency}\t{T1}\t{Omega_1}\n')
        except:
            QMessageBox.warning(self, "Save failed", f"Sorry, couldn't save the separate file in .sef for magnetization.", QMessageBox.Ok)
            return
        #TODO: try to find reasons for that maybe? 
        # When the load happends, there is no key T1 1 created in dictionary
        # So ppb should read value from the table

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
        # print(self.Smooth)

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

        Real_apod   = Cal._calculate_apodization(self.Real_freq_phased, Frequency)
        M2, T2 = Cal._calculate_M2(Real_apod, Frequency)

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
