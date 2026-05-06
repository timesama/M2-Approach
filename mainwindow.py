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
import Calculator as Cal # Mathematical procedures
from controllers import (
    SETabController, DQTabController, DQTempTabController,
    T1T2TabController, DQMQTabController, GSTabController, GeneralSEDQController
)
from dialogs.open_files_dialog import OpenFilesDialog
import dialogs.open_files_dialog as open_files_dialog_module
from dialogs.notification_dialog import NotificationDialog
from dialogs.group_window import GroupWindow
from dialogs.save_files_dialog import SaveFilesDialog
from dialogs.phasing_manual import PhasingManual
from widgets.table_copy_enabler import TableCopyEnabler
from app_state import AppState
from controllers.table_columns import T1Columns, GSColumns, DQTempColumns
from utils.ui_busy import busy_cursor
import logging

logger = logging.getLogger(__name__)

pg.CONFIG_OPTIONS['background'] = 'w'
pg.CONFIG_OPTIONS['foreground'] = 'k'

# IF you have ever wondered how bad code in Python looks like:
# here is the generous example of the masterpiece in bad coding.
# but it works and I don't care


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
        self.phased_spectra_SE = {}
        self.phased_spectra_DQ = {}
        self.tab = None
        self.se_controller = SETabController(ui=self.ui, state=self.app_state, parent=self)
        self.dq_controller = DQTabController(
            ui=self.ui,
            state=self.app_state,
            parent=self,
        )
        self.dq_temp_controller = DQTempTabController(ui=self.ui, state=self.app_state, parent=self)
        self.t1t2_controller = T1T2TabController(ui=self.ui, state=self.app_state, parent=self)
        self.dqmq_controller = DQMQTabController(ui=self.ui, state=self.app_state, parent=self)
        self.gs_controller = GSTabController(ui=self.ui, state=self.app_state, parent=self)
        self.general_se_dq_controller = GeneralSEDQController(ui=self.ui, state=self.app_state, parent=self, se_controller=self.se_controller, dq_controller=self.dq_controller)

        self.state()

        # Connect buttons to their respective slots
        self.ui.pushButton_DefaultFolder.clicked.connect(self.default_folder)
        self.ui.commandLinkButton.clicked.connect(self.open_url)
        self.check_for_updates()

        self.ui.tabWidget.currentChanged.connect(self.state)

        self.ui.btn_Save.clicked.connect(self.save_data)
        self.ui.DQMQ_Button_Save.clicked.connect(self.save_data)
        self.ui.btn_Save_6.clicked.connect(self.save_data)
        self.ui.btn_Save_3.clicked.connect(self.save_data)
        self.ui.btn_Save_4.clicked.connect(self.save_data)

        self.ui.btn_Load.clicked.connect(self.load_data)
        self.ui.btn_Load_2.clicked.connect(self.load_data)
        self.ui.btn_Load_3.clicked.connect(self.load_data)
        self.ui.DQMQ_Button_Load.clicked.connect(self.load_data)
        self.ui.btn_Load_5.clicked.connect(self.load_data)
        self.ui.btn_Phasing.clicked.connect(self.general_se_dq_controller.open_phasing_manual)

        self.ui.btn_SelectFiles.clicked.connect(self.open_select_dialog)
        self.ui.btn_Add.clicked.connect(self.add_select_dialog)
        self.ui.btn_SelectFiles_T1.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.DQMQ_Button_SelectFiles.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.btn_SelectFilesDQ.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.btn_SelectFiles_GS.clicked.connect(self.open_select_comparison_files_dialog)

        self.ui.btn_ClearTable.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_2.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_3.clicked.connect(self.clear_list)
        self.ui.SE_Button_ClearTable.clicked.connect(self.clear_list)
        self.ui.btn_ClearTable_5.clicked.connect(self.clear_list)

        self.ui.btn_DeleteRow.clicked.connect(self.delete_row)
        self.ui.SE_Button_DeleteRow.clicked.connect(self.delete_row)
        self.ui.btn_DeleteRow_2.clicked.connect(self.delete_row)
        self.ui.btn_DeleteRow_3.clicked.connect(self.delete_row)

        self.ui.btn_Start.clicked.connect(self.general_se_dq_controller.analysis)

        self.ui.DQMQ_Button_PlotOriginal.clicked.connect(
            self.dqmq_controller.plot_original
        )
        self.ui.radioButton_Log.clicked.connect(self.dq_controller.plot_fit)
        self.ui.DQMQ_Button_PlotNorm.clicked.connect(self.dqmq_controller.plot_norm)
        self.ui.DQMQ_Button_CalculateIntegralSum.clicked.connect(
            self.dqmq_controller.calculate_integral_sum
        )
        self.ui.DQMQ_Button_CalculateDres.clicked.connect(
            self.dqmq_controller.calculate_dres
        )
        self.ui.DQMQ_DoubleSpinBox_IntegralShift.editingFinished.connect(
            self.dqmq_controller.update_integral_sum_shift
        )
        self._connect_dqmq_workflow_signals()
        self.ui.btn_Plot1.clicked.connect(self.t1t2_controller.plot_relaxation_time)
        self.ui.btn_Plot_GS.clicked.connect(self.gs_controller.plot_sqrt_time)

        self.ui.pushButton_GroupT1T2.clicked.connect(self.open_group_window)
        self.ui.pushButton_GroupSD.clicked.connect(self.open_group_window)

        # Graph setup
        self.setup_graph(self.ui.FFTWidget, "Frequency, MHz", "Amplitude, a.u", "FFT")
        self.setup_graph(self.ui.FidWidget, "Time, μs", "Amplitude", "NMR Signal")
        self.setup_graph(self.ui.SE_PlotWidget_Main, "Temperature, °C", "Choose", "")
        self.setup_graph(self.ui.DQ_Widget_1, "DQ Filtering Time", "T₂*", "")
        self.setup_graph(self.ui.DQ_Widget_2, "X axis", "Norm. DQ Intensity", "")
        self.setup_graph(self.ui.DQ_Widget_4, "T₂*", "Norm. DQ Intensity", "FunctionFit")
        self.setup_graph(self.ui.DQ_Widget_5, "X axis", "Center", "")
        self.setup_graph(self.ui.T1_Widget_1, "Time, ms", "Signal", "")
        self.setup_graph(self.ui.T1_Widget_2, "X axis", "τ, ms", "")
        self.setup_graph(self.ui.DQMQ_PlotWidget_Signal, "Time", "NMR signal", "")
        self.setup_graph(self.ui.DQMQ_PlotWidget_Dres, "Dres/2π, KHz", "P(Dres)", "")
        self.setup_graph(self.ui.GS_Widget_1, "√Time, √us", "Signal", "")
        self.setup_graph(self.ui.GS_Widget_2, "X axis", "√Time, √us", "")
        self.setup_graph(self.ui.DQ_Widget_polyFit, "T₂*", "Norm. DQ Intensity", "PolyFit")
        self.se_controller.connect_signals()

        # Table Headers
        self.copy_enabler = TableCopyEnabler(self)

        self.ui.table_T1.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=T1Columns.X_AXIS: self.renameSection(self.ui.table_T1, index=T1Columns.X_AXIS)
)
        self.ui.table_GS.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=GSColumns.X_AXIS: self.renameSection(self.ui.table_GS, index=GSColumns.X_AXIS)
)
        self.ui.table_DQ_2.horizontalHeader().sectionDoubleClicked.connect(
    lambda index=DQTempColumns.NAME: self.renameSection(self.ui.table_DQ_2, index=DQTempColumns.NAME)
)
        self._apply_table_header_order()
        # Connect table signals to slots
        self.ui.table_DQ.currentItemChanged.connect(self.dq_controller.update_graphs)
        self.ui.table_DQ.itemSelectionChanged.connect(self.dq_controller.update_graphs)
        self.ui.table_T1.itemSelectionChanged.connect(self.t1t2_controller.plot_relaxation_time)
        self.ui.T1T2_FitWith1ExpButton.clicked.connect(self.t1t2_controller.change_exponential_order)
        self.ui.T1T2_FitWith2ExpButton.clicked.connect(self.t1t2_controller.change_exponential_order)
        self.ui.T1T2_FitWith3ExpButton.clicked.connect(self.t1t2_controller.change_exponential_order)

        self.ui.radioButton_16.clicked.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.radioButton_17.clicked.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.T1T2_fit_from.valueChanged.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.T1T2_fit_to.valueChanged.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.DSB_ExpFitting1.valueChanged.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.DSB_ExpFitting2.valueChanged.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.DSB_ExpFitting3.valueChanged.connect(self.t1t2_controller.calculate_relaxation_time)

        self.ui.checkBox_3.clicked.connect(self.gs_controller.calculate_sqrt_time)
        self.ui.radioButton_short.clicked.connect(self.gs_controller.calculate_sqrt_time)
        self.ui.radioButton_medium.clicked.connect(self.gs_controller.calculate_sqrt_time)
        self.ui.radioButton_long.clicked.connect(self.gs_controller.calculate_sqrt_time)
        self.ui.GS_fit_from_1.valueChanged.connect(self.gs_controller.calculate_sqrt_time)
        self.ui.GS_fit_to_1.valueChanged.connect(self.gs_controller.calculate_sqrt_time)

        self.ui.GS_beta.valueChanged.connect(self.gs_controller.calculate_sqrt_time)
        self.ui.GS_r2.valueChanged.connect(self.gs_controller.calculate_sqrt_time)
        self.ui.GS_m2.valueChanged.connect(self.gs_controller.calculate_sqrt_time)

        # Connect combobox signals to slots
        self.ui.comboBox_FunctionDQ.activated.connect(self.dq_controller.plot_fit)
        self.ui.T1T2_ChooseFileComboBox.activated.connect(self.t1t2_controller.calculate_relaxation_time)
        self.ui.comboBox_7.activated.connect(self.gs_controller.calculate_sqrt_time)


        # Connect change events
        self.ui.dq_min.valueChanged.connect(self.dq_controller.update_graphs)
        self.ui.dq_max.valueChanged.connect(self.dq_controller.update_graphs)
        self.ui.comboBox_4.activated.connect(self.update_file)

        # Disable buttons initially
        self.disable_buttons()
        self.ui.SE_GroupBox_EAct.setHidden(True)

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
            self.general_se_dq_controller.process_file_data(file_path, file_path_gly, file_path_empty, i)
        except Exception:
            return

        # Update general figures
        if self.tab == 'DQ':
            self.highlight_row(self.ui.table_DQ, i)
            self.dq_controller.update_graphs()
        elif self.tab == 'SE':
            self.highlight_row(self.ui.SE_Table_Data, i)
            self.update_se_graphs()

        #TODO sometime I should add the highlight of the certain point on graph, but I am too lazy

    def clear_list(self):
        if self.tab == 'SE':
            self.app_state.dq_files = []
            self.app_state.dqmq_files = []
            self.selected_files = []
            self.selected_files_gly = []
            self.selected_files_empty = []
            self.ui.SE_Table_Data.setRowCount(0)
            self.ui.SE_PlotWidget_Main.clear()
            self.ui.FFTWidget.clear()
            self.ui.FidWidget.clear()
            self.ui.btn_Start.setStyleSheet("background-color: none")
            self.group_data_SE = {}
            self.phased_spectra_SE = {}
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
            self.phased_spectra_DQ = {}
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
            self.dqmq_controller.reset_cached_results()
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
            table = self.ui.SE_Table_Data
            combobox = self.ui.comboBox_4
            files = self.selected_files
        elif self.tab == 'DQ':
            table = self.ui.table_DQ
            combobox = self.ui.comboBox_4
            files = self.selected_files_DQ_single
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

        table.removeRow(row)
        if combobox.count() > row:
            combobox.removeItem(row)

        try:
            if files and len(files) > row:
                files.pop(row)
            if self.tab in ('SE', 'DQ') and self.selected_files_gly and len(self.selected_files_gly) > row:
                self.selected_files_gly.pop(row)
            if self.tab in ('SE', 'DQ') and self.selected_files_empty and len(self.selected_files_empty) > row:
                self.selected_files_empty.pop(row)
        except Exception:
            QMessageBox.warning(self, "Delete warning", "Row was removed from table, but file lists may be out of sync.", QMessageBox.Ok)

        if self.tab == 'SE':
            self.se_controller.update_graphs()
        elif self.tab == 'DQ':
            self.dq_controller.update_graphs()

        if self.tab in ('SE', 'DQ') and combobox.count() > 0:
            combobox.setCurrentIndex(min(row, combobox.count() - 1))
            self.update_file()
        elif self.tab in ('SE', 'DQ'):
            self.ui.FidWidget.clear()
            self.ui.FFTWidget.clear()

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

        # self.ui.SE_Table_Data.selectRow(5)
        # self.ui.SE_Table_Data.currentRow()

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
                self.dq_temp_controller.update_DQ_comparison()
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
                self.gs_controller.update_GS_table()
            elif self.tab == 'DQMQ':
                self.selected_DQMQfile = dlg.selectedFiles()
                self.app_state.dqmq_files = self.selected_DQMQfile
                self.dqmq_controller.state.dqmq_files = self.selected_DQMQfile
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

    def open_group_window(self):
        self.group_window = GroupWindow()

        if self.tab == 'SE':
            self.group_window.copy_table_data(self.ui.SE_Table_Data)
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
            self.gs_controller.plot_sqrt_time()


    def _connect_dqmq_workflow_signals(self):
        analysis_parameter_widgets = [
            "DQMQ_DoubleSpinBox_FitFrom",
            "DQMQ_DoubleSpinBox_FitTo",
            "DQMQ_DoubleSpinBox_Power",
            "DQMQ_DoubleSpinBox_Noise",
            "DQMQ_DoubleSpinBox_TimeShift",
            "DQMQ_DoubleSpinBox_SmoothFrom",
            "DQMQ_DoubleSpinBox_SmoothTo",
            "DQMQ_DoubleSpinBox_SmoothWindow",
        ]
        for widget_name in analysis_parameter_widgets:
            self._connect_dqmq_signal(
                widget_name,
                "editingFinished",
                self.dqmq_controller.on_analysis_parameter_editing_finished,
            )

        visibility_checkboxes = [
            "DQMQ_CheckBox_DQ",
            "DQMQ_CheckBox_Reference",
            "DQMQ_CheckBox_MQ",
            "DQMQ_CheckBox_NDQ",
            "DQMQ_CheckBox_Difference",
            "DQMQ_CheckBox_TailFitting",
            "DQMQ_CheckBox_IntegralSum",
            "DQMQ_CheckBox_DresFitting",
        ]
        for widget_name in visibility_checkboxes:
            self._connect_dqmq_signal(
                widget_name,
                "toggled",
                self.dqmq_controller.on_visibility_checkbox_changed,
            )

        dres_parameter_widgets = [
            "DQMQ_DoubleSpinBox_DresK",
            "DQMQ_DoubleSpinBox_DresCenter1",
            "DQMQ_DoubleSpinBox_DresWidth1",
            "DQMQ_DoubleSpinBox_DresCenter2",
            "DQMQ_DoubleSpinBox_DresWidth2",
            "DQMQ_DoubleSpinBox_DresFraction1",
            "DQMQ_DoubleSpinBox_DresWeibullBeta",
        ]
        for widget_name in dres_parameter_widgets:
            self._connect_dqmq_signal(
                widget_name,
                "editingFinished",
                self.dqmq_controller.mark_dres_stale,
            )

        self._connect_dqmq_signal(
            "DQMQ_ComboBox_Kernel",
            "currentIndexChanged",
            self.dqmq_controller.mark_dres_stale,
        )
        self._connect_dqmq_signal(
            "DQMQ_RadioButton_OneDres",
            "toggled",
            self.dqmq_controller.mark_dres_stale,
        )
        self._connect_dqmq_signal(
            "DQMQ_RadioButton_TwoDres",
            "toggled",
            self.dqmq_controller.mark_dres_stale,
        )

    def _connect_dqmq_signal(self, widget_name, signal_name, slot):
        widget = getattr(self.ui, widget_name, None)
        if widget is None:
            logger.warning(
                "DQMQ widget %s is missing; signal not connected",
                widget_name,
            )
            return

        signal = getattr(widget, signal_name, None)
        if signal is None:
            logger.warning(
                "DQMQ widget %s has no %s signal; signal not connected",
                widget_name,
                signal_name,
            )
            return

        signal.connect(lambda *args: slot())

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
        self.ui.comboBox_FunctionDQ.setEnabled(False)
        self.ui.radioButton_Log.setEnabled(False)
        self.ui.btn_Add.setEnabled(False)
        self.ui.btn_Plot1.setEnabled(False)
        self.ui.DQMQ_Button_PlotOriginal.setEnabled(False)
        self.ui.DQMQ_Button_PlotNorm.setEnabled(False)

    def enable_buttons(self):
        self.ui.btn_SelectFiles.setEnabled(True)
        self.ui.btn_Start.setEnabled(True)
        self.ui.btn_Save.setEnabled(True)
        self.ui.radioButton.setEnabled(True)
        self.ui.btn_Load.setEnabled(True)
        self.ui.dq_min.setEnabled(True)
        self.ui.dq_max.setEnabled(True)
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

    def extract_info(self, pattern):
        if pattern:
            info = pattern.group(1)
        else:
            info = '0'
        return info

    def update_xaxis(self, table, index):

        if self.tab == 'SE':
            figure = self.ui.SE_PlotWidget_Main
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

    def _apply_table_header_order(self):
        self.ui.table_T1.setHorizontalHeaderLabels([
            "X axis", "tau 1", "A 1", "tau 2", "A 2", "tau 3", "A 3", "File name", "Folder"
        ])
        self.ui.table_GS.setHorizontalHeaderLabels([
            "X axis", "sqrt time", "d, nm", "File name", "Folder"
        ])
        self.ui.table_DQ_2.setHorizontalHeaderLabels([
            "Name", "Center Gauss", "Center Lorenz", "Center Voigt", "Center y",
            "FWHM Gauss", "FWHM Lorenz", "FWHM Voigt", "Folder"
        ])
        self.ui.table_T1.resizeColumnsToContents()
        self.ui.table_GS.resizeColumnsToContents()
        self.ui.table_DQ_2.resizeColumnsToContents()

    # Working with graphs
    def update_graphs(self, x, y1, y2, y3, graph):
        graph.clear()
        graph.plot(x, y1, pen=mkPen('k', width=3))
        graph.plot(x, y2, pen=mkPen('r', width=2))
        graph.plot(x, y3, pen=mkPen('b', width=2))

    # Working with SE graphs
    def update_se_graphs(self):
        self.se_controller.update_graphs()


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
            table = self.ui.SE_Table_Data
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
            table = self.ui.DQMQ_Table_Data
            files = self.selected_DQMQfile
            pattern = r'.*_(.*)'
            try:
                selected_folder = os.path.split(os.path.dirname(files[0]))[1]
                match = re.search(pattern, selected_folder)
                default_name = 'DQMQ_data_' + match.group(1)
            except (IndexError, AttributeError, TypeError):
                default_name = 'DQMQ_data_'
        dialog = SaveFilesDialog(self)
        dialog.save_data_as_csv(self, table, files, default_name)
        if self.tab == 'DQMQ' and dialog.last_saved_file_path:
            self.dqmq_controller.save_integral_sum_result(dialog.last_saved_file_path)
            self.dqmq_controller.save_dres_result(dialog.last_saved_file_path)
        if self.tab in ('SE', 'DQ') and dialog.last_saved_file_path:
            files_json = os.path.splitext(dialog.last_saved_file_path)[0] + '_files.json'
            phased_data = self.phased_spectra_SE if self.tab == 'SE' else self.phased_spectra_DQ
            with open(files_json, 'w') as f:
                json.dump({"files": files, "phased": phased_data}, f)


    def load_data(self):
        dlg = OpenFilesDialog(self)
        self.ui.btn_Save.setEnabled(False)
        if dlg.exec():
            tableName = dlg.selectedFiles()
            self.load_table_from_csv(tableName)

    def load_table_from_csv(self, tableName):
        with busy_cursor():
            self._load_table_from_csv_impl(tableName)

    def _load_table_from_csv_impl(self, tableName):

        self.clear_list()
        self.enable_buttons()
        self.ui.comboBox_4.clear()

        file_path = tableName[0]
        logger.info("Loading saved table: %s", file_path)
        try:
            files_list = file_path.strip().split('.')[0] + '_files.json'
            with open(files_list, 'r') as file:
                files_payload = json.load(file)
                if isinstance(files_payload, dict):
                    files = files_payload.get("files", [])
                    phased = files_payload.get("phased", {})
                else:
                    files = files_payload
                    phased = {}
        except Exception as e:
            QMessageBox.warning(self, "File missing", f"Didn't find file list, only the tabular result is available", QMessageBox.Ok)
            files = None
            phased = {}
            self.ui.comboBox_4.setEnabled(False)
            self.ui.btn_Phasing.setEnabled(False)

        try:
            phased_path = os.path.splitext(file_path)[0] + '_phased.json'
            with open(phased_path, 'r') as file:
                phased = json.load(file)
        except Exception:
            phased = {}

        if self.tab == 'SE':
            table = self.ui.SE_Table_Data
            self.selected_files = files
            self.phased_spectra_SE = phased
        elif self.tab == 'DQ':
            table = self.ui.table_DQ
            self.selected_files_DQ_single = files
            self.app_state.dq_files = files
            self.phased_spectra_DQ = phased
        elif self.tab == 'DQ_Temp':
            table = self.ui.table_DQ_2
            self.selected_DQfiles = files
        elif self.tab == 'T1T2':
            table = self.ui.table_T1
            self.selected_T1files = files
        elif self.tab == 'DQMQ':
            table = self.ui.DQMQ_Table_Data
            self.selected_DQMQfile = files
            self.app_state.dqmq_files = files
        elif self.tab == 'GS':
            table = self.ui.table_GS
            self.selected_GSfiles = files

        with open(file_path, 'r') as f:
            lines = f.readlines()
            if self._is_probably_old_table_format(lines):
                self._warn_old_table_format()
                return
            table.setRowCount(len(lines))
            for row, line in enumerate(lines):
                values = line.strip().split(',')
                for col, value in enumerate(values):
                    item = QTableWidgetItem(value)
                    table.setItem(row, col, item)

        if self.tab in ('SE', 'DQ'):
            if files:
                for path in files:
                    self.ui.comboBox_4.addItem(os.path.basename(path))
            else:
                for row in range(table.rowCount()):
                    self.ui.comboBox_4.addItem(f"File #{row+1}")

        if self.tab == 'SE':
            self.update_se_graphs()
            if self.ui.comboBox_4.count() > 0:
                self.ui.comboBox_4.setCurrentIndex(0)
                self.update_file()
        elif self.tab == 'DQ':
            self.dq_controller.update_graphs()
            if self.ui.comboBox_4.count() > 0:
                self.ui.comboBox_4.setCurrentIndex(0)
                self.update_file()
        elif self.tab == 'DQ_Temp':
            self.dq_temp_controller.update_DQ_comparison()
        elif self.tab == 'T1T2':
            self.t1t2_controller.update_T12_table()
        elif self.tab == 'DQMQ':
            self.ui.DQMQ_Button_PlotOriginal.setEnabled(True)
            self.ui.DQMQ_Button_PlotNorm.setEnabled(True)
            self.ui.DQMQ_DoubleSpinBox_FitFrom.setEnabled(True)
            self.ui.DQMQ_DoubleSpinBox_FitTo.setEnabled(True)
            self.ui.DQMQ_DoubleSpinBox_Power.setEnabled(True)
            self.dqmq_controller.plot_nDQ_on_Load()

        elif self.tab == 'GS':
            self.gs_controller.update_GS_table()

    def _warn_old_table_format(self):
        QMessageBox.warning(
            self,
            "Old saved file format",
            "You appear to be loading a saved file from an older version of Relaxyzer.\n\n"
            "The table column order has changed in the current version, so this file cannot be loaded safely.\n\n"
            "Please open the saved file manually, rearrange the columns to match the new table order, save it again, "
            "and then try loading it once more.",
            QMessageBox.Ok
        )

    def _is_probably_old_table_format(self, lines):
        if not lines:
            return False
        first_values = lines[0].strip().split(',')
        if self.tab == 'T1T2':
            if len(first_values) > T1Columns.FOLDER:
                return self._looks_like_path(first_values[T1Columns.X_AXIS]) or self._looks_numeric(first_values[T1Columns.FOLDER])
        if self.tab == 'GS':
            if len(first_values) > GSColumns.FOLDER:
                return self._looks_like_path(first_values[GSColumns.X_AXIS]) or self._looks_numeric(first_values[GSColumns.FOLDER])
        if self.tab == 'DQ_Temp':
            if len(first_values) > DQTempColumns.FOLDER:
                return self._looks_like_path(first_values[DQTempColumns.NAME]) and not self._looks_like_path(first_values[DQTempColumns.FOLDER])
        return False

    @staticmethod
    def _looks_numeric(value):
        try:
            float(value)
            return True
        except Exception:
            return False

    @staticmethod
    def _looks_like_path(value):
        return ('/' in value) or ('\\' in value) or value.endswith(('.txt', '.dat', '.csv'))

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
