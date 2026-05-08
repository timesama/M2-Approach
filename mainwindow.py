# This Python file uses the following encoding: utf-8
import csv
import json
import logging
import os
import re
import winreg
from webbrowser import open as open_application

import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
import requests
from PySide6.QtCore import QCoreApplication, Qt
from PySide6.QtGui import QColor, QIcon
from PySide6.QtWidgets import QApplication, QFileDialog, QDialog, QInputDialog, QMainWindow, QMessageBox, QScrollArea, QTableWidgetItem
from pyqtgraph import mkColor, mkPen
from ui_Form import Ui_NMR
from controllers import (
    SETabController, DQTabController, DQTempTabController,
    T1T2TabController, DQMQTabController, GSTabController, GeneralSEDQController, RecFIDController
)
from dialogs.open_files_dialog import OpenFilesDialog
import dialogs.open_files_dialog as open_files_dialog_module
from dialogs.group_window import GroupWindow
from dialogs.save_files_dialog import SaveFilesDialog
from widgets.table_copy_enabler import TableCopyEnabler
from app_state import AppState
from controllers.table_columns import T1Columns, GSColumns, DQTempColumns
from utils.ui_busy import busy_cursor
logger = logging.getLogger(__name__)

pg.CONFIG_OPTIONS['background'] = 'w'
pg.CONFIG_OPTIONS['foreground'] = 'k'

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
        self.recfid_controller = RecFIDController(ui=self.ui, state=self.app_state, parent=self)
        self.general_se_dq_controller = GeneralSEDQController(ui=self.ui, state=self.app_state, parent=self, se_controller=self.se_controller, dq_controller=self.dq_controller)

        self.state()

        self.ui.Settings_Button_DefaultFolder.clicked.connect(self.default_folder)
        self.ui.Settings_Button_OpenProjectPage.clicked.connect(self.open_url)
        self.check_for_updates()

        self.ui.tabWidget.currentChanged.connect(self.state)

        self.ui.btn_Save.clicked.connect(self.save_data)
        self.ui.DQMQ_Button_Save.clicked.connect(self.save_data)
        self.ui.DQTemp_Button_Save.clicked.connect(self.save_data)
        self.ui.T1T2_Button_Save.clicked.connect(self.save_data)
        self.ui.GS_Button_Save.clicked.connect(self.save_data)
        if hasattr(self.ui, "pushButtonSaveRecFID"):
            self.ui.pushButtonSaveRecFID.clicked.connect(self.save_data)

        self.ui.btn_Load.clicked.connect(self.load_data)
        self.ui.DQTemp_Button_Load.clicked.connect(self.load_data)
        self.ui.T1T2_Button_Load.clicked.connect(self.load_data)
        self.ui.DQMQ_Button_Load.clicked.connect(self.load_data)
        self.ui.GS_Button_Load.clicked.connect(self.load_data)
        self.ui.btn_Phasing.clicked.connect(self.general_se_dq_controller.open_phasing_manual)

        self.ui.btn_SelectFiles.clicked.connect(self.open_select_dialog)
        self.ui.btn_Add.clicked.connect(self.add_select_dialog)
        self.ui.T1T2_Button_SelectFiles.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.DQMQ_Button_SelectFiles.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.DQTemp_Button_SelectFiles.clicked.connect(self.open_select_comparison_files_dialog)
        self.ui.GS_Button_SelectFiles.clicked.connect(self.open_select_comparison_files_dialog)

        self.ui.DQTemp_Button_ClearTable.clicked.connect(self.clear_list)
        self.ui.T1T2_Button_ClearTable.clicked.connect(self.clear_list)
        self.ui.GS_Button_ClearTable.clicked.connect(self.clear_list)
        self.ui.SE_Button_ClearTable.clicked.connect(self.clear_list)
        self.ui.DQ_Button_ClearTable.clicked.connect(self.clear_list)

        self.ui.btn_DeleteRow.clicked.connect(self.delete_row)
        self.ui.SE_Button_DeleteRow.clicked.connect(self.delete_row)
        self.ui.GS_Button_DeleteRow.clicked.connect(self.delete_row)
        self.ui.DQ_Button_DeleteRow.clicked.connect(self.delete_row)

        self.ui.btn_Start.clicked.connect(self.general_se_dq_controller.analysis)

        self.ui.DQMQ_Button_PlotOriginal.clicked.connect(
            self.dqmq_controller.plot_original
        )
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

        self.ui.T1T2_Button_Group.clicked.connect(self.open_group_window)
        self.ui.GS_Button_Group.clicked.connect(self.open_group_window)

        # Graph setup
        self.setup_graph(self.ui.FFTWidget, "Frequency, MHz", "Amplitude, a.u", "FFT")
        self.setup_graph(self.ui.FidWidget, "Time, μs", "Amplitude", "NMR Signal")
        self.setup_graph(self.ui.SE_PlotWidget_Main, "Temperature, °C", "Choose", "")
        self.setup_graph(self.ui.DQ_PlotWidget_T2, "DQ Filtering Time", "T₂*", "")
        self.setup_graph(self.ui.DQ_PlotWidget_NormIntensity, "X axis", "Norm. DQ Intensity", "")
        self.setup_graph(self.ui.DQTemp_PlotWidget_T2Distribution, "T₂*", "Norm. DQ Intensity", "FunctionFit")
        self.setup_graph(self.ui.DQTemp_PlotWidget_CenterVsXAxis, "X axis", "Center", "")
        self.setup_graph(self.ui.T1T2_PlotWidget_RawSignal, "Time, ms", "Signal", "")
        self.setup_graph(self.ui.T1T2_PlotWidget_RelaxationTime, "X axis", "τ, ms", "")
        self.setup_graph(self.ui.DQMQ_PlotWidget_Signal, "Time", "NMR signal", "")
        self.setup_graph(self.ui.DQMQ_PlotWidget_Dres, "Dres/2π, KHz", "P(Dres)", "")
        self.setup_graph(self.ui.GS_PlotWidget_RawSignal, "√Time, √us", "Signal", "")
        self.setup_graph(self.ui.GS_PlotWidget_SqrtTime, "X axis", "√Time, √us", "")
        self.setup_graph(self.ui.DQTemp_PlotWidget_PolyFit, "T₂*", "Norm. DQ Intensity", "PolyFit")
        self.recfid_controller.initialize_plots()
        self.se_controller.connect_signals()
        self.dq_controller.connect_signals()
        self.dq_temp_controller.connect_signals()
        self.t1t2_controller.connect_signals()
        self.gs_controller.connect_signals()
        self.recfid_controller.connect_signals()

        # Table Headers
        self.copy_enabler = TableCopyEnabler(self)

        self._apply_table_header_order()

        # Connect change events
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
        """Check GitHub releases and prompt when a newer app version exists."""
        current_version = '0.2.2'
        url = 'https://api.github.com/repos/timesama/M2-Approach/releases/latest'
        try:
            response = requests.get(url)
            response.raise_for_status()
            latest_release = response.json()
            latest_version = latest_release['tag_name']

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
            logger.warning("Failed to check for updates: %s", e)

    def open_url(self):
        """Open the GitHub releases page used by the Settings tab."""
        open_application('https://github.com/timesama/M2-Approach/releases')

    def update_file(self):
        """Load the selected SE/DQ file into the shared FID/FFT preview widgets."""
        # self.ui.btn_Phasing.setEnabled(True)

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

        if self.ui.Settings_CheckBox_Glycerol.isChecked() and self.ui.Settings_CheckBox_Baseline.isChecked():
            file_path_gly = self.selected_files_gly[i-1]
            file_path_empty = self.selected_files_empty[i-1]
        elif self.ui.Settings_CheckBox_Glycerol.isChecked():
            file_path_gly = self.selected_files_gly[i-1]
            file_path_empty=[]
        elif self.ui.Settings_CheckBox_Baseline.isChecked():
            file_path_gly = []
            file_path_empty = self.selected_files_empty[i-1]
        else:
            file_path_gly = []
            file_path_empty = []

        if self.ui.Settings_CheckBox_SmoothFft.isChecked():
            window_first_value = self.ui.Settings_DoubleSpinBox_SmoothWindowFrom.value()
            window_last_value = self.ui.Settings_DoubleSpinBox_SmoothWindowTo.value()
            self.window_array = np.linspace(window_first_value, window_last_value, len(files), dtype=np.int32)
        else:
            self.window_array = np.array([])

        try:
            self.general_se_dq_controller.process_file_data(file_path, file_path_gly, file_path_empty, i)
        except Exception:
            return

        if self.tab == 'DQ':
            self.highlight_row(self.ui.DQ_Table_Data, i)
            self.dq_controller.update_graphs()
        elif self.tab == 'SE':
            self.highlight_row(self.ui.SE_Table_Data, i)
            self.update_se_graphs()


    def clear_list(self):
        """Clear data, plots, and file lists for the active tab."""
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
            self.ui.DQ_Table_Data.setRowCount(0)
            self.ui.DQ_PlotWidget_T2.clear()
            self.ui.DQ_PlotWidget_NormIntensity.clear()
            self.ui.DQ_TextEdit_FitResult.setText("")
            self.ui.FFTWidget.clear()
            self.ui.FidWidget.clear()
            self.ui.btn_Start.setStyleSheet("background-color: none")
            self.phased_spectra_DQ = {}
        elif self.tab == 'DQ_Temp':
            self.selected_DQfiles = []
            self.dq_t2 = {}
            self.ui.DQTemp_Table_Results.setRowCount(0)
            self.ui.DQTemp_PlotWidget_T2Distribution.clear()
            self.ui.DQTemp_PlotWidget_CenterVsXAxis.clear()
            self.ui.DQTemp_PlotWidget_PolyFit.clear()
        elif self.tab == 'T1T2':
            self.selected_T1files = []
            self.tau_dictionary = {}
            self.ui.T1T2_Table_Results.setRowCount(0)
            self.ui.T1T2_PlotWidget_RawSignal.clear()
            self.ui.T1T2_PlotWidget_RelaxationTime.clear()
            self.group_data_T1T2 = {}
        elif self.tab == 'DQMQ':
            self.selected_DQMQfile = []
            self.app_state.dqmq_files = []
            self.dqmq_controller.reset_cached_results()
        elif self.tab == 'GS':
            self.selected_GSfiles = []
            self.GS_dictionary = {}
            self.ui.GS_Table_Results.setRowCount(0)
            self.ui.GS_PlotWidget_RawSignal.clear()
            self.ui.GS_PlotWidget_SqrtTime.clear()
            self.group_data_SD = {}
        elif self.tab == 'Extra':
            return

        if self.tab == 'T1T2':
            combobox = self.ui.T1T2_ComboBox_ChooseFile
        elif self.tab in ('SE', 'DQ'):
            combobox = self.ui.comboBox_4
        elif self.tab == 'GS':
            combobox = self.ui.GS_ComboBox_ChooseFile

        if self.tab != 'DQMQ' and  self.tab != 'DQ_Temp':
            while combobox.count()>0:
                combobox.removeItem(0)

        self.window_array = np.array([])

    def delete_row(self):
        """Delete the selected row and keep the active tab file list in sync."""
        if self.tab == 'SE':
            table = self.ui.SE_Table_Data
            combobox = self.ui.comboBox_4
            files = self.selected_files
        elif self.tab == 'DQ':
            table = self.ui.DQ_Table_Data
            combobox = self.ui.comboBox_4
            files = self.selected_files_DQ_single
        elif self.tab =='T1T2':
            table = self.ui.T1T2_Table_Results
            combobox = self.ui.T1T2_ComboBox_ChooseFile
            files = self.selected_T1files
        elif self.tab == 'GS':
            table = self.ui.GS_Table_Results
            combobox = self.ui.GS_ComboBox_ChooseFile
            files = self.selected_GSfiles
        else:
            return

        row = table.currentRow()
        if row == -1:
            if self.tab == "DQ":
                QMessageBox.warning(self, "No DQ row selected", "No DQ row selected.", QMessageBox.Ok)
            elif self.tab == "T1T2":
                QMessageBox.warning(self, "No T1/T2 row selected", "No T1/T2 row selected.", QMessageBox.Ok)
            elif self.tab == "GS":
                QMessageBox.warning(
                    self,
                    "No spin diffusion row selected",
                    "No spin diffusion row selected.",
                    QMessageBox.Ok,
                )
            else:
                QMessageBox.warning(self, "Cricket sounds", "Select the row.", QMessageBox.Ok)
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
        """Open the multi/single-file picker for comparison-style tabs."""
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
                while self.ui.T1T2_ComboBox_ChooseFile.count()>0:
                    self.ui.T1T2_ComboBox_ChooseFile.removeItem(0)
                T1fileNames = dlg.selectedFiles()
                self.selected_T1files.extend(T1fileNames)
                self.t1t2_controller.update_T12_table()
            elif self.tab == 'GS':
                while self.ui.GS_ComboBox_ChooseFile.count() > 0:
                    self.ui.GS_ComboBox_ChooseFile.removeItem(0)
                GSfileNames = dlg.selectedFiles()
                self.selected_GSfiles.extend(GSfileNames)
                self.gs_controller.update_GS_table()
            elif self.tab == 'DQMQ':
                self.selected_DQMQfile = dlg.selectedFiles()
                self.app_state.dqmq_files = self.selected_DQMQfile
                self.dqmq_controller.state.dqmq_files = self.selected_DQMQfile
                self.dqmq_controller.dq_mq_analysis()

    def open_select_dialog(self):
        """Open the primary SE/DQ file picker."""
        open_files_dialog_module.State_multiple_files = True
        dlg = OpenFilesDialog(self)

        if dlg.exec():
            self.clear_list()
            files = []
            fileNames = dlg.selectedFiles()
            files.extend(fileNames)
            # self.ui.btn_Start.setEnabled(True)
            self.ui.btn_Start.setStyleSheet("background-color: green")
            # self.ui.btn_Add.setEnabled(True)

        if self.tab == 'SE':
            self.selected_files = files
        else:
            self.selected_files_DQ_single = files
            self.app_state.dq_files = files

    def add_select_dialog(self):
        """Append files to the active SE/DQ selection."""
        open_files_dialog_module.State_multiple_files = True
        dlg = OpenFilesDialog(self)

        if self.tab == 'SE':
            files = self.selected_files
        else:
            files = self.selected_files_DQ_single

        if dlg.exec():
            fileNames = dlg.selectedFiles()
            files.extend(fileNames)
            # self.ui.btn_Start.setEnabled(True)
            self.ui.btn_Start.setStyleSheet("background-color: green")
            # self.ui.btn_Add.setEnabled(True)

        if self.tab == 'SE':
            self.selected_files = files
        else:
            self.selected_files_DQ_single = files
            self.app_state.dq_files = files

    def open_group_window(self):
        """Open the grouping dialog for tabs that support grouped plots."""
        self.group_window = GroupWindow()

        if self.tab == 'SE':
            self.group_window.copy_table_data(self.ui.SE_Table_Data)
            self.group_data_SE = {}
        elif self.tab == 'T1T2':
            if self.ui.T1T2_Table_Results.rowCount() == 0:
                QMessageBox.warning(
                    self,
                    "No T1/T2 data",
                    "No T1/T2 data available. Select T1/T2 files first.",
                    QMessageBox.Ok,
                )
                return
            self.group_window.copy_table_data(self.ui.T1T2_Table_Results)
            self.group_data_T1T2 = {}
        elif self.tab == 'GS':
            if self.ui.GS_Table_Results.rowCount() == 0:
                QMessageBox.warning(
                    self,
                    "No spin diffusion data",
                    "No spin diffusion data available. Select spin diffusion files first.",
                    QMessageBox.Ok,
                )
                return
            self.group_window.copy_table_data(self.ui.GS_Table_Results)
            self.group_data_SD = {}

        if self.group_window.exec_() == QDialog.Accepted:
            data = self.group_window.group_dict
        else:
            logger.info("Group window was cancelled.")
            return

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
        """Connect DQMQ edit signals without recalculating while users type."""
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
            "DQMQ_CheckBox_MQTailDiff",
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
        """Connect a DQMQ widget signal when the generated UI exposes it."""
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
        """Update the cached tab name when the user changes tabs."""
        current_tab_index = self.ui.tabWidget.currentIndex()
        current_widget = self.ui.tabWidget.currentWidget()
        current_name = current_widget.objectName() if current_widget is not None else ""

        if current_name == "g_RecFID":
            self.tab = 'RecFID'
        elif current_tab_index == 0:
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
        pass
        # self.ui.btn_Start.setEnabled(False)
        # self.ui.btn_Save.setEnabled(False)
        # self.ui.btn_Phasing.setEnabled(False)
        # self.ui.btn_Add.setEnabled(False)
        # self.ui.DQMQ_Button_PlotOriginal.setEnabled(False)
        # self.ui.DQMQ_Button_PlotNorm.setEnabled(False)

    def enable_buttons(self):
        pass
        # self.ui.btn_SelectFiles.setEnabled(True)
        # self.ui.btn_Start.setEnabled(True)
        # self.ui.btn_Save.setEnabled(True)
        # self.ui.radioButton.setEnabled(True)
        # self.ui.btn_Load.setEnabled(True)
        # self.ui.comboBox_4.setEnabled(True)
        # self.ui.btn_Add.setEnabled(True)

    def default_folder(self):
        """Store the Settings tab default folder in the user registry."""
        folder_path = QFileDialog.getExistingDirectory(self, "Select Default Directory")

        if folder_path:
            try:
                key = winreg.CreateKey(winreg.HKEY_CURRENT_USER, r"Software\MyApp")
                winreg.SetValueEx(key, "SelectedFolder", 0, winreg.REG_SZ, folder_path)
                winreg.CloseKey(key)
                logger.info("Default folder saved to registry: %s", folder_path)
            except Exception as e:
                logger.warning("Failed to save the default folder to the registry: %s", e)

    def renameSection(self, table, index):
        """Rename a table header and mirror it to the linked x-axis label."""
        current_header = table.horizontalHeaderItem(index).text()
        new_header, ok = QInputDialog.getText(self, "Rename Column", 
                                            f"Enter new name for the column '{current_header}':")
        if ok and new_header:
            table.horizontalHeaderItem(index).setText(new_header)
            self.update_xaxis(table, index)

    def update_xaxis(self, table, index):
        """Update the active comparison plot x-axis label from a table header."""
        if self.tab == 'SE':
            figure = self.ui.SE_PlotWidget_Main
        elif self.tab == 'T1T2':
            figure = self.ui.T1T2_PlotWidget_RelaxationTime
        elif self.tab == 'GS':
            figure = self.ui.GS_PlotWidget_SqrtTime
        elif self.tab == 'DQ_Temp':
            figure = self.ui.DQTemp_PlotWidget_CenterVsXAxis
        else:
            return

        name = table.horizontalHeaderItem(index).text()
        figure.getAxis('bottom').setLabel(name)

    def _apply_table_header_order(self):
        """Restore table headers that Qt Designer does not keep in controller order."""
        self.ui.T1T2_Table_Results.setHorizontalHeaderLabels([
            "X axis", "tau 1", "A 1", "tau 2", "A 2", "tau 3", "A 3", "File name", "Folder"
        ])
        self.ui.GS_Table_Results.setHorizontalHeaderLabels([
            "X axis", "sqrt time", "d, nm", "File name", "Folder"
        ])
        self.ui.DQTemp_Table_Results.setHorizontalHeaderLabels([
            "Name", "Center Gauss", "Center Lorenz", "Center Voigt", "Center y",
            "FWHM Gauss", "FWHM Lorenz", "FWHM Voigt", "Folder"
        ])
        self.ui.T1T2_Table_Results.resizeColumnsToContents()
        self.ui.GS_Table_Results.resizeColumnsToContents()
        self.ui.DQTemp_Table_Results.resizeColumnsToContents()

    def update_graphs(self, x, y1, y2, y3, graph):
        """Draw the shared raw time/frequency preview curves."""
        graph.clear()
        graph.plot(x, y1, pen=mkPen('k', width=3))
        graph.plot(x, y2, pen=mkPen('r', width=2))
        graph.plot(x, y3, pen=mkPen('b', width=2))

    def update_se_graphs(self):
        """Refresh the SE result plot from the SE controller."""
        self.se_controller.update_graphs()

    def write_collective_dictionary(self, dictionary, save_path_name):
        """Write dictionary list values as aligned CSV columns."""
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
        """Save the active tab results table and associated source file list."""
        if self.tab == 'SE':
            table = self.ui.SE_Table_Data
            files = self.selected_files
            default_name = os.path.split(os.path.dirname(files[0]))[1] + '_SE_MSE_FID'
            if table.rowCount() == 0 or not files:
                QMessageBox.warning(
                    self,
                    "No data loaded",
                    "Load data files first.",
                    QMessageBox.Ok,
                )
                return

        elif self.tab == 'DQ':
            table = self.ui.DQ_Table_Data
            files = self.selected_files_DQ_single
            default_name = 'Table_DQ_' + os.path.split(os.path.dirname(files[0]))[1]
            if table.rowCount() == 0 or not files:
                QMessageBox.warning(
                    self,
                    "No DQ data loaded",
                    "Load DQ data files first.",
                    QMessageBox.Ok,
                )
                return

        elif self.tab == 'DQ_Temp':
            table = self.ui.DQTemp_Table_Results
            files = self.selected_DQfiles
            if table.rowCount() == 0 or not files:
                QMessageBox.warning(
                    self,
                    "No DQ data",
                    "Load DQ comparison files first.",
                    QMessageBox.Ok,
                )
                return

            path = os.path.dirname(files[0]) + '/Table_DQ_comparison_parametrs'
            self.write_collective_dictionary(self.dq_comparison_distribution, path)
            default_name = 'Table_DQ_comparison'

        elif self.tab == 'T1T2':
            table = self.ui.T1T2_Table_Results
            files = self.selected_T1files
            if table.rowCount() == 0 or not files:
                QMessageBox.warning(
                    self,
                    "No T1/T2 data",
                    "Load T1/T2 files first.",
                    QMessageBox.Ok,
                )
                return
            default_name = os.path.split(os.path.dirname(files[0]))[1] + '_T'

        elif self.tab == 'GS':
            table = self.ui.GS_Table_Results
            files = self.selected_GSfiles
            if table.rowCount() == 0 or not files:
                QMessageBox.warning(
                    self,
                    "No spin diffusion data",
                    "Load spin diffusion files first.",
                    QMessageBox.Ok,
                )
                return
            default_name = os.path.split(os.path.dirname(files[0]))[1] + '_SpinDiffusion'

        elif self.tab == 'DQMQ':
            table = self.ui.DQMQ_Table_Data
            files = self.selected_DQMQfile

            default_name = os.path.split(os.path.dirname(files[0]))[1] + '_DQMQ_data'

            if table.rowCount() == 0 or not files:
                QMessageBox.warning(
                    self,
                    "No DQMQ data",
                    "Load DQMQ file first.",
                    QMessageBox.Ok,
                )
                return

        elif self.tab == 'RecFID':
            self.recfid_controller.save_results()
            return

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
        """Load a saved results table for the active tab."""
        dlg = OpenFilesDialog(self)
        # self.ui.btn_Save.setEnabled(False)
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
            # self.ui.comboBox_4.setEnabled(False)
            # self.ui.btn_Phasing.setEnabled(False)

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
            table = self.ui.DQ_Table_Data
            self.selected_files_DQ_single = files
            self.app_state.dq_files = files
            self.phased_spectra_DQ = phased
        elif self.tab == 'DQ_Temp':
            table = self.ui.DQTemp_Table_Results
            self.selected_DQfiles = files
        elif self.tab == 'T1T2':
            table = self.ui.T1T2_Table_Results
            self.selected_T1files = files
        elif self.tab == 'DQMQ':
            table = self.ui.DQMQ_Table_Data
            self.selected_DQMQfile = files
            self.app_state.dqmq_files = files
        elif self.tab == 'GS':
            table = self.ui.GS_Table_Results
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
            # self.ui.DQMQ_Button_PlotOriginal.setEnabled(True)
            # self.ui.DQMQ_Button_PlotNorm.setEnabled(True)
            # self.ui.DQMQ_DoubleSpinBox_FitFrom.setEnabled(True)
            # self.ui.DQMQ_DoubleSpinBox_FitTo.setEnabled(True)
            # self.ui.DQMQ_DoubleSpinBox_Power.setEnabled(True)
            self.dqmq_controller.plot_nDQ_on_Load()

        elif self.tab == 'GS':
            self.gs_controller.update_GS_table()

    def _warn_old_table_format(self):
        """Warn when a saved table appears to use an obsolete column order."""
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
        """Export FID/FFT preview images and CSV data beside source files."""
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
