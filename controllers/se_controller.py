import re

import numpy as np
from PySide6.QtWidgets import QMessageBox, QDialog, QTableWidgetItem

import Calculator as Cal
from controllers.base_tab_controller import BaseTabController
from dialogs.group_window import GroupWindow


class SETabController(BaseTabController):
    def connect_signals(self):
        self.ui.SE_Button_EAct.clicked.connect(self.plot_Arr)
        self.ui.SE_Button_Group.clicked.connect(self.open_group_window)
        self.ui.SE_Table_Data.itemSelectionChanged.connect(self.update_graphs)
        self.ui.SE_Table_Data.horizontalHeader().sectionDoubleClicked.connect(
            lambda index: self.parent.renameSection(self.ui.SE_Table_Data, index)
        )
        self.ui.SE_ComboBox_YAxis.activated.connect(self.update_graphs_from_user)
        self.ui.SE_ComboBox_YAxis.currentIndexChanged.connect(
            lambda _index: self.ui.SE_Table_Data.clearSelection()
        )
        self.ui.SE_ComboBox_YAxis.currentIndexChanged.connect(lambda _index: self.update_graphs())
        self.ui.SE_RadioButton_EActKCal.clicked.connect(self.plot_Arr)
        self.ui.SE_RadioButton_EActKJ.clicked.connect(self.plot_Arr)
        self.ui.SE_CheckBox_EActToKelvin.clicked.connect(self.plot_Arr)
        self.ui.SE_DoubleSpinBox_EActStart.editingFinished.connect(self.plot_Arr)
        self.ui.SE_DoubleSpinBox_EActEnd.editingFinished.connect(self.plot_Arr)
        self.ui.SE_Button_EActDone.clicked.connect(self.hide_Eact)

    def process_processed_file(self, i, filename, amp, m2, t2, file_path):
        match = re.search(r".*_(-?\s*\d+\.?\d*).*.dat", filename)
        temperature = match.group(1) if match else "0"

        short_from = self.ui.Settings_DoubleSpinBox_SCShortStart.value() * 2
        short_to = self.ui.Settings_DoubleSpinBox_SCShortEnd.value() * 2
        long_from = self.ui.Settings_DoubleSpinBox_SCLongStart.value() * 2
        long_to = self.ui.Settings_DoubleSpinBox_SCLongEnd.value() * 2
        absolute = self.ui.Settings_RadioButton_SCAbsolute.isChecked()

        times = [int(short_from), int(short_to), int(long_from), int(long_to)]
        sfc = Cal.calculate_SC(amp, times, absolute)

        self.ui.SE_Table_Data.setRowCount(i)
        self.parent.fill_table(self.ui.SE_Table_Data, temperature, sfc, m2, t2, i)
        self.ui.SE_Table_Data.setItem(i - 1, 4, QTableWidgetItem(filename))

        if self.ui.radioButton.isChecked():
            self.parent.save_figures(file_path, filename)
        self.update_graphs()

    def update_graphs_from_user(self):
        if self.ui.SE_Table_Data.rowCount() == 0:
            self._status("No SE data available.")
            QMessageBox.warning(
                self.parent,
                "No SE data",
                "No SE data available. Load/analyze SE files first.",
                QMessageBox.Ok,
            )
            return

        self.update_graphs()

    def update_graphs(self):
        table = self.ui.SE_Table_Data
        text = self.ui.SE_ComboBox_YAxis.currentText()

        if table.rowCount() == 0:
            self.ui.SE_PlotWidget_Main.clear()
            return

        if not text:
            self.ui.SE_PlotWidget_Main.getAxis("left").setLabel("Not Set")
            return

        x = self.read_column_values(table, 0)
        if text == "SC":
            y = self.read_column_values(table, 1)
            self.ui.SE_PlotWidget_Main.getAxis("left").setLabel("SC")
        elif text == "M₂":
            y = self.read_column_values(table, 2)
            self.ui.SE_PlotWidget_Main.getAxis("left").setLabel("M₂")
        elif text == "T₂*":
            y = self.read_column_values(table, 3)
            self.ui.SE_PlotWidget_Main.getAxis("left").setLabel("T₂*")
        else:
            self.ui.SE_PlotWidget_Main.getAxis("left").setLabel("Not Set")
            self._status("Cannot update SE plot: selected Y column is missing.")
            QMessageBox.warning(
                self.parent,
                "SE plot unavailable",
                "Cannot update SE plot because the selected Y column is missing.",
                QMessageBox.Ok,
            )
            return

        self.ui.SE_PlotWidget_Main.clear()
        self.ui.SE_PlotWidget_Main.plot(
            x,
            y,
            pen=None,
            symbol="o",
            symbolPen=None,
            symbolBrush=(255, 0, 0, 255),
            symbolSize=10,
        )

        self._plot_group_overlays()
        self.highlight_selected_point()

    def _plot_group_overlays(self):
        group_data = getattr(self.parent, "group_data_SE", {})
        colors = getattr(self.parent, "tab10_colors", [])
        if not group_data:
            return

        y_col = self.ui.SE_ComboBox_YAxis.currentIndex() + 1
        for i, (_group_number, group_rows) in enumerate(group_data.items()):
            group_x = []
            group_y = []
            for row_data in group_rows:
                try:
                    group_x.append(float(row_data[0]))
                    group_y.append(float(row_data[y_col]))
                except (ValueError, IndexError):
                    continue

            sorted_points = sorted(zip(group_x, group_y), key=lambda point: point[0])
            if len(sorted_points) > 1:
                xs, ys = zip(*sorted_points)
                color = colors[i % len(colors)] if colors else "r"
                self.ui.SE_PlotWidget_Main.plot(
                    xs,
                    ys,
                    pen={"color": color, "width": 2},
                    symbol="o",
                    symbolBrush=color,
                    symbolPen=None,
                    symbolSize=8,
                )

    def plot_Arr(self):
        self.ui.SE_GroupBox_EAct.setHidden(False)
        table = self.ui.SE_Table_Data
        graph = self.ui.SE_PlotWidget_Main

        if table.rowCount() == 0:
            self._status("No SE data available.")
            QMessageBox.warning(
                self.parent,
                "No SE data",
                "Cannot calculate activation energy because the SE table is empty.",
                QMessageBox.Ok,
            )
            return

        temperature = np.array(self.read_column_values(table, 0))
        t2 = np.array(self.read_column_values(table, 3))
        sorted_indices = np.argsort(temperature)
        temperature = temperature[sorted_indices]
        t2 = t2[sorted_indices]

        if self.ui.SE_CheckBox_EActToKelvin.isChecked():
            temperature = temperature + 273.15

        starting_point = int(self.ui.SE_DoubleSpinBox_EActStart.value())
        ending_point = -int(self.ui.SE_DoubleSpinBox_EActEnd.value())

        if ending_point == 0:
            ending_point = None

        try:
            temperature = temperature[starting_point:ending_point]
            t2 = t2[starting_point:ending_point]

            valid_mask = np.isfinite(temperature) & np.isfinite(t2) & (t2 > 0)
            temperature = temperature[valid_mask]
            t2 = t2[valid_mask]

            if len(temperature) < 2:
                raise ValueError("Not enough valid points for Arrhenius fit (need at least 2 positive T₂* values).")

            reciprocal_temperature, ln_t2 = Cal.calculate_Arrhenius_ax(temperature, t2)
            temp_fit, fitted_curve, eact, r2 = Cal.calculate_Eact(
                reciprocal_temperature,
                ln_t2,
                self.ui.SE_RadioButton_EActKJ.isChecked(),
            )

            graph.clear()
            graph.plot(
                reciprocal_temperature,
                ln_t2,
                pen=None,
                symbol="o",
                symbolPen=None,
                symbolBrush=(255, 0, 0, 255),
                symbolSize=10,
            )
            graph.plot(temp_fit, fitted_curve, pen="b")
            self.parent.setup_graph(graph, "1000/T, 𝐾⁻¹", "ln(τ)", "")
            self.ui.SE_TextEdit_EActResult.setText(f"Eact = {eact}\nR² {r2}")
            self._status("Fit completed.")
        except Exception as exc:
            self._status(f"Could not fit Arrhenius plot: {exc}")
            QMessageBox.warning(self.parent, "Arrhenius fit error", str(exc), QMessageBox.Ok)

    def hide_Eact(self):
        self.ui.SE_GroupBox_EAct.setHidden(True)
        self.ui.SE_PlotWidget_Main.clear()
        self.parent.setup_graph(self.ui.SE_PlotWidget_Main, "", "", "")
        self.parent.update_xaxis(self.ui.SE_Table_Data, 0)
        self._status("Activation energy view closed.")

    def open_group_window(self):
        table = self.ui.SE_Table_Data
        if table.rowCount() == 0:
            self._status("No SE data available.")
            QMessageBox.warning(
                self.parent,
                "No SE data",
                "No SE data available. Load/analyze SE files first.",
                QMessageBox.Ok,
            )
            return

        group_window = GroupWindow()
        group_window.copy_table_data(table)
        self.parent.group_data_SE = {}

        if group_window.exec_() != QDialog.Accepted:
            return

        self.parent.group_data_SE = group_window.group_dict
        self.update_graphs()
        self._status("Updated SE groups.")

    def highlight_selected_point(self):
        row = self.ui.SE_Table_Data.currentRow()
        if row < 0:
            return

        x_item = self.ui.SE_Table_Data.item(row, 0)
        y_col = self.ui.SE_ComboBox_YAxis.currentIndex() + 1
        y_item = self.ui.SE_Table_Data.item(row, y_col)
        if x_item is None or y_item is None:
            return

        try:
            x = float(x_item.text())
            y = float(y_item.text())
        except ValueError:
            return

        self.ui.SE_PlotWidget_Main.plot(
            [x],
            [y],
            pen=None,
            symbol="o",
            symbolBrush=(255, 255, 0, 255),
            symbolPen="k",
            symbolSize=13,
        )
