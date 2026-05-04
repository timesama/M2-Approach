from controllers.base_tab_controller import BaseTabController
import numpy as np
from PySide6.QtWidgets import QMessageBox
import Calculator as Cal


class SETabController(BaseTabController):
    def process_processed_file(self, i, filename, amp, m2, t2, file_path):
        match = __import__("re").search(r'.*_(-?\s*\d+\.?\d*).*.dat', filename)
        temperature = match.group(1) if match else '0'
        short_from = self.ui.SC_short.value() * 2
        short_to = self.ui.SC_short_2.value() * 2
        long_from = self.ui.SC_long.value() * 2
        long_to = self.ui.SC_long_2.value() * 2
        absolute = self.ui.radioButton_absolute.isChecked()
        times = [int(short_from), int(short_to), int(long_from), int(long_to)]
        sfc = Cal.calculate_SC(amp, times, absolute)
        self.ui.table_SE.setRowCount(i)
        self.parent.fill_table(self.ui.table_SE, temperature, sfc, m2, t2, i)
        self.ui.table_SE.setItem(i - 1, 4, __import__("PySide6.QtWidgets", fromlist=["QTableWidgetItem"]).QTableWidgetItem(filename))
        if self.ui.radioButton.isChecked():
            self.parent.save_figures(file_path, filename)
        self.update_graphs()

    def update_graphs(self):
        x = self.read_column_values(self.ui.table_SE, 0)
        text = self.ui.comboBox_SE_chooseY.currentText()

        if text == "SC":
            y = self.read_column_values(self.ui.table_SE, 1)
            self.ui.SEWidget.getAxis('left').setLabel("SC")
        elif text == "M₂":
            y = self.read_column_values(self.ui.table_SE, 2)
            self.ui.SEWidget.getAxis('left').setLabel("M₂")
        elif text == "T₂*":
            y = self.read_column_values(self.ui.table_SE, 3)
            self.ui.SEWidget.getAxis('left').setLabel("T₂*")
        else:
            self.ui.SEWidget.getAxis('left').setLabel("Not Set")
            return

        self.ui.SEWidget.clear()
        self.ui.SEWidget.plot(
            x, y, pen=None,
            symbol='o', symbolPen=None,
            symbolBrush=(255, 0, 0, 255), symbolSize=10
        )

        group_data = getattr(self.parent, "group_data_SE", {})
        colors = getattr(self.parent, "tab10_colors", [])
        if group_data:
            y_col = self.ui.comboBox_SE_chooseY.currentIndex() + 1
            for i, (_group_number, group_rows) in enumerate(group_data.items()):
                group_x, group_y = [], []
                for row_data in group_rows:
                    try:
                        group_x.append(float(row_data[0]))
                        group_y.append(float(row_data[y_col]))
                    except (ValueError, IndexError):
                        continue
                sorted_points = sorted(zip(group_x, group_y), key=lambda p: p[0])
                if len(sorted_points) > 1:
                    xs, ys = zip(*sorted_points)
                    color = colors[i % len(colors)] if colors else 'r'
                    self.ui.SEWidget.plot(
                        xs, ys,
                        pen={'color': color, 'width': 2},
                        symbol='o', symbolBrush=color, symbolPen=None, symbolSize=8
                    )

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
            Temperature = Temperature + 273.15
        starting_point = int(self.ui.Eact_start.value())
        ending_point = -(int(self.ui.Eact_end.value()))
        if ending_point == 0:
            ending_point = None
        try:
            Temperature = Temperature[starting_point:ending_point]
            T2 = T2[starting_point:ending_point]

            valid_mask = np.isfinite(Temperature) & np.isfinite(T2) & (T2 > 0)
            Temperature = Temperature[valid_mask]
            T2 = T2[valid_mask]

            if len(Temperature) < 2:
                raise ValueError("Not enough valid points for Arrhenius fit (need at least 2 positive T₂* values).")

            reciprocal_temperature, lnT2 = Cal.calculate_Arrhenius_ax(Temperature, T2)
            Temp_fit, fitted_curve, Eact, R2 = Cal.calculate_Eact(reciprocal_temperature, lnT2, self.ui.radioButton_8.isChecked())
            graph.clear()
            graph.plot(reciprocal_temperature, lnT2, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 255), symbolSize=10)
            graph.plot(Temp_fit, fitted_curve, pen='b')
            self.parent.setup_graph(graph, "1000/T, 𝐾⁻¹", "ln(τ)", "")
            self.ui.textEdit_EAct.setText(f"Eact = {Eact}\nR² {R2}")
        except Exception as e:
            QMessageBox.warning(self.parent, "Arrhenius fit error", str(e), QMessageBox.Ok)

    def hide_Eact(self):
        self.ui.groupBox_EAct.setHidden(True)
        self.ui.SEWidget.clear()
        self.parent.setup_graph(self.ui.SEWidget, "", "", "")
        self.parent.update_xaxis(self.ui.table_SE, 0)
