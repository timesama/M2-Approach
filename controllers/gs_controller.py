import logging
import os

from PySide6.QtWidgets import QMessageBox, QTableWidgetItem

from calculations import gs_signal
from controllers.base_tab_controller import BaseTabController
from controllers.table_columns import GSColumns
from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


class GSTabController(BaseTabController):
    def connect_signals(self):
        self.ui.GS_Table_Results.horizontalHeader().sectionDoubleClicked.connect(
            lambda index: self.parent.renameSection(self.ui.GS_Table_Results, index)
        )
        self.ui.GS_Table_Results.itemSelectionChanged.connect(self.plot_sqrt_time)
        self.ui.GS_Button_Plot.clicked.connect(self.plot_sqrt_time_from_user)
        self.ui.GS_CheckBox_UseSqrtTime.clicked.connect(self.calculate_sqrt_time)
        self.ui.GS_RadioButton_Short.clicked.connect(self.calculate_sqrt_time)
        self.ui.GS_RadioButton_Medium.clicked.connect(self.calculate_sqrt_time)
        self.ui.GS_RadioButton_Long.clicked.connect(self.calculate_sqrt_time)
        self.ui.GS_DoubleSpinBox_FitFrom.editingFinished.connect(self.calculate_sqrt_time)
        self.ui.GS_DoubleSpinBox_FitTo.editingFinished.connect(self.calculate_sqrt_time)
        self.ui.GS_DoubleSpinBox_Beta.editingFinished.connect(self.calculate_sqrt_time)
        self.ui.GS_DoubleSpinBox_R2.editingFinished.connect(self.calculate_sqrt_time)
        self.ui.GS_DoubleSpinBox_M2.editingFinished.connect(self.calculate_sqrt_time)
        self.ui.GS_ComboBox_ChooseFile.activated.connect(lambda *_args: self.calculate_sqrt_time())

    def update_GS_table(self):
        if self.parent.selected_GSfiles:
            self._status(f"Loading {len(self.parent.selected_GSfiles)} spin diffusion file(s)...")
        with busy_cursor():
            self._update_GS_table_impl()

    def _update_GS_table_impl(self):
        selected_files = self.parent.selected_GSfiles
        if not selected_files:
            self._warn_no_data("No spin diffusion data available. Select spin diffusion files first.")
            return

        logger.info("GS loading %d files", len(selected_files))
        table = self.ui.GS_Table_Results
        combobox = self.ui.GS_ComboBox_ChooseFile
        dictionary = self.parent.GS_dictionary
        preserved_rows = self._preserved_rows(table)
        dictionary.clear()
        table.setRowCount(0)
        combobox.clear()
        failed_files = []

        for file_path in selected_files:
            current_file = os.path.basename(file_path)
            row = table.rowCount()
            logger.info("GS processing file %d/%d: %s", row + 1, len(selected_files), current_file)

            try:
                sqrt_time, short_signal, medium_signal, long_signal = gs_signal.read_spin_diffusion_file(file_path)
                x_axis = gs_signal.x_axis_from_filename(file_path, row)
            except Exception:
                logger.exception("GS file could not be read: %s", file_path)
                failed_files.append(current_file)
                continue

            preserved_row = preserved_rows.get(file_path)
            effective_x = (
                preserved_row.get(GSColumns.X_AXIS, str(x_axis))
                if preserved_row is not None
                else str(x_axis)
            )
            dictionary[file_path] = {
                "X Axis": [effective_x],
                "sqrtTime": list(sqrt_time),
                "short": list(short_signal),
                "medium": list(medium_signal),
                "long": list(long_signal),
            }
            self._add_table_row(table, combobox, row, file_path, current_file, effective_x, preserved_row)

        if failed_files:
            self._warn_failed_files("Some spin diffusion files could not be read. They were skipped.", failed_files)

        combobox.setCurrentIndex(-1)
        table.resizeColumnsToContents()
        if table.rowCount():
            self._status(f"Loaded {table.rowCount()} spin diffusion row(s).")
        else:
            self._status("Could not load spin diffusion files.")

    def _add_table_row(self, table, combobox, row, folder, file_name, default_x_axis, preserved_row=None):
        table.insertRow(row)

        if preserved_row is not None:
            for col in range(table.columnCount()):
                value = preserved_row.get(col, "")
                table.setItem(row, col, QTableWidgetItem(str(value)))

        table.setItem(row, GSColumns.FOLDER, QTableWidgetItem(folder))
        table.setItem(row, GSColumns.FILE_NAME, QTableWidgetItem(file_name))

        if preserved_row is None:
            table.setItem(row, GSColumns.X_AXIS, QTableWidgetItem(str(default_x_axis)))

        if combobox is not None:
            combobox.addItem(file_name)

    def _preserved_rows(self, table):
        preserved_rows = {}

        for row in range(table.rowCount()):
            file_item = table.item(row, GSColumns.FOLDER)
            if file_item is None:
                continue

            file_key = file_item.text().strip()
            if not file_key:
                continue

            row_data = {}
            for col in range(table.columnCount()):
                item = table.item(row, col)
                row_data[col] = item.text() if item is not None else ""

            preserved_rows[file_key] = row_data

        return preserved_rows

    def calculate_sqrt_time(self, *_args, show_warning=True):
        with busy_cursor():
            self._calculate_sqrt_time_impl(show_warning=show_warning)

    def _calculate_sqrt_time_impl(self, show_warning=True):
        table = self.ui.GS_Table_Results
        figure = self.ui.GS_PlotWidget_RawSignal
        dictionary = self.parent.GS_dictionary
        idx = self.ui.GS_ComboBox_ChooseFile.currentIndex()

        if table.rowCount() == 0 or not dictionary:
            if show_warning:
                self._warn_no_data("No spin diffusion data available. Select spin diffusion files first.")
            return

        if idx == -1:
            if show_warning:
                self._status("No rows selected.")
                QMessageBox.warning(
                    self.parent,
                    "No spin diffusion row selected",
                    "No spin diffusion row selected.",
                    QMessageBox.Ok,
                )
            return

        item = table.item(idx, GSColumns.FOLDER)
        if item is None or item.text() not in dictionary:
            if show_warning:
                self._status("Could not fit spin diffusion data: non-numeric values.")
                QMessageBox.warning(
                    self.parent,
                    "No spin diffusion data",
                    "Cannot calculate spin diffusion result because required columns contain non-numeric values.",
                    QMessageBox.Ok,
                )
            return

        key = item.text()
        dictionary_entry = dictionary[key]
        time_original = gs_signal.transformed_time(
            dictionary_entry["sqrtTime"],
            self.ui.GS_CheckBox_UseSqrtTime.isChecked(),
        )
        if len(time_original) == 0:
            if show_warning:
                self._warn_invalid_range()
            return

        self._update_fit_limits(time_original)
        signals = gs_signal.signal_arrays(dictionary_entry)
        source = self._selected_signal_source()
        signal_original = gs_signal.selected_signal(signals, source)
        fit_from = self.ui.GS_DoubleSpinBox_FitFrom.value()
        fit_to = self.ui.GS_DoubleSpinBox_FitTo.value()
        fit_time, fit_signal = gs_signal.fit_range(time_original, signal_original, fit_from, fit_to)

        if not gs_signal.is_valid_fit_range(time_original, fit_time, fit_signal):
            if show_warning:
                self._warn_invalid_range()
            return

        try:
            fit_time_curve, fit_signal_curve, sqrt_time, r2_value, diffusion_distance = gs_signal.fit_spin_diffusion(
                fit_time,
                fit_signal,
                self.ui.GS_DoubleSpinBox_Beta.value(),
                self.ui.GS_DoubleSpinBox_R2.value(),
                self.ui.GS_DoubleSpinBox_M2.value(),
            )
        except Exception:
            logger.exception("GS fit failed: index=%d", idx)
            figure.clear()
            self._status("Could not fit spin diffusion data.")
            QMessageBox.warning(
                self.parent,
                "Fitting failed",
                "Cannot fit spin diffusion data because the selected range is invalid.",
                QMessageBox.Ok,
            )
            return

        self.ui.GS_TextEdit_FitResult.setText(f"R² {r2_value}")
        table.setItem(idx, GSColumns.SQRT_TIME, QTableWidgetItem(str(sqrt_time)))
        table.setItem(idx, GSColumns.D_NM, QTableWidgetItem(str(diffusion_distance)))
        table.resizeColumnsToContents()
        figure.clear()
        figure.plot(
            time_original,
            signal_original,
            pen=None,
            symbolPen=None,
            symbol="o",
            symbolBrush="r",
            symbolSize=10,
        )
        figure.plot(fit_time_curve, fit_signal_curve, pen="b")
        dictionary_entry["sqrtT"] = sqrt_time
        dictionary_entry["d"] = diffusion_distance
        logger.info("GS fit completed: sqrtT=%s d=%s", sqrt_time, diffusion_distance)
        if show_warning:
            self._status("Fit completed.")

    def _update_fit_limits(self, time_original):
        self.ui.GS_DoubleSpinBox_FitFrom.setMinimum(time_original[0])
        self.ui.GS_DoubleSpinBox_FitFrom.setMaximum(time_original[-1])
        self.ui.GS_DoubleSpinBox_FitTo.setMinimum(time_original[min(15, len(time_original) - 1)])
        self.ui.GS_DoubleSpinBox_FitTo.setMaximum(time_original[-1])

    def _selected_signal_source(self):
        if self.ui.GS_RadioButton_Short.isChecked():
            return "short"
        if self.ui.GS_RadioButton_Medium.isChecked():
            return "medium"

        return "long"

    def plot_sqrt_time_from_user(self):
        self.plot_sqrt_time(show_warning=True)

    def plot_sqrt_time(self, show_warning=False):
        table = self.ui.GS_Table_Results
        graph = self.ui.GS_PlotWidget_SqrtTime
        graph.clear()

        if table.rowCount() < 1:
            if show_warning:
                self._status("Cannot update spin diffusion plot because the table is empty.")
                QMessageBox.warning(
                    self.parent,
                    "No spin diffusion data",
                    "Cannot update spin diffusion plot because the table is empty.",
                    QMessageBox.Ok,
                )
            return

        x_axis = []
        sqrt_time_values = []
        fallback_x_axis = 1
        invalid_result_count = 0

        for row in range(table.rowCount()):
            x_item = table.item(row, GSColumns.X_AXIS)
            sqrt_time_item = table.item(row, GSColumns.SQRT_TIME)

            try:
                x_axis.append(float(x_item.text()))
            except Exception:
                x_axis.append(fallback_x_axis)

            try:
                sqrt_time_values.append(float(sqrt_time_item.text()))
            except Exception:
                sqrt_time_values.append(0)
                invalid_result_count += 1

            fallback_x_axis += 1

        if invalid_result_count == table.rowCount() and show_warning:
            self._status("Could not plot spin diffusion data: non-numeric values.")
            QMessageBox.warning(
                self.parent,
                "No spin diffusion data",
                "Cannot calculate spin diffusion result because required columns contain non-numeric values.",
                QMessageBox.Ok,
            )
            return

        graph.plot(
            x_axis,
            sqrt_time_values,
            pen=None,
            symbolPen=None,
            symbol="o",
            symbolBrush="r",
            symbolSize=10,
        )
        self._plot_group_overlays()
        self.highlight_selected_sqrt_time_point()
        if show_warning:
            self._status("Plot updated.")

    def _plot_group_overlays(self):
        group_data = getattr(self.parent, "group_data_SD", {})
        colors = getattr(self.parent, "tab10_colors", [])
        if not group_data:
            return

        column = self._selected_gs_plot_column()

        for i, (_group_number, group_rows) in enumerate(group_data.items()):
            group_x = []
            group_y = []
            for row_data in group_rows:
                try:
                    group_x.append(float(row_data[GSColumns.X_AXIS]))
                    group_y.append(float(row_data[column]))
                except (ValueError, IndexError):
                    continue

            sorted_points = sorted(zip(group_x, group_y), key=lambda point: point[0])
            if len(sorted_points) > 1:
                xs, ys = zip(*sorted_points)
                color = colors[i % len(colors)] if colors else "r"
                self.ui.GS_PlotWidget_SqrtTime.plot(
                    xs,
                    ys,
                    pen={"color": color, "width": 2},
                    symbol="o",
                    symbolBrush=color,
                    symbolPen=None,
                    symbolSize=8,
                )

    def highlight_selected_sqrt_time_point(self):
        row = self.ui.GS_Table_Results.currentRow()
        if row < 0:
            return

        column = self._selected_gs_plot_column()
        x_item = self.ui.GS_Table_Results.item(row, GSColumns.X_AXIS)
        y_item = self.ui.GS_Table_Results.item(row, column)
        if x_item is None or y_item is None:
            return

        try:
            x = float(x_item.text())
            y = float(y_item.text())
        except ValueError:
            return

        self.ui.GS_PlotWidget_SqrtTime.plot(
            [x],
            [y],
            pen=None,
            symbol="o",
            symbolBrush=(255, 255, 0, 255),
            symbolPen="k",
            symbolSize=13,
        )

    def _selected_gs_plot_column(self):
        return GSColumns.SQRT_TIME

    def _warn_no_data(self, message):
        self._status(message)
        QMessageBox.warning(self.parent, "No spin diffusion data", message, QMessageBox.Ok)

    def _warn_invalid_range(self):
        self._status("Could not fit spin diffusion data: invalid range.")
        QMessageBox.warning(
            self.parent,
            "Invalid spin diffusion range",
            "Cannot fit spin diffusion data because the selected range is invalid.",
            QMessageBox.Ok,
        )

    def _warn_failed_files(self, message, failed_files):
        preview = "\n".join(failed_files[:5])
        if len(failed_files) > 5:
            preview += f"\n...and {len(failed_files) - 5} more"
        self._status(message)
        QMessageBox.warning(self.parent, "Spin diffusion load warning", f"{message}\n\n{preview}", QMessageBox.Ok)
