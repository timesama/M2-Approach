# controller for fitting FID

import logging
import os
import re

import numpy as np
import pandas as pd

from PySide6.QtWidgets import QMessageBox, QHBoxLayout, QLabel, QDoubleSpinBox, QCheckBox, QTableWidgetItem

from scipy.optimize import curve_fit
from PySide6.QtCore import Qt
from pyqtgraph import mkPen

import Calculator as Cal
from calculations import modelfit_signal
from controllers.base_tab_controller import BaseTabController

from utils.ui_busy import busy_cursor

logger = logging.getLogger(__name__)


MODELFIT_PLOT_STYLES = {
    "Original": {
        "color": (0, 0, 0),
        "width": 4,
    },
    "Cumulative": {
        "color": (255, 0, 0),
        "width": 3,
    },
    "Contribution1": {
        "color": (0, 255, 255),
        "width": 3,
    },
    "Contribution2": {
        "color": (255, 170, 0),
        "width": 3,
    },
    "Contribution3": {
        "color": (170, 170, 255),
        "width": 3,
    },
    "Contribution4": {
        "color": (170, 255, 127),
        "width": 2,
    },
}

class ModelFitController(BaseTabController):

    def __init__(self, ui, state, parent=None):
        super().__init__(ui, state, parent)
        self.file_settings = {}
        self.active_file_path = None

        self.parameter_widgets = {
            1: {},
            2: {},
            3: {},
            4: {},
        }

    def apply_plot_style_to_checkboxes(self):
        label_styles = {
            "ModelFit_Label_1": "Contribution1",
            "ModelFit_Label_2": "Contribution2",
            "ModelFit_Label_3": "Contribution3",
            "ModelFit_Label_4": "Contribution4",
        }

        gradient = []

        for widget_name, style_name in label_styles.items():
            widget = getattr(self.ui, widget_name, None)
            if widget is None:
                continue
            color = MODELFIT_PLOT_STYLES[style_name]["color"]
            widget.setStyleSheet(f"background-color: rgb{color}; color: black;")

            color_list = list(color)
            color_list.insert(3, 255)
            color_tuple = tuple(color_list)

            gradient.append(color_tuple)

        checkbox_styles = {
            "ModelFit_CheckBox_Show_Original": "Original",
            "ModelFit_CheckBox_Show_Cumulative": "Cumulative",
        }

        for widget_name, style_name in checkbox_styles.items():
            widget = getattr(self.ui, widget_name, None)
            if widget is None:
                continue
            color = MODELFIT_PLOT_STYLES[style_name]["color"]
            widget.setStyleSheet(f"color: rgb{color};")

        cumulative_checkbox = self.ui.ModelFit_CheckBox_Show_Contributions

        cumulative_text_gradient = f"color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba{gradient[0]}, stop:0.33 rgba{gradient[1]}, stop:0.66 rgba{gradient[2]}, stop:1 rgba{gradient[3]});"
        cumulative_checkbox.setStyleSheet(cumulative_text_gradient)



        # color = MODELFIT_PLOT_STYLES[style_name]["color"]

    def connect_signals(self):
        self.apply_plot_style_to_checkboxes()
        self.initialize_model_rows()
        self.connect_model_rows()

        self.ui.ModelFit_ComboBox_ChooseFile.activated.connect(self.on_file_selected)

        self.ui.ModelFit_Button_Fit.clicked.connect(self.run_fitting)

        self.ui.ModelFit_CheckBox_Show_Original.toggled.connect(self.plot)
        self.ui.ModelFit_CheckBox_Show_Cumulative.toggled.connect(self.plot)
        self.ui.ModelFit_CheckBox_Show_Contributions.toggled.connect(self.plot)

### when the model is chosen
    def initialize_model_rows(self):

        model_names = list(modelfit_signal.MODEL_SPECS.keys())

        comboboxes = [
            self.ui.ModelFit_Combobox_Function_1,
            self.ui.ModelFit_Combobox_Function_2,
            self.ui.ModelFit_Combobox_Function_3,
            self.ui.ModelFit_Combobox_Function_4,
        ]

        for combo in comboboxes:
            combo.clear()
            combo.addItems(model_names)
            combo.setCurrentIndex(-1)

    def connect_model_rows(self):

        for row in range(1, 5):
            combo = getattr(self.ui, f"ModelFit_Combobox_Function_{row}")

            combo.currentTextChanged.connect(
                lambda model_name, row=row:
                    self.on_model_changed(row, model_name)
            )

    def on_model_changed(self, row, model_name):

        equation_label = getattr(self.ui, f"ModelFit_Label_Equation_{row}")

        if not model_name:
            equation_label.setText("Equation")
            self.build_parameter_widgets(row, "")
            return

        spec = modelfit_signal.MODEL_SPECS.get(model_name)

        if spec is None:
            equation_label.setText("Equation")
            return

        equation_label.setText(spec["equation"])

        self.build_parameter_widgets(row, model_name)

    def build_parameter_widgets(self, row, model_name):

        container = getattr(self.ui, f"ModelFit_Widget_Equation_{row}")

        layout = container.layout()

        if layout is None:
            layout = QHBoxLayout(container)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(4)

        # Delete whatever model was displayed previously
        while layout.count():
            item = layout.takeAt(0)

            widget = item.widget()

            if widget is not None:
                widget.deleteLater()

        self.parameter_widgets[row] = {}

        if not model_name:
            return

        spec = modelfit_signal.MODEL_SPECS.get(model_name)

        if spec is None:
            return

        for parameter in spec["parameters"]:

            parameter_name = parameter["name"]

            # parameter symbol
            label = QLabel(parameter["label"])

            # value / initial guess
            spinbox = QDoubleSpinBox()
            spinbox.setFixedSize(100, 30)
            spinbox.setDecimals(3)

            spinbox.setRange(parameter["minimum"], parameter["maximum"])
            spinbox.setValue(parameter["default"])

            # fixed/free switch
            fixed_checkbox = QCheckBox()
            fixed_checkbox.setText("")
            fixed_checkbox.setFixedSize(30, 30)
            fixed_checkbox.setToolTip("Fix Value")

            layout.addWidget(label)
            layout.addWidget(spinbox)
            layout.addWidget(fixed_checkbox)

            self.parameter_widgets[row][parameter_name] = {
                "label": label,
                "spinbox": spinbox,
                "fixed": fixed_checkbox,
            }

        layout.addStretch()

#### When there are files
    def on_files_load(self):
        self.populate_combobox()

    def populate_combobox(self):

        combo = self.ui.ModelFit_ComboBox_ChooseFile
        while combo.count() > 0:
            combo.removeItem(0)

        for file_path in self.parent.selected_ModelFit_files:
            filename = os.path.basename(file_path)
            combo.addItem(filename, file_path)

            logger.info("Update files list: %s file", filename)

        combo.setCurrentIndex(-1)

### On file selected
    def on_file_selected(self):

        selected_file_path = self.ui.ModelFit_ComboBox_ChooseFile.currentData()
        if not selected_file_path:
            return

        if self.active_file_path and self.active_file_path != selected_file_path:
            self.file_settings[self.active_file_path] = self._capture_current_settings()

        self.active_file_path = selected_file_path

        selected_filename = os.path.basename(selected_file_path)
        logger.info("Open file: %s", selected_filename)

        if not os.path.isfile(selected_file_path):
            QMessageBox.warning(
                self.parent,
                "Can't open the file",
                "The selected file no longer exists.",
                QMessageBox.Ok,
            )

            logger.warning(
                "File no longer exists: %s",
                selected_file_path
            )

            self._status("Selected file does not exist.")
            return

        dictionary = self.parent.ModelFit_dictionary

        with busy_cursor():
            self._restore_current_settings()

            # regular pre-processing
            try:
                if selected_file_path not in dictionary:
                    Time, Real_Signal = self.run_analysis(selected_file_path)

                    dictionary[selected_file_path] = {
                        "Time": np.asarray(Time),
                        "Signal": np.asarray(Real_Signal),
                    }
                else:
                    Time = np.asarray(dictionary[selected_file_path]["Time"])
                    Real_Signal = np.asarray(dictionary[selected_file_path]["Signal"])

                self.plot()

            except Exception as exc:
                logger.exception(
                    "Preprocessing failed for %s",
                    selected_file_path
                )

                QMessageBox.warning(
                    self.parent,
                    "Can't open the file",
                    "Pre-processing failed for the selected file.",
                    QMessageBox.Ok,
                )

                self._status("Preprocessing failed.")
                return

### Preprocessing on file open
    def run_analysis(self, file_path):
        """ Run preprocessing for poor with double frequency adjustment. Return time and real signal """
        time, re_signal, im_signal = Cal.analysis_time_domain(file_path, [], False)

        frequency = Cal._calculate_frequency_scale(time)
        real_adjusted, _ = Cal._adjust_frequency(frequency, re_signal, im_signal)

        return time, real_adjusted

#### Plot and populate the table
    def plot(self):

        graph = self.ui.ModelFit_PlotWidget_PlotData
        graph.clear()

        file_path = self.ui.ModelFit_ComboBox_ChooseFile.currentData()

        if not file_path:
            return

        data = self.parent.ModelFit_dictionary.get(file_path)

        if data is None:
            return

        time = np.asarray(data["Time"])
        signal = np.asarray(data["Signal"])

        if self.ui.ModelFit_CheckBox_Show_Original.isChecked():
            graph.plot(
                    time,
                    signal,
                    pen=self._modelfit_pen("Original"),
                    name="Original"
                )

        if (self.ui.ModelFit_CheckBox_Show_Cumulative.isChecked() and "Cumulative" in data):
            graph.plot(
                time,
                np.asarray(data["Cumulative"]),
                pen=self._modelfit_pen("Cumulative"),
                name="Cumulative",
            )

        contribution_checkbox = getattr(self.ui, "ModelFit_CheckBox_Show_Contributions", None)

        show_contributions = (contribution_checkbox is None or contribution_checkbox.isChecked())

        if show_contributions:

            for row in range(1, 5):

                key = f"Contribution{row}"

                if key not in data:
                    continue

                graph.plot(
                    time,
                    np.asarray(data[key]),
                    pen=self._modelfit_pen(
                        key,
                        dashed=True
                    ),
                    name=key,
                )

    def _modelfit_pen(self, style_name, *, dashed=False):
        style = MODELFIT_PLOT_STYLES[style_name]
        qt_style = Qt.DashLine if dashed else Qt.SolidLine

        return mkPen(
            style["color"],
            width=style["width"],
            style=qt_style,
        )

    def _populate_fit_table(self, parameter_results, dependencies, r_squared, reduced_chi_squared, amplitude_ratios):

        table = self.ui.ModelFit_Table_Data
        table.clearContents()
        number_rows = (len(parameter_results) + 6)
        table.setRowCount(number_rows)

        row_index = 0

        for result in parameter_results:

            key = (result["row"], result["name"])

            parameter_text = (
                f"{result['row']} "
                f"{result['model']} "
                f"{result['label']}"
                )

            table.setItem(row_index, 0, QTableWidgetItem(parameter_text))
            table.setItem(row_index, 1, QTableWidgetItem(f"{result['value']:.6g}"))

            if result["fixed"]:
                error_text = "—"

            elif np.isfinite(result["error"]):
                error_text = (f"{result['error']:.6g}")

            else:
                error_text = "NaN"

            table.setItem(row_index, 2, QTableWidgetItem(error_text))
            table.setItem(row_index, 3, QTableWidgetItem(dependencies[key]))

            row_index += 1

        for component_row in range(1, 5):
            table.setItem(row_index, 0, QTableWidgetItem(f"A{component_row} ratio"))
            table.setItem(row_index, 1, QTableWidgetItem(f"{amplitude_ratios[component_row]:.6g}"))
            row_index += 1

        # R2
        table.setItem(row_index, 0, QTableWidgetItem("R²"))
        table.setItem(row_index, 1, QTableWidgetItem(f"{r_squared:.6g}"))

        row_index += 1

        # reduced chi squared
        table.setItem(row_index, 0, QTableWidgetItem("Reduced χ²"))
        table.setItem(row_index, 1, QTableWidgetItem(f"{reduced_chi_squared:.6g}"))

        table.resizeColumnsToContents()

#### Save settings for file - restore -reset etc
    def _restore_current_settings(self):
        file_path = self.ui.ModelFit_ComboBox_ChooseFile.currentData()

        self.reset_fit_ui()

        settings = self.file_settings.get(file_path)
        if settings is None:
            return

        for row, row_settings in enumerate(settings["rows"], start=1):
            use_box = getattr(self.ui, f"ModelFit_CheckBox_UseFunction_{row}")
            combo = getattr(self.ui, f"ModelFit_Combobox_Function_{row}")

            use_box.setChecked(row_settings["enabled"])

            model_name = row_settings["model"]

            if model_name:
                combo.setCurrentText(model_name)
                self.on_model_changed(row, model_name)

                for name, parameter in row_settings["parameters"].items():
                    widgets = self.parameter_widgets[row].get(name)

                    if widgets is None:
                        continue

                    widgets["spinbox"].setValue(parameter["value"])
                    widgets["fixed"].setChecked(parameter["fixed"])

        plot = settings["plot"]

        self.ui.ModelFit_CheckBox_Show_Original.setChecked(plot["original"])
        self.ui.ModelFit_CheckBox_Show_Cumulative.setChecked(plot["cumulative"])
        self.ui.ModelFit_CheckBox_Show_Contributions.setChecked(plot["contributions"])

        self._restore_table(settings["table"])

    def _restore_table(self, rows):
        table = self.ui.ModelFit_Table_Data

        table.clearContents()
        table.setRowCount(len(rows))

        for row_index, row_data in enumerate(rows):
            for column_index, value in enumerate(row_data):
                table.setItem(row_index, column_index, QTableWidgetItem(value))

        table.resizeColumnsToContents()

    def _save_current_settings(self):
        file_path = self.ui.ModelFit_ComboBox_ChooseFile.currentData()

        if file_path:
            self.file_settings[file_path] = self._capture_current_settings()

    def _capture_current_settings(self):
        rows = []

        for row in range(1, 5):
            combo = getattr(self.ui, f"ModelFit_Combobox_Function_{row}")
            use_box = getattr(self.ui, f"ModelFit_CheckBox_UseFunction_{row}")

            parameters = {}

            for name, widgets in self.parameter_widgets[row].items():
                parameters[name] = {
                    "value": widgets["spinbox"].value(),
                    "fixed": widgets["fixed"].isChecked(),
                }

            rows.append({
                "enabled": use_box.isChecked(),
                "model": combo.currentText(),
                "parameters": parameters,
            })

        return {
            "rows": rows,
            "plot": {
                "original": self.ui.ModelFit_CheckBox_Show_Original.isChecked(),
                "cumulative": self.ui.ModelFit_CheckBox_Show_Cumulative.isChecked(),
                "contributions": self.ui.ModelFit_CheckBox_Show_Contributions.isChecked(),
            },
            "table": self._capture_table(),
        }

    def _capture_table(self):
        table = self.ui.ModelFit_Table_Data

        return [
            [
                table.item(row, column).text() if table.item(row, column) else ""
                for column in range(table.columnCount())
            ]
            for row in range(table.rowCount())
        ]

    def reset_fit_ui(self):
        for row in range(1, 5):
            use_box = getattr(self.ui, f"ModelFit_CheckBox_UseFunction_{row}")
            combo = getattr(self.ui, f"ModelFit_Combobox_Function_{row}")

            use_box.setChecked(False)
            combo.setCurrentIndex(-1)
            self.build_parameter_widgets(row, "")

            getattr(self.ui, f"ModelFit_Label_Equation_{row}").setText("Equation")

        self.ui.ModelFit_CheckBox_Show_Original.setChecked(True)
        self.ui.ModelFit_CheckBox_Show_Cumulative.setChecked(False)
        self.ui.ModelFit_CheckBox_Show_Contributions.setChecked(False)

        self.ui.ModelFit_Table_Data.clearContents()
        self.ui.ModelFit_Table_Data.setRowCount(0)

### Actual fit
    def run_fitting(self):
        file_path = self.ui.ModelFit_ComboBox_ChooseFile.currentData()

        if not file_path:
            QMessageBox.warning(self.parent, "No data", "Choose a file first.", QMessageBox.Ok)
            return

        data = self.parent.ModelFit_dictionary.get(file_path)
        if data is None:
            return

        try:
            components = self._collect_components()

            if not components:
                raise ValueError("Enable at least one fitting function.")

            result = modelfit_signal.fit_components(
                np.asarray(data["Time"], dtype=float),
                np.asarray(data["Signal"], dtype=float),
                components,
            )

            self._apply_fit_result(file_path, result)

        except Exception as exc:
            logger.exception("ModelFit fitting failed")
            QMessageBox.warning(self.parent, "Fitting failed", str(exc), QMessageBox.Ok)
            self._status("Model fitting failed.")

    def _collect_components(self):

        components = []

        for row in range(1, 5):

            include_checkbox = getattr(self.ui, f"ModelFit_CheckBox_UseFunction_{row}")

            if not include_checkbox.isChecked():
                continue

            combo = getattr(self.ui, f"ModelFit_Combobox_Function_{row}")

            model_name = combo.currentText()

            if not model_name:
                raise ValueError(f"Function {row} is enabled but no model is selected.")

            spec = modelfit_signal.MODEL_SPECS[model_name]

            parameters = []

            for parameter_spec in spec["parameters"]:

                name = parameter_spec["name"]

                widgets = self.parameter_widgets[row][name]

                parameters.append({
                    "name": name,
                    "label": parameter_spec["label"],
                    "value": widgets["spinbox"].value(),
                    "fixed": widgets["fixed"].isChecked(),
                    "minimum": parameter_spec["minimum"],
                    "maximum": parameter_spec["maximum"],
                })

            components.append({
                "row": row,
                "model": model_name,
                "function": spec["function"],
                "parameters": parameters,
            })

        return components

    def _apply_fit_result(self, file_path, result):
        data = self.parent.ModelFit_dictionary[file_path]

        # Store fitted curves
        data["Cumulative"] = np.asarray(result["cumulative"])

        for row in range(1, 5):
            data.pop(f"Contribution{row}", None)

        for row, contribution in result["contributions"].items():
            data[f"Contribution{row}"] = np.asarray(contribution)

        # Store fit metadata/results
        data["Fit"] = {
            "success": True,
            "parameters": result["parameters"],
            "covariance": result["covariance"],
            "correlation": result["correlation"],
            "dependency": result["dependency"],
            "r_squared": result["r_squared"],
            "reduced_residual_variance": result["reduced_residual_variance"],
            "amplitude_ratios": result["amplitude_ratios"],
        }

        # Update dynamic parameter spinboxes with fitted values
        for parameter in result["parameters"]:
            row = parameter["row"]
            name = parameter["name"]

            widgets = self.parameter_widgets[row].get(name)

            if widgets is None:
                continue

            widgets["spinbox"].setValue(parameter["value"])

        # Update results table
        self._populate_fit_table(
            result["parameters"],
            result["dependency"],
            result["r_squared"],
            result["reduced_residual_variance"],
            result["amplitude_ratios"],
        )

        # Remember the complete GUI state for this file
        self._save_current_settings()

        # Refresh plot with cumulative/contributions
        self.plot()

        self._status("Model fitting completed.")

### Save and delete
    def delete_file(self):
        combo = self.ui.ModelFit_ComboBox_ChooseFile
        index = combo.currentIndex()
        file_path = combo.currentData()

        if index < 0 or not file_path:
            QMessageBox.warning(self.parent, "No file selected", "Select a file to delete.", QMessageBox.Ok)
            return

        filename = os.path.basename(file_path)

        # Remove numerical data/results
        self.parent.ModelFit_dictionary.pop(file_path, None)

        # Remove remembered GUI settings
        self.file_settings.pop(file_path, None)

        # Remove from loaded-file list
        if file_path in self.parent.selected_ModelFit_files:
            self.parent.selected_ModelFit_files.remove(file_path)

        # Remove from combobox
        combo.removeItem(index)
        combo.setCurrentIndex(-1)

        # Nothing is active anymore
        self.active_file_path = None

        # Reset GUI
        self.reset_fit_ui()
        self.ui.ModelFit_PlotWidget_PlotData.clear()

        logger.info("Deleted ModelFit file: %s", filename)
        self._status(f"Deleted {filename}.")

    def save_results_excel(self, base_file_path):
        """
        Save the currently selected ModelFit result into one Excel workbook.

        Sheet order:
            1. Original Data
            2. Fitting Data
            3. Fitting Metadata
        """
        file_path = self.ui.ModelFit_ComboBox_ChooseFile.currentData()

        if not file_path:
            raise ValueError("No ModelFit file selected.")

        data = self.parent.ModelFit_dictionary.get(file_path)

        if data is None:
            raise ValueError("No ModelFit data available.")

        if "Fit" not in data:
            raise ValueError("The selected file has not been fitted.")

        root, _ = os.path.splitext(base_file_path)
        save_path = f"{root}.xlsx"

        # 1. Original Data
        original_df = pd.DataFrame({
            "Time": data["Time"],
            "Amplitude": data["Signal"],
        })

        # 2. Fitting Data
        fitting_data = {
            "Time": data["Time"],
            "Cumulative": data["Cumulative"],
        }

        for row in range(1, 5):
            key = f"Contribution{row}"

            if key in data:
                fitting_data[f"Contribution {row}"] = data[key]

        fitting_df = pd.DataFrame(fitting_data)

        # 3. Fitting Metadata
        metadata_df = self._fit_table_dataframe()

        try:
            with pd.ExcelWriter(save_path, engine="openpyxl") as writer:
                original_df.to_excel(
                    writer,
                    sheet_name="Original Data",
                    index=False,
                )

                fitting_df.to_excel(
                    writer,
                    sheet_name="Fitting Data",
                    index=False,
                )

                metadata_df.to_excel(
                    writer,
                    sheet_name="Fitting Metadata",
                    index=False,
                )

        except Exception:
            logger.exception("ModelFit Excel export failed")
            raise

        return save_path

    def _fit_table_dataframe(self):
        table = self.ui.ModelFit_Table_Data

        headers = []

        for column in range(table.columnCount()):
            header = table.horizontalHeaderItem(column)
            headers.append(header.text() if header else f"Column {column + 1}")

        rows = []

        for row in range(table.rowCount()):
            row_data = []

            for column in range(table.columnCount()):
                item = table.item(row, column)
                row_data.append(item.text() if item else "")

            rows.append(row_data)

        return pd.DataFrame(rows, columns=headers)
