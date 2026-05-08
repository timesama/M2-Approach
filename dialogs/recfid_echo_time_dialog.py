"""Manual echo-time assignment dialog for RecFID SE files."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QHBoxLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout


class RecFIDEchoTimeDialog(QDialog):
    """Two-column filename/echo-time editor used when SE echo times cannot be parsed."""

    def __init__(self, files, echo_times=None, parent=None):
        super().__init__(parent)
        self.setObjectName("RecFIDManual_Dialog_EchoTimes")
        self.setWindowTitle("Manual Echo Time")
        self.resize(700, 400)

        echo_times = echo_times or []
        self.RecFIDManual_Table_EchoTimes = QTableWidget(len(files), 2, self)
        self.RecFIDManual_Table_EchoTimes.setObjectName("RecFIDManual_Table_EchoTimes")
        self.RecFIDManual_Table_EchoTimes.setHorizontalHeaderLabels(["Filename", "Echo Time"])

        for row, file_path in enumerate(files):
            file_item = QTableWidgetItem(Path(file_path).name)
            file_item.setData(Qt.UserRole, str(file_path))
            file_item.setFlags(file_item.flags() & ~Qt.ItemIsEditable)
            self.RecFIDManual_Table_EchoTimes.setItem(row, 0, file_item)
            value = echo_times[row] if row < len(echo_times) else ""
            self.RecFIDManual_Table_EchoTimes.setItem(row, 1, QTableWidgetItem(str(value)))
        self.RecFIDManual_Table_EchoTimes.resizeColumnsToContents()

        self.RecFIDManual_Button_Accept = QPushButton("Accept")
        self.RecFIDManual_Button_Accept.setObjectName("RecFIDManual_Button_Accept")
        self.RecFIDManual_Button_Accept.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)

        button_layout = QHBoxLayout()
        button_layout.addStretch(1)
        button_layout.addWidget(self.RecFIDManual_Button_Accept)
        button_layout.addWidget(cancel_button)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Enter echo times for each loaded SE file."))
        layout.addWidget(self.RecFIDManual_Table_EchoTimes)
        layout.addLayout(button_layout)

    def echo_times(self):
        """Return manually entered echo times as floats."""
        values = []
        for row in range(self.RecFIDManual_Table_EchoTimes.rowCount()):
            item = self.RecFIDManual_Table_EchoTimes.item(row, 1)
            if item is None or item.text().strip() == "":
                raise ValueError("Every SE file needs an echo time.")
            values.append(float(item.text()))
        return values
