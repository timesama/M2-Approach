import json
import os
import sys
import winreg
from PySide6.QtWidgets import QFileDialog, QMessageBox


class SaveFilesDialog(QFileDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptMode(QFileDialog.AcceptSave)
        self.last_saved_file_path = None

    def save_data_as_csv(self, directory, table, files, default_filename):
        try:
            key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, r"Software\MyApp", 0, winreg.KEY_READ)
            directory, _ = winreg.QueryValueEx(key, "SelectedFolder")
            winreg.CloseKey(key)
        except Exception:
            directory = os.path.dirname(sys.argv[0])

        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File As", directory + '/' + default_filename, "CSV files (*.csv)", options=options)
        if file_path:
            self.last_saved_file_path = file_path
            with open(file_path, 'w') as f:
                for row in range(table.rowCount()):
                    row_values = []
                    for col in range(table.columnCount()):
                        item = table.item(row, col)
                        row_values.append(item.text() if item is not None else "")
                    f.write(','.join(row_values) + '\n')
            files_list_path = os.path.splitext(file_path)[0] + '_files.json'
            with open(files_list_path, 'w') as file_list:
                json.dump(files, file_list)
