import json
import os
import sys
import winreg
from PySide6.QtWidgets import QFileDialog, QMessageBox


class SaveFilesDialog(QFileDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptMode(QFileDialog.AcceptSave)

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

    def save_file_in_sef(self, wtf, dictionary, tau, n, begin, save):
        if not save:
            return
        try:
            key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, r"Software\MyApp", 0, winreg.KEY_READ)
            directory, _ = winreg.QueryValueEx(key, "SelectedFolder")
            winreg.CloseKey(key)
        except Exception:
            directory = os.path.dirname(sys.argv[0])
        try:
            options = QFileDialog.Options()
            default_add = '_recalculated_tau_' + str(n)
            file_path, _ = QFileDialog.getSaveFileName(self, "Save File As", directory + '/' + begin + default_add, "SEF files (*.sef)", options=options)
            with open(file_path, 'w') as f:
                f.write('STELAR Export File\n\n')
                f.write('_BRLX______\t_T1________\t_R1________\t____%err___\t+-err(R1)__\tZone\tFile\n\n')
                for key in dictionary:
                    Frequency = dictionary[key]['X Axis']
                    T1 = dictionary[key][tau]
                    Omega_1 = 0 if T1 == 0 else 1000 / T1
                    f.write(f'{Frequency}\t{T1}\t{Omega_1}\n')
        except Exception:
            QMessageBox.warning(self, "Save failed", "Sorry, couldn't save the separate file in .sef for magnetization.", QMessageBox.Ok)
