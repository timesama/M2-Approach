import os
import sys
import winreg
from PySide6.QtWidgets import QFileDialog

State_multiple_files = None


class OpenFilesDialog(QFileDialog):
    def __init__(self, parent=None):
        global State_multiple_files
        super().__init__(parent)

        if State_multiple_files:
            self.setFileMode(QFileDialog.ExistingFiles)  # Allow selecting multiple files
        else:
            pass

        self.setNameFilter(str("Data (*.dat *.txt *.csv *.sef)"))

        directory = self.get_initial_directory()
        self.setDirectory(directory)

        self.selected_files = []  # dictionary to store selected file paths

    def get_initial_directory(self):
        """Retrieve the initial directory from the registry."""
        try:
            key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, r"Software\MyApp", 0, winreg.KEY_READ)
            directory, _ = winreg.QueryValueEx(key, "SelectedFolder")
            winreg.CloseKey(key)
            return directory
        except FileNotFoundError:
            return os.path.dirname(sys.argv[0])
        except Exception as e:
            print(f"Couldn't read the initial directory: {e}")
            return os.path.dirname(sys.argv[0])

    def on_file_selected(self):
        options = QFileDialog.Options()
        files, _ = QFileDialog.getOpenFileNames(self, "Load Files", "", "Data Files (*.dat *.txt *.csv)", options=options)
        if files:
            self.selected_files.extend(files)
