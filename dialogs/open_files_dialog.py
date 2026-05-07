import logging
import os
import sys
import winreg

from PySide6.QtWidgets import QFileDialog

logger = logging.getLogger(__name__)

State_multiple_files = None


class OpenFilesDialog(QFileDialog):
    def __init__(self, parent=None):
        global State_multiple_files
        super().__init__(parent)

        if State_multiple_files:
            self.setFileMode(QFileDialog.ExistingFiles)
        else:
            self.setFileMode(QFileDialog.ExistingFile)

        self.setNameFilter(str("Data (*.dat *.txt *.csv)"))

        directory = self.get_initial_directory()
        self.setDirectory(directory)

        self.selected_files = []

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
            logger.warning("Could not read the initial directory: %s", e)
            return os.path.dirname(sys.argv[0])
