from PySide6.QtCore import Signal
from PySide6.QtWidgets import QDialog
from ui_Notification import Ui_Note


class NotificationDialog(QDialog, Ui_Note):
    stateChanged = Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.buttonBox.accepted.connect(self.accepted_handler)
        self.buttonBox.rejected.connect(self.rejected_handler)
        self.flag_ok = False

    def accepted_handler(self):
        self.flag_ok = True
        self.stateChanged.emit(self.flag_ok)
        self.accept()

    def rejected_handler(self):
        self.flag_ok = False
        self.stateChanged.emit(self.flag_ok)
        self.reject()
