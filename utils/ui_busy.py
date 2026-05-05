from contextlib import contextmanager
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication


@contextmanager
def busy_cursor():
    QApplication.setOverrideCursor(Qt.WaitCursor)
    QApplication.processEvents()
    try:
        yield
    finally:
        QApplication.restoreOverrideCursor()
