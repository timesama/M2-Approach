import sys
from PySide6.QtWidgets import QApplication
from mainwindow import MainWindow
from utils.logging_config import setup_logging

from PySide6.QtGui import QPalette, QColor, Qt
from PySide6.QtCore import QSettings


def main():
    setup_logging()
    app = QApplication(sys.argv)

    app.setStyle("Fusion")

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()


