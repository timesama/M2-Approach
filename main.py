import sys
from PySide6.QtWidgets import QApplication
from mainwindow import MainWindow
from utils.logging_config import setup_logging

from PySide6.QtGui import QPalette, QColor, Qt
from PySide6.QtCore import QSettings
import pyqtgraph as pg


def main():
    setup_logging()
    app = QApplication(sys.argv)

    app.setStyle("Fusion")
    palette = app.palette()
    pg.setConfigOption("foreground", palette.color(QPalette.Text))
    if app.styleHints().colorScheme() == Qt.ColorScheme.Light:
        palette.setColor(QPalette.Window, QColor("white"))
        palette.setColor(QPalette.Base, QColor("white"))
        app.setPalette(palette)

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()


