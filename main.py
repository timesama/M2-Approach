import sys
from PySide6.QtWidgets import QApplication
from mainwindow import MainWindow
from utils.logging_config import setup_logging

from PySide6.QtGui import QPalette, QColor, Qt
from PySide6.QtCore import QSettings


def main():
    setup_logging()
    app = QApplication(sys.argv)

    apply_theme(app)

    window = MainWindow()
    window.show()
    sys.exit(app.exec())

def is_windows_in_dark_theme():
    registry_path = r"HKEY_CURRENT_USER\Software\Microsoft\Windows\CurrentVersion\Themes\Personalize"
    settings = QSettings(registry_path, QSettings.Format.NativeFormat)

    return settings.value

def apply_theme(app: QApplication):
    app.setStyle("Fusion")
    palette= QPalette()

    if is_windows_in_dark_theme():
        palette.setColor(QPalette.ColorRole.Window, QColor(53,53,53))
    else:
        pass
        # pg.CONFIG_OPTIONS['background'] = 'w'
        # pg.CONFIG_OPTIONS['foreground'] = 'k'

if __name__ == '__main__':
    main()


