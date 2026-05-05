import sys
from PySide6.QtWidgets import QApplication
from main_window import MainWindow
from utils.logging_config import setup_logging


def main():
    setup_logging()
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
