import sys
from PySide6.QtWidgets import QApplication
from main_window import MainWindow
from logging_system import install_function_logger


def main():
    install_function_logger()
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
