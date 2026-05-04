from PySide6.QtWidgets import QDialog, QTableWidgetItem
from ui_GroupForm import Ui_GroupForm


class GroupWindow(QDialog):
    def __init__(self, parent=None, data=None, array=None):
        super().__init__(parent)
        self.ui = Ui_GroupForm()
        self.ui.setupUi(self)

        self.ui.pushButton_group.clicked.connect(self.create_group)
        self.ui.pushButton_clear.clicked.connect(self.clear)
        self.ui.pushButton_cancel.clicked.connect(self.close)
        self.ui.pushButton_done.clicked.connect(self.done_clicked)

        self.group_counter = 1
        self.group_dict = {}

    def done_clicked(self):
        print("Done clicked. Sending data back...")
        self.accept()  # This will close the dialog and return QDialog.Accepted

    def copy_table_data(self, source_table):
        rows = source_table.rowCount()
        cols = source_table.columnCount()

        self.ui.tableWidget.setRowCount(rows)
        self.ui.tableWidget.setColumnCount(cols)

        for col in range(cols):
            header = source_table.horizontalHeaderItem(0)
            if header:
                self.ui.tableWidget.setHorizontalHeaderItem(0, QTableWidgetItem(header.text()))

        for row in range(rows):
            for col in range(cols):
                item = source_table.item(row, col)
                if item:
                    self.ui.tableWidget.setItem(row, col, QTableWidgetItem(item.text()))

        current_cols = self.ui.tableWidget.columnCount()
        self.ui.tableWidget.setColumnCount(current_cols + 1)
        self.ui.tableWidget.setHorizontalHeaderItem(current_cols, QTableWidgetItem("Group"))

        self.ui.tableWidget.setSortingEnabled(True)

    def create_group(self):
        selected_rows = set()
        for item in self.ui.tableWidget.selectedItems():
            selected_rows.add(item.row())
        if not selected_rows:
            print("No rows selected.")
            return
        selected_rows = sorted(selected_rows)

        group_data = []
        for row in selected_rows:
            row_data = []
            for col in range(self.ui.tableWidget.columnCount() - 1):
                item = self.ui.tableWidget.item(row, col)
                row_data.append(item.text() if item else "")
            group_data.append(row_data)
            self.ui.tableWidget.setItem(row, self.ui.tableWidget.columnCount() - 1, QTableWidgetItem(str(self.group_counter)))

        self.group_dict[self.group_counter] = group_data
        self.group_counter += 1

    def clear(self):
        self.group_counter = 1
        self.group_dict = {}

        rows = self.ui.tableWidget.rowCount()
        cols = self.ui.tableWidget.columnCount()

        for row in range(rows):
            self.ui.tableWidget.setItem(row, cols - 1, QTableWidgetItem(""))
