from PySide6.QtCore import QEvent
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import QApplication, QTableWidget, QWidget


class TableCopyEnabler(QWidget):
    def __init__(self, root_widget):
        super().__init__()
        self.enable_table_copying(root_widget)

    def enable_table_copying(self, widget):
        """Install copy handling on every table below the root widget."""
        for table in widget.findChildren(QTableWidget):
            table.setSelectionBehavior(QTableWidget.SelectItems)
            table.setSelectionMode(QTableWidget.ExtendedSelection)
            table.installEventFilter(self)

    def eventFilter(self, obj, event):
        if isinstance(obj, QTableWidget) and event.type() == QEvent.KeyPress:
            if event.matches(QKeySequence.Copy):
                self.copy_table_selection(obj)
                return True
        return super().eventFilter(obj, event)

    def copy_table_selection(self, table):
        """Copy the selected table ranges to the clipboard as tab-separated text."""
        selected_ranges = table.selectedRanges()
        if not selected_ranges:
            return

        copied_text = ""
        for r in selected_ranges:
            for row in range(r.topRow(), r.bottomRow() + 1):
                row_data = []
                for col in range(r.leftColumn(), r.rightColumn() + 1):
                    item = table.item(row, col)
                    row_data.append(item.text() if item else "")
                copied_text += "\t".join(row_data) + "\n"

        QApplication.clipboard().setText(copied_text.strip())
