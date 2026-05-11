import logging

logger = logging.getLogger(__name__)


class BaseTabController:
    def __init__(self, ui, state, parent=None):
        self.ui = ui
        self.state = state
        self.parent = parent

    def show_status(self, message, timeout=8000):
        """Show a concise user-facing footer/status-bar message and log it."""
        if self.parent is not None and hasattr(self.parent, "show_status"):
            self.parent.show_status(message, timeout=timeout)
            return
        if self.parent is not None and hasattr(self.parent, "statusBar"):
            self.parent.statusBar().showMessage(message, timeout)
        logger.info(message)

    def _status(self, message, timeout=8000):
        """Backward-compatible alias for controllers that used the RecFID helper."""
        self.show_status(message, timeout=timeout)

    def read_column_values(self, table, column):
        values = []
        for row in range(table.rowCount()):
            item = table.item(row, column)
            if item is None:
                continue
            try:
                values.append(float(item.text()))
            except ValueError:
                continue
        return values
