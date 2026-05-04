class BaseTabController:
    def __init__(self, ui, state, parent=None):
        self.ui = ui
        self.state = state
        self.parent = parent

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
