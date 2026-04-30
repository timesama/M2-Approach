"""SE tab UI responsibilities.

This module keeps plotting and UI-only behavior for the SE tab separate
from the main window orchestration.
"""


class SETabUIController:
    def __init__(self, main_window):
        self.main_window = main_window

    def update_graphs(self):
        mw = self.main_window

        x = mw.read_column_values(mw.ui.table_SE, 0)
        text = mw.ui.comboBox_SE_chooseY.currentText()

        if text == "SC":
            y = mw.read_column_values(mw.ui.table_SE, 1)
            mw.ui.SEWidget.getAxis('left').setLabel("SC")

        elif text == "M₂":
            y = mw.read_column_values(mw.ui.table_SE, 2)
            mw.ui.SEWidget.getAxis('left').setLabel("M₂")

        elif text == "T₂*":
            y = mw.read_column_values(mw.ui.table_SE, 3)
            mw.ui.SEWidget.getAxis('left').setLabel("T₂*")

        else:  # Set y
            y = [0]
            x = [0]
            mw.ui.SEWidget.getAxis('left').setLabel("Not Set")
            return

        mw.ui.SEWidget.clear()
        mw.ui.SEWidget.plot(
            x, y, pen=None,
            symbol='o', symbolPen=None,
            symbolBrush=(255, 0, 0, 255),
            symbolSize=10
        )

        if mw.group_data_SE:
            y_index = mw.ui.comboBox_SE_chooseY.currentIndex()
            y_col = y_index + 1

            for i, (_group_number, group_rows) in enumerate(mw.group_data_SE.items()):
                group_x = []
                group_y = []

                for row_data in group_rows:
                    try:
                        x_val = float(row_data[0])
                        y_val = float(row_data[y_col])
                        group_x.append(x_val)
                        group_y.append(y_val)
                    except (ValueError, IndexError):
                        continue

                sorted_points = sorted(zip(group_x, group_y), key=lambda p: p[0])
                if len(sorted_points) > 1:
                    xs, ys = zip(*sorted_points)

                    color = mw.tab10_colors[i % len(mw.tab10_colors)]

                    mw.ui.SEWidget.plot(
                        xs, ys,
                        pen={'color': color, 'width': 2},
                        symbol='o',
                        symbolBrush=color,
                        symbolPen=None,
                        symbolSize=8
                    )
