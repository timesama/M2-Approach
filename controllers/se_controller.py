from controllers.base_tab_controller import BaseTabController


class SETabController(BaseTabController):
    def update_graphs(self):
        x = self.read_column_values(self.ui.table_SE, 0)
        text = self.ui.comboBox_SE_chooseY.currentText()

        if text == "SC":
            y = self.read_column_values(self.ui.table_SE, 1)
            self.ui.SEWidget.getAxis('left').setLabel("SC")
        elif text == "M₂":
            y = self.read_column_values(self.ui.table_SE, 2)
            self.ui.SEWidget.getAxis('left').setLabel("M₂")
        elif text == "T₂*":
            y = self.read_column_values(self.ui.table_SE, 3)
            self.ui.SEWidget.getAxis('left').setLabel("T₂*")
        else:
            self.ui.SEWidget.getAxis('left').setLabel("Not Set")
            return

        self.ui.SEWidget.clear()
        self.ui.SEWidget.plot(
            x, y, pen=None,
            symbol='o', symbolPen=None,
            symbolBrush=(255, 0, 0, 255), symbolSize=10
        )

        group_data = getattr(self.parent, "group_data_SE", {})
        colors = getattr(self.parent, "tab10_colors", [])
        if group_data:
            y_col = self.ui.comboBox_SE_chooseY.currentIndex() + 1
            for i, (_group_number, group_rows) in enumerate(group_data.items()):
                group_x, group_y = [], []
                for row_data in group_rows:
                    try:
                        group_x.append(float(row_data[0]))
                        group_y.append(float(row_data[y_col]))
                    except (ValueError, IndexError):
                        continue
                sorted_points = sorted(zip(group_x, group_y), key=lambda p: p[0])
                if len(sorted_points) > 1:
                    xs, ys = zip(*sorted_points)
                    color = colors[i % len(colors)] if colors else 'r'
                    self.ui.SEWidget.plot(
                        xs, ys,
                        pen={'color': color, 'width': 2},
                        symbol='o', symbolBrush=color, symbolPen=None, symbolSize=8
                    )
