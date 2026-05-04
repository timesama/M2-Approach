from SE_tab import SETabUIController


class SETabController:
    def __init__(self, main_window):
        self.main_window = main_window
        self.ui_controller = SETabUIController(main_window)

    def update_graphs(self):
        self.ui_controller.update_graphs()
