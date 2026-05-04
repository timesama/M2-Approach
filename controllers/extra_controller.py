class ExtraTabController:
    def __init__(self, main_window):
        self.main_window = main_window

    def hide_activation_energy(self):
        self.main_window.hide_Eact()
