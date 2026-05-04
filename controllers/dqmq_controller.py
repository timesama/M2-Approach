class DQMQTabController:
    def __init__(self, main_window):
        self.main_window = main_window

    def plot_original(self):
        self.main_window.plot_original()

    def plot_norm(self):
        self.main_window.plot_norm()

    def plot_diff(self):
        self.main_window.plot_diff()

    def plot_nDQ(self):
        self.main_window.plot_nDQ()
