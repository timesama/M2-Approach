class GSTabController:
    def __init__(self, main_window):
        self.main_window = main_window

    def calculate(self):
        self.main_window.calculate_sqrt_time()

    def plot(self):
        self.main_window.plot_sqrt_time()
