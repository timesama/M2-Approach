class T1T2TabController:
    def __init__(self, main_window):
        self.main_window = main_window

    def calculate(self):
        self.main_window.calculate_relaxation_time()

    def plot(self):
        self.main_window.plot_relaxation_time()
