class DQTabController:
    def __init__(self, main_window):
        self.main_window = main_window

    def update_graphs(self):
        if len(self.main_window.selected_files_DQ_single) > 1:
            self.main_window.linearization()
            self.main_window.plot_fit()

    def plot_t2_graph(self):
        self.main_window.dq_t2_graph()
