# controller for fitting FID

import os

from controllers.base_tab_controller import BaseTabController


class ModelFitController(BaseTabController):

    def connect_signals(self):
        pass

    def on_files_load(self):
        print('files have been laoded')

    def on_file_chosen(self):
        # Prepcosseing as in DQ
        # read saved properties - if there was fit. Otherwise - default empty everything - no fit

        # Plot the original data
        # AND - if there was fit done - if there were things chosen for this fit, what to plot -original,cumulative,contributions.
        print('run the pre-pocessing steps, as in DQ')

    def run_fitting(self):
        print('fitting')
        # read settings: how many functions, what are the initial parameters, what parameters are fixed
        # plot accroding to the things that are included in the plot: original\cumulative\contributions
        # fill in the table: amount of rows - amount of parameters (fixed not fixed doesnt matter): columns are parameter | value | standard error | dependency + add in the end R2 and Chi reduced - only values for them
        # save settings on successful fit

        # if not successful - window with error, clear-default parameters, but leave the plot (so like on file chosen with no fit)

    def delete_file(self):
        print('TODO: feature delete file in progress')
        # clean the file from the dictionary: clear the settings of the file: clear the table: clear the modul construction goup - return everything there to default : clean the plot: chosoe -1 index on the combobox

    def save_results(self):
        print('TODO: save results feature in progress')
        # Only the chosen file data
        # excel with tabs
        # 1. original data : Time - Amplitude
        # 2. Fitting data: Time - Amplitude cumulative, amplitude of contribution(s)
        # 3. Fitting Metadata - kist simply save the table data