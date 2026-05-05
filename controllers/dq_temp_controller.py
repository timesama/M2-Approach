import os
import re

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem
from scipy.optimize import curve_fit

import Calculator as Cal
from controllers.base_tab_controller import BaseTabController


class DQTempTabController(BaseTabController):
    def update_DQ_comparison(self):
        table = self.ui.table_DQ_2
        table.setRowCount(len(self.parent.selected_DQfiles))
        for row, parent_folder in enumerate(self.parent.selected_DQfiles, start=0):
            foldername = os.path.dirname(parent_folder)
            filename = os.path.basename(parent_folder)
            try:
                data = np.loadtxt(parent_folder, delimiter=',')
                if data.shape[1] < 3:
                    QMessageBox.warning(self.parent, "Invalid Data", f"The file {foldername} does not have at least 3 columns and will be removed from the table and file list.", QMessageBox.Ok)
                    table.removeRow(row)
                    del self.parent.selected_DQfiles[row]
            except Exception as e:
                QMessageBox.warning(self.parent, "Invalid Data", f"The file {foldername} {e} Removed from the table and file list.", QMessageBox.Ok)
                table.removeRow(row)
                del self.parent.selected_DQfiles[row]
            item = QTableWidgetItem(filename)
            item_name = QTableWidgetItem()
            pattern = r'Table_DQ_(-?[0-9]+).*.csv'
            try:
                item_xaxis = float(re.search(pattern, filename).group(1))
                item_name.setData(Qt.EditRole, item_xaxis)
            except Exception:
                item_name.setData(Qt.EditRole, float(row + 1))
            table.setItem(row, 0, item)
            table.setItem(row, 1, item_name)
        table.resizeColumnsToContents()
        self.launch()

    def launch(self):
        try:
            self.parent.dq_t2 = {}
            for row, parent_folder in enumerate(self.parent.selected_DQfiles, start=0):
                data = np.loadtxt(parent_folder, delimiter=',')
                self.parent.dq_t2[row] = data[:, [4, 5]]
            self.update_DQ_comparison_plot()
        except Exception as e:
            QMessageBox.warning(self.parent, "Corrupted file", f"Couldn't analyse the {os.path.dirname(parent_folder)} because {e}", QMessageBox.Ok)

    def update_DQ_comparison_plot(self):
        cmap = pg.ColorMap([0, len(self.parent.dq_t2)], [pg.mkColor('b'), pg.mkColor('r')])
        self.parent.dq_comparison_distribution = {'File name': [], 'X axis': [], 'Center': [], 'FWHM': [], 'Lorentz ratio': [], 'Fitting type': [], 'T2 limit': []}
        legend = self.ui.DQ_Widget_4.addLegend(offset=(0, 0))
        legend_fc = self.ui.DQ_Widget_5.addLegend()
        self.ui.DQ_Widget_5.clear(); self.ui.DQ_Widget_4.clear(); self.ui.DQ_Widget_polyFit.clear()
        if legend is not None: legend.clear(); legend.setPen((0, 0, 0))
        if legend_fc is not None: legend_fc.clear(); legend_fc.setPen((0, 0, 0))
        legend.anchor(itemPos=(1, 0), parentPos=(1, 0))
        center_g, center_l, center_v, center_d, comparison_par = [], [], [], [], []
        for row, (key, data) in zip(range(self.ui.table_DQ_2.rowCount()), self.parent.dq_t2.items()):
            file_name_item = self.ui.table_DQ_2.item(row, 1)
            file_item = self.ui.table_DQ_2.item(row, 0)
            if file_name_item is not None:
                file_name = file_name_item.text()
                if file_name != 'hide':
                    try: _cp = float(file_name)
                    except Exception:
                        _cp = 0; self.ui.table_DQ_2.setItem(row, 1, QTableWidgetItem('0'))
                    comparison_par.append(_cp)
                else:
                    continue
            t2_lin, dq_norm = data[:, 0], data[:, 1]
            p = [1, 5, 5, 0]; b = ([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1, np.inf]); b1 = ([0, 0, 0, -10], [np.inf, np.inf, np.inf, np.inf])
            t2_fit = np.arange(0, np.max(t2_lin) + 0.001, 0.01)
            params, _ = curve_fit(Cal.gaussian, t2_lin, dq_norm, p0=p, bounds=b1); cen_g, fwhm_g = params[1], params[2]; center_g.append(cen_g)
            params, _ = curve_fit(Cal.lorenz, t2_lin, dq_norm, p0=p, bounds=b1); cen_l, fwhm_l = params[1], params[2]; center_l.append(cen_l)
            params, _ = curve_fit(Cal.voigt, t2_lin, dq_norm, bounds=b); cen_v, fwhm_v = params[1], params[2]; center_v.append(cen_v)
            y_fit = Cal.voigt(t2_fit, *params)
            t2d, ndqd, center_derivative = Cal.derivative_peak_find(t2_lin, dq_norm); center_d.append(center_derivative)
            color = tuple(cmap.map(key))
            self.ui.DQ_Widget_4.plot(t2_lin, dq_norm, pen=None, symbolPen=None, symbol='o', symbolBrush=color, symbolSize=10, name=file_name)
            self.ui.DQ_Widget_4.plot(t2_fit, y_fit, pen=color)
            self.ui.DQ_Widget_polyFit.plot(t2_lin, dq_norm, pen=None, symbolPen=None, symbol='o', symbolBrush=color, symbolSize=10)
            self.ui.DQ_Widget_polyFit.plot(t2d, ndqd, pen=color)
            def set_num(r,c,v):
                it=QTableWidgetItem(); it.setData(Qt.EditRole, round(float(v),2)); self.ui.table_DQ_2.setItem(r,c,it)
            set_num(row,2,cen_g); set_num(row,3,cen_l); set_num(row,4,cen_v); set_num(row,5,center_derivative)
            set_num(row,6,fwhm_g); set_num(row,7,fwhm_l); set_num(row,8,fwhm_v)
        self.ui.DQ_Widget_5.plot(comparison_par, center_g, pen='r', symbolPen=None, symbol='o', symbolBrush='r', name='Gaus')
        self.ui.DQ_Widget_5.plot(comparison_par, center_l, pen='b', symbolPen=None, symbol='o', symbolBrush='b', name='Lorenz')
        self.ui.DQ_Widget_5.plot(comparison_par, center_v, pen='k', symbolPen=None, symbol='o', symbolBrush='k', name='Voigt')
        self.ui.DQ_Widget_5.plot(comparison_par, center_d, pen='g', symbolPen=None, symbol='o', symbolBrush='g', name='Derivative')
