import os,re
import numpy as np
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem
from scipy.signal import savgol_filter
import Calculator as Cal
from controllers.base_tab_controller import BaseTabController
from dialogs.open_files_dialog import OpenFilesDialog
from dialogs.phasing_manual import PhasingManual, Frequency, Re_spectra, Im_spectra


class GeneralSEDQController(BaseTabController):
    def __init__(self, ui, state, parent=None, se_controller=None, dq_controller=None):
        super().__init__(ui, state, parent)
        self.se_controller = se_controller
        self.dq_controller = dq_controller

    def analysis(self):
        mw=self.parent
        while self.ui.comboBox_4.count()>0: self.ui.comboBox_4.removeItem(0)
        if mw.tab=='SE': files=mw.selected_files; self.ui.SEWidget.clear(); self.ui.comboBox_SE_chooseY.setCurrentIndex(-1)
        else: files=mw.selected_files_DQ_single; self.ui.DQ_Widget_1.clear(); self.ui.DQ_Widget_2.clear(); self.ui.DQ_Widget_4.clear(); self.ui.textEdit_4.setText(''); self.ui.comboBox_FunctionDQ.setCurrentIndex(-1)
        if len(files)==0: return
        if self.ui.checkBox_glycerol.isChecked() and mw.selected_files_gly==[]:
            self.open_select_dialog_glycerol()
        if self.ui.checkBox_baseline.isChecked() and mw.selected_files_empty==[]:
            self.open_select_dialog_baseline()
        mw.disable_buttons(); self.ui.btn_SelectFiles.setEnabled(False); self.ui.btn_Load.setEnabled(False); self.ui.radioButton.setEnabled(False); self.ui.comboBox_4.setCurrentIndex(-1)
        mw.window_array = np.linspace(self.ui.SmoothWindowFrom.value(), self.ui.SmoothWindowTo.value(), len(files), dtype=np.int32) if self.ui.checkBox_Smooth.isChecked() else np.array([])
        for i,file_path in enumerate(files,start=1):
            try:
                f_g = mw.selected_files_gly[i-1] if self.ui.checkBox_glycerol.isChecked() else []
                f_e = mw.selected_files_empty[i-1] if self.ui.checkBox_baseline.isChecked() else []
                self.process_file_data(file_path,f_g,f_e,i)
            except Exception:
                self.analysis_error(file_path, files)
        self.update_legends_and_dq_graphs(); self.ui.btn_Start.setStyleSheet('background-color: none')

    def analysis_error(self, file_path, files):
        mw=self.parent
        QMessageBox.warning(mw, "Invalid Data", f"Couldn't read the {os.path.basename(file_path)}. Deleting the file.", QMessageBox.Ok)
        files.remove(file_path); self.ui.btn_SelectFiles.setEnabled(True); self.ui.btn_Start.setStyleSheet("background-color: none")
        self.ui.FidWidget.clear(); self.ui.FFTWidget.clear(); self.ui.btn_Phasing.setEnabled(False); mw.enable_buttons()

    def update_legends_and_dq_graphs(self):
        self.parent.enable_buttons()
        if self.parent.tab=='DQ': self.dq_controller.update_graphs()

    def process_file_data(self,file_path,file_path_gly,file_path_empty,i):
        global Frequency, Re_spectra, Im_spectra
        mw=self.parent; filename=os.path.basename(file_path)
        subtract=self.ui.checkBox_baseline.isChecked()
        data=np.loadtxt(file_path); x,y,z=data[:,0],data[:,1],data[:,2]
        Time,Re,Im=Cal.analysis_time_domain(file_path,file_path_empty,subtract)
        if self.ui.checkBox_glycerol.isChecked():
            Time_r,Re_r,Im_r=Cal.analysis_time_domain(file_path_gly,[],False); Time,Re,Im=Cal.magnet_inhomogenity_correction(Time,Time_r,Re,Re_r,Im,Im_r)
        if self.ui.checkBox_long_component.isChecked(): Time,Re,Im=Cal.subtract_long_component(Time,Re,Im)
        Amp=Cal._calculate_amplitude(Re,Im); mw.update_graphs(Time,Amp,Re,Im,self.ui.FidWidget)
        Time_fid,Fid=Cal.final_analysis_time_domain(Time,Re,Im,2**16); Frequency=Cal._calculate_frequency_scale(Time_fid)
        FFT=np.fft.fftshift(np.fft.fft(Fid))
        if mw.window_array.size!=0:
            window=mw.window_array[i-1]; FFT=np.array(savgol_filter(np.real(FFT),window,1)+1j*savgol_filter(np.imag(FFT),window,1))
        Amp_spectra,Re_spectra,Im_spectra=Cal._simple_baseline_correction(FFT); Real_apod=Cal._calculate_apodization(Re_spectra,Frequency)
        mw.update_graphs(Frequency,Amp_spectra,Re_spectra,Im_spectra,self.ui.FFTWidget)
        if self.ui.comboBox_4.findText(filename)==-1: self.ui.comboBox_4.addItem(filename)
        if self.ui.comboBox_4.currentIndex()==-1:
            M2,T2=Cal._calculate_M2(Real_apod,Frequency)
            if mw.tab=='SE':
                m=re.search(r'.*_(-?\s*\d+\.?\d*).*.dat',filename); temperature=m.group(1) if m else '0'
                times=[int(self.ui.SC_short.value()*2),int(self.ui.SC_short_2.value()*2),int(self.ui.SC_long.value()*2),int(self.ui.SC_long_2.value()*2)]
                SFC=Cal.calculate_SC(Amp,times,self.ui.radioButton_absolute.isChecked())
                mw.ui.table_SE.setRowCount(i); mw.fill_table(mw.ui.table_SE,temperature,SFC,M2,T2,i); mw.ui.table_SE.setItem(i-1,4,QTableWidgetItem(filename))
            elif mw.tab=='DQ':
                m=re.search(r'_(\d+\.\d+)_',filename); dq_time=m.group(1) if m else '0'; DQ=Cal.calculate_DQ_intensity(x,Cal._calculate_amplitude(y,z))
                mw.ui.table_DQ.setRowCount(i); mw.fill_table(mw.ui.table_DQ,dq_time,DQ,M2,T2,i)

    def after_phasing(self):
        global Frequency, Re_spectra, Im_spectra
        mw=self.parent; i=self.ui.comboBox_4.currentIndex(); Real_apod=Cal._calculate_apodization(Re_spectra,Frequency); Amp_spectra=Cal._calculate_amplitude(Re_spectra,Im_spectra)
        mw.update_graphs(Frequency,Amp_spectra,Re_spectra,Im_spectra,self.ui.FFTWidget); M2,T2=Cal._calculate_M2(Real_apod,Frequency)
        table=self.ui.table_SE if mw.tab=='SE' else self.ui.table_DQ
        table.setItem(i,2,QTableWidgetItem(str(round(M2,6)))); table.setItem(i,3,QTableWidgetItem(str(round(T2,6))))
        if mw.tab=='SE': self.se_controller.update_graphs()
        elif mw.tab=='DQ': self.dq_controller.update_graphs()

    def open_phasing_manual(self):
        mw=self.parent; mw.phasing_manual_window=PhasingManual(); mw.phasing_manual_window.read_data(); mw.phasing_manual_window.show(); mw.phasing_manual_window.closed.connect(self.after_phasing)

    def open_select_dialog_glycerol(self):
        dlg=OpenFilesDialog(self.parent); dlg.setWindowTitle('Select Glycerol Files')
        if dlg.exec(): self.parent.selected_files_gly.extend(dlg.selectedFiles())
        self.ui.btn_Start.setEnabled(True); self.ui.btn_Add.setEnabled(True)

    def open_select_dialog_baseline(self):
        dlg=OpenFilesDialog(self.parent); dlg.setWindowTitle('Select Baseline Files')
        if dlg.exec(): self.parent.selected_files_empty.extend(dlg.selectedFiles())
        self.ui.btn_Start.setEnabled(True); self.ui.btn_Add.setEnabled(True)
