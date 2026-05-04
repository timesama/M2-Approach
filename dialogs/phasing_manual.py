import numpy as np
from scipy.signal import savgol_filter
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QDialog

import Calculator as Cal
from ui_PhasingManual import Ui_Phasing as Ui_PhasingManual

Frequency = None
Re_spectra = None
Im_spectra = None


class PhasingManual(QDialog):
    closed = Signal()
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_PhasingManual()
        self.ui.setupUi(self)
        g = self.ui.PhasingGraph
        g.getAxis('bottom').setLabel("Frequency, MHz")
        g.getAxis('left').setLabel("Amplitude, a.u.")
        g.setTitle("Phasing")
        g.addLegend()
        self.ui.pushButton_2.clicked.connect(self.zero)
        self.ui.verticalSlider_a.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_b.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_c.valueChanged.connect(self.value_changed)
        self.ui.verticalSlider_d.valueChanged.connect(self.value_changed)
        self.ui.dial.valueChanged.connect(self.smoothing_changed)
        self.ui.pushButton_3.clicked.connect(self.manual_read)
        self.ui.pushButton.clicked.connect(self.save_data)
        self.Real_freq_phased = None

    def read_data(self): self.zero()
    def save_data(self): self.close()
    def closeEvent(self, event): self.closed.emit(); super().closeEvent(event)
    def zero(self):
        self.a=self.b=self.c=self.d=self.Smooth=0
        for s,b in [(self.ui.verticalSlider_a,self.ui.Box_a),(self.ui.verticalSlider_b,self.ui.Box_b),(self.ui.verticalSlider_c,self.ui.Box_c),(self.ui.verticalSlider_d,self.ui.Box_d)]: s.setValue(0); b.setValue(0)
        self.ui.dial.setValue(0); self.ui.Box_smooth.setValue(0)
        if self.Real_freq_phased is not None: self.process_data()
    def process_data(self): self.Real_freq_phased=self.calculate_phase(); self.update_plot(); self.update_text()
    def calculate_phase(self):
        phi=self.a+self.b*Frequency+self.c*Frequency**2+self.d*Frequency**3
        y=Re_spectra*np.cos(np.deg2rad(phi))-Im_spectra*np.sin(np.deg2rad(phi))
        if self.Smooth>1: y=savgol_filter(y, window_length=int(self.Smooth), polyorder=1)
        return y
    def update_plot(self):
        self.ui.PhasingGraph.clear(); self.ui.PhasingGraph.plot(Frequency, Re_spectra, pen='r', name='Original'); self.ui.PhasingGraph.plot(Frequency, self.Real_freq_phased, pen='b', name='Phased')
    def update_text(self):
        Integral=np.trapz(self.Real_freq_phased); delta=np.mean(self.Real_freq_phased[:100])-np.mean(self.Real_freq_phased[-100:])
        self.ui.Integral.setText(f"Integral: {round(Integral,3)}"); self.ui.Delta.setText(f"Delta: {round(delta,7)}")
        M2,T2=Cal._calculate_M2(Cal._calculate_apodization(self.Real_freq_phased, Frequency), Frequency)
        self.ui.M2.setText(f"M₂: {round(M2,5)}"); self.ui.T2.setText(f"T₂*: {round(T2,3)}")
    def value_changed(self):
        self.a=self.ui.verticalSlider_a.value(); self.b=self.ui.verticalSlider_b.value(); self.c=self.ui.verticalSlider_c.value(); self.d=self.ui.verticalSlider_d.value()
        self.ui.Box_a.setValue(self.a); self.ui.Box_b.setValue(self.b); self.ui.Box_c.setValue(self.c); self.ui.Box_d.setValue(self.d)
        if Frequency is not None and Re_spectra is not None and Im_spectra is not None: self.process_data()
    def smoothing_changed(self): self.Smooth=self.ui.dial.value(); self.ui.Box_smooth.setValue(self.Smooth); self.process_data()
    def manual_read(self):
        self.a=self.ui.Box_a.value(); self.b=self.ui.Box_b.value(); self.c=self.ui.Box_c.value(); self.d=self.ui.Box_d.value(); self.Smooth=self.ui.Box_smooth.value(); self.process_data()
