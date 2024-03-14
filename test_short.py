# This Python file uses the following encoding: utf-8
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def analysis_time_domain(file_path):
    # 1. Read data
    Time, Real, Imag = read_data(file_path)

    # 2. Crop time below zero
    T_cr, R_cr, I_cr = crop_time_zero(Time, Real, Imag)

    # 3. Phase the data
    R_ph, I_ph = time_domain_phase(R_cr, I_cr)

    # 4. Adjust Frequency
    # 4.1 Calculate Freq
    Frequency = calculate_frequency_scale(T_cr)
    # 4.2 Shift Freq
    R_sh, I_sh = adjust_frequency(Frequency, R_ph, I_ph)

    return T_cr, R_sh, I_sh

def long_component(Time_s, Time_r, Re_s, Re_r, Im_s, Im_r):
    # r stands for reference, s stands for sample
    Amp_r = calculate_amplitude(Re_s, Im_r)
    Amp_s = calculate_amplitude(Re_s, Im_s)

    # 1. Crop the arrays together (they should be of the same length, but I know, I know...)
    if len(Time_s) > len(Time_r):
        Time_  =   Time_s[:len(Time_r)]
        Amp_s   =   Amp_s[:len(Time_r)]
        Re_s    =   Re_s[:len(Time_r)]
        Im_s    =   Im_s[:len(Time_r)]
    else:  
        Time_  =   Time_r[:len(Time_s)]
        Amp_r   =   Amp_r[:len(Time_s)]
        Re_r    =   Re_r[:len(Time_s)]
        Im_r    =   Im_r[:len(Time_s)]

    # 2. Normalize (reference) components to Amplitude of the reference
    Re_gly_norm = Re_s/Amp_r
    Im_gly_norm = Im_s/Amp_r
    #Amp_gly_norm = Amp_s/Amp_r

    # 3. Cut the ranges for fitting
    minimum = find_nearest(Time_, 50)
    maximum = find_nearest(Time_, 250)

    Time_range = Time_[minimum:maximum]
    Re_gly_norm_range = Re_gly_norm[minimum:maximum]
    Im_gly_norm_range = Im_gly_norm[minimum:maximum]
    #Amp_gly_norm_range = Amp_gly_norm[minimum:maximum]

    # 7. Fit data to exponential decay
    coeff = [0.9, 400, 0.1]
    popt, pcov      = curve_fit(decaying_exponential, Time_range, Re_gly_norm_range, p0=coeff, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
    popt_, pcov_    = curve_fit(decaying_exponential, Time_range, Im_gly_norm_range, p0=[1, 400, 0])
    #popt__, pcov__ = curve_fit(decaying_exponential, Time_range, Amp_gly_norm_range,p0=coeff, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
    
    # 8. Set the ranges for subtraction
    Time_cropped = Time_[0:maximum]
    Real_cropped    = Re_gly_norm[0:maximum]
    Im_cropped      = Im_gly_norm[0:maximum]
    #Amp_cropped = Amp_gly_norm[0:maximum]

    # 9. Calculate the curves fitted to data within the desired range
    Real_fitted = decaying_exponential(Time_cropped, *popt)
    Im_fitted   = decaying_exponential(Time_cropped, *popt_)
    #Amp_fitted = decaying_exponential(Time_cropped, *popt__)

    # 10. Subtract
    Real_subtracted = Real_cropped - Real_fitted
    Im_subtracted = Im_cropped - Im_fitted

    # 11. Normalize
    Re_n, Im_n = normalize(Real_subtracted, Im_subtracted)

    return Time_cropped, Re_n, Im_n

def FFT_processing(Time, Real, Imaginary):
    # 5. Apodize the time-domain
    Re_ap, Im_ap = apodization(Time, Real, Imaginary)
    
    # 6. Add zeros
    Tim, Fid = add_zeros(Time, Re_ap, Im_ap, 16383)

    # 7. FFT
    FFT = np.fft.fftshift(np.fft.fft(Fid))
    # Frequency 
    Frequency = calculate_frequency_scale(Tim)

    # 8. Simple baseline
    Amp, Re, Im = simple_baseline_correction(FFT)

    # 9. Apodization
    Real_apod = calculate_apodization(Re, Frequency)

    # 10. M2 & T2
    M2, T2 = calculate_M2(Real_apod, Frequency)

    return M2, T2

def read_data(file_path):

    data = np.loadtxt(file_path)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]

    return x, y, z

def crop_time_zero(Time, Real, Imaginary):
    if Time[0] < 0:
        Time_start = 0
        Time_crop_idx = np.where(Time >= Time_start)[0][0]
        Time_cropped = Time[Time_crop_idx:]
        Real_cropped = Real[Time_crop_idx:]
        Imaginary_cropped = Imaginary[Time_crop_idx:]
        return Time_cropped, Real_cropped, Imaginary_cropped
    else:
        return Time, Real, Imaginary

def time_domain_phase(Real, Imaginary):
    delta = np.zeros(360)
    
    for phi in range(360):
        Re_phased = Real * np.cos(np.deg2rad(phi)) - Imaginary * np.sin(np.deg2rad(phi))
        Im_phased = Real * np.sin(np.deg2rad(phi)) + Imaginary * np.cos(np.deg2rad(phi))
        Magnitude_phased = calculate_amplitude(Re_phased, Im_phased)
        
        Re_cut = Re_phased[:50]
        Ma_cut = Magnitude_phased[:50]
        
        delta[phi] = np.mean(Ma_cut - Re_cut)
    
    idx = np.argmin(delta)
    print(idx)

    Re = Real * np.cos(np.deg2rad(idx)) - Imaginary * np.sin(np.deg2rad(idx))
    Im = Real * np.sin(np.deg2rad(idx)) + Imaginary * np.cos(np.deg2rad(idx))

    return Re, Im

def adjust_frequency(Frequency, Re, Im):
    # Create complex FID
    Fid_unshifted = np.array(Re + 1j * Im)

    # FFT
    FFT = np.fft.fftshift(np.fft.fft(Fid_unshifted))

    # Check the length of FFT and Frequency (it is always the same, this is just in case)
    if len(Frequency) != len(FFT):
        Frequency = np.linspace(Frequency[0], Frequency[-1], len(FFT))

    # Find index of max spectrum (amplitude)
    index_max = np.argmax(FFT)

    # Find index of zero (frequency)
    index_zero = find_nearest(Frequency, 0)

    # Find difference
    delta_index = index_max - index_zero

    # Shift the spectra (amplitude) by the difference in indices
    FFT_shifted = np.concatenate((FFT[delta_index:], FFT[:delta_index]))

    # iFFT
    Fid_shifted = np.fft.ifft(np.fft.fftshift(FFT_shifted))

    # Define Real, Imaginary and Amplitude
    Re_shifted = np.real(Fid_shifted)
    Im_shifted = np.imag(Fid_shifted)

    return Re_shifted, Im_shifted

def normalize(Real, Imaginary):
    Amplitude = np.sqrt(Real ** 2 + Imaginary ** 2)
    Amplitude_max = np.max(Amplitude)
    Amp = Amplitude/Amplitude_max
    Re = Real/Amplitude_max
    Im = Imaginary/Amplitude_max
    return Re, Im

def apodization(Time, Real, Imaginary):
    Amplitude = calculate_amplitude(Real, Imaginary)
    coeffs = np.polyfit(Time, Amplitude, 1)  # Fit an exponential decay function
    c = np.polyval(coeffs, Time)
    d = np.argmin(np.abs(c - 1e-5))
    sigma = Time[d]
    apodization_function = np.exp(-(Time / sigma) ** 4)
    Re_ap = Real * apodization_function
    Im_ap = Imaginary * apodization_function
    return Re_ap, Im_ap

def add_zeros(Time, Real, Imaginary, number_of_points):
    length_diff = number_of_points - len(Time)
    amount_to_add = np.zeros(length_diff+1)

    Re_zero = np.concatenate((Real, amount_to_add))
    Im_zero = np.concatenate((Imaginary, amount_to_add))

    dt = Time[1] - Time[0]
    Time_to_add = Time[-1] + np.arange(1, length_diff + 1) * dt

    Time = np.concatenate((Time, Time_to_add))
    Fid = np.array(Re_zero + 1j * Im_zero)

    return Time, Fid

def simple_baseline_correction(FFT):
    twentyperc = int(round(len(FFT) * 0.02))
    Baseline = np.mean(np.real(FFT[:twentyperc]))
    FFT_corrected = FFT - Baseline
    Re = np.real(FFT_corrected)
    Im = np.imag(FFT_corrected)
    Amp = calculate_amplitude(Re, Im)
    return Amp, Re, Im

def calculate_apodization(Real, Freq):
    # Find sigma at 2% from the max amplitude of the spectra
    Maximum = np.max(np.abs(Real))
    idx_max = np.argmax(np.abs(Real))
    ten_percent = Maximum * 0.02

    b = np.argmin(np.abs(Real[idx_max:] - ten_percent))
    sigma_ap = Freq[idx_max + b]

    apodization_function_s = np.exp(-(Freq / sigma_ap) ** 6)

    Real_apod = Real * apodization_function_s
    
    return Real_apod

def calculate_M2(FFT_real, Frequency):
    # NoClass
    # Take the integral of the REAL PART OF FFT by counts
    Integral = np.trapz(np.real(FFT_real))
    
    # Normalize FFT to the Integral value
    Fur_normalized = np.real(FFT_real) / Integral
    
    # Calculate the integral of normalized FFT to receive 1
    Integral_one = np.trapz(Fur_normalized)
    
    # Multiplication (the power ^n will give the nth moment (here it is n=2)
    Multiplication = (Frequency ** 2) * Fur_normalized
    
    # Calculate the integral of multiplication - the nth moment
    # The (2pi)^2 are the units to transform from rad/sec to Hz
    # ppbly it should be (2pi)^n for generalized moment calculation
    M2 = (np.trapz(Multiplication)) * 4 * np.pi ** 2
    
    # Check the validity
    if np.abs(np.mean(Multiplication[0:10])) > 10 ** (-6):
        print('Apodization is wrong!')

    if M2 < 0:
        M2 = 0
        T2 = 0
    else:
        T2 = np.sqrt(2/M2)
    
    return M2, T2

def decaying_exponential(x, a, b, c):
    return a * np.exp(-x/b) + c

def calculate_SFC(Amplitude):
    S = np.mean(Amplitude[1:4])
    L = np.mean(Amplitude[50:70])
    SFC = (S-L)/S
    return SFC

def calculate_amplitude(Real, Imaginary):
    # NoClass
    Amp = np.sqrt(Real ** 2 + Imaginary ** 2)
    return Amp

def calculate_frequency_scale(Time):
    numberp = len(Time)

    dt = Time[1] - Time[0]
    f_range = 1 / dt
    f_nyquist = f_range / 2
    df = 2 * (f_nyquist / numberp)
    Freq = np.arange(-f_nyquist, f_nyquist + df, df)

    return Freq

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

file_path = r'C:\Mega\NMR\003_Temperature\2024_03_01_Chocolates\2024_03_07_Chocolate_85per_MSE_SE\SE_Ch85_70_c.dat'
#Time, Amp, Re, Im = process_file_data(r'C:\Mega\NMR\003_Temperature\2023_12_05_CPMG_SE_Temperature_choco_85_2\SE\SE_Chocolate85percent_20_c.dat')
#Time, Amp, Re, Im = process_file_data(r'C:\Mega\NMR\003_Temperature\2023_12_21_SE_Temperature_PS35000\SE_PS_70_c.dat')

file_path_gly = r'C:\Mega\NMR\003_Temperature\2024_03_08_SE_Temperature_Glycerol\SE_Glycerol_70_c.dat'

Time_s, Re_s, Im_s = analysis_time_domain(file_path)
Time_r, Re_r, Im_r = analysis_time_domain(file_path_gly)
Time_cr, Re_n, Im_n = long_component(Time_s, Time_r, Re_s, Re_r, Im_s, Im_r)
M2, T2 = FFT_processing(Time_cr, Re_n, Im_n)

print(M2)
print(T2)

print('end')