# This Python file uses the following encoding: utf-8
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import os, re
from scipy.interpolate import griddata
from matplotlib import cm


# Find nearest value in array
def find_nearest(array, value):
    idx = (np.abs(np.asarray(array) - value)).argmin()
    return idx

# Gaussian functions
def gauss(x, A, sigma, y0):
    return A * np.exp(-x**2 / (2 * sigma**2)) + y0

def gauss_fixed_sum(x, C, sigma, y0):
    A = C - y0
    return A * np.exp(-x**2 / (2 * sigma**2)) + y0

def gauss_const_ampl(amplitude):
    return lambda x, sigma, y0: gauss(x, amplitude-y0, sigma, y0)

#Polynom functions
# 4th power
def polynom4(x, A, c, g):
    return A + c * x**2 + g * x**4

def polynom4_const_ampl(amplitude):
    return lambda x, c, g: polynom4(x, amplitude, c, g)

# 6th power
def polynom6(x, A, c, g, h):
    return A +  c * x**2 + g * x**4 + h * x**6

def polynom6_const_ampl(amplitude):
    return lambda x, c, g, h: polynom6(x, amplitude, c, g, h)

# 8th power
def polynom8(x, A, c, g, h, j):
    return A +  c * x**2 + g * x**4 + h * x**6 + j * x**8

def polynom8_const_ampl(amplitude):
    return lambda x, c, g, h, j: polynom8(x, amplitude, c, g, h, j)

# Exp function
def decaying_exponential(x, a, b, c):
    return a * np.exp(-x/b) + c

# Read data from file
def read_data(file_path):
    data = np.loadtxt(file_path)
    return data[:, 0], data[:, 1], data[:, 2]

# Calculate frequency from time
def calculate_frequency_scale(Time):
    numberp = len(Time)

    dt = Time[1] - Time[0]
    f_range = 1 / dt
    f_nyquist = f_range / 2
    df = 2 * (f_nyquist / numberp)
    Freq = np.arange(-f_nyquist, f_nyquist + df, df)
    Freq = Freq[:-1]

    return Freq

# Calculate amplitude
def calculate_amplitude(real, imaginary):
    return np.sqrt(real**2 + imaginary**2)

# Smooth noisy signal
def smooth_noisy_signal(Signal, SmoothOrder, SmoothWindow):

    Signal_smooth = savgol_filter(Signal, SmoothWindow, SmoothOrder)

    return Signal_smooth

# Adjust phase
def adjust_phase(Real, Imaginary):
    delta = np.zeros(360)

    for phi in range(360):
        Re_phased = Real * np.cos(np.deg2rad(phi)) - Imaginary * np.sin(np.deg2rad(phi))
        Im_phased = Real * np.sin(np.deg2rad(phi)) + Imaginary * np.cos(np.deg2rad(phi))
        Magnitude_phased = calculate_amplitude(Re_phased, Im_phased)

        Re_cut = Re_phased[:5]
        Ma_cut = Magnitude_phased[:5]

        delta[phi] = np.mean(Ma_cut - Re_cut)

    idx = np.argmin(delta)

    Re = Real * np.cos(np.deg2rad(idx)) - Imaginary * np.sin(np.deg2rad(idx))
    Im = Real * np.sin(np.deg2rad(idx)) + Imaginary * np.cos(np.deg2rad(idx))

    return Re, Im

def voigt(x, amp, cen, wid, frac, y0):
    lorentzian = amp *(2 * wid) / (np.pi * (4 * (x - cen)**2 + wid**2))
    gaussian = amp *(np.exp((-4 * np.log(2) * (x - cen)**2) / wid**2)) / (wid * np.sqrt(np.pi / (4 * np.log(2))))
    return  (frac * lorentzian + (1 - frac) * gaussian) + y0

def fit_fft_with_voigh(Frequency, FFT):
    spectrum = np.abs(FFT)
    # coarse center guess from max
    idx_peak = np.argmax(spectrum)
    mu0 = Frequency[idx_peak]

    # window in Hz
    window = (Frequency[-1] - Frequency[0]) * 0.1  # 10% of full span
    mask = (Frequency >= (mu0 - window/2)) & (Frequency <= (mu0 + window/2))

    x = Frequency[mask]
    y = spectrum[mask]

    #  initial guesses
    cen0 = Frequency[idx_peak]   # center guess
    amp0 = np.max(spectrum)
    wid0 = (Frequency[-1] - Frequency[0]) / 50  # +-200 khz
    frac0 = 0.5     # start with equal mix
    y00 = 0

    p0 = [amp0, cen0, wid0, frac0, y00]

    # amplitude, center, width, fraction, y0
    lower = [0, Frequency[0], (Frequency[1]-Frequency[0]), 0, -10]
    upper = [np.inf, Frequency[-1], (Frequency[-1]-Frequency[0]), 1, 10]

    #  fit
    popt, pcov = curve_fit(voigt, x, y, bounds=(lower, upper), p0=p0, maxfev=20000)

    return popt

# Adjust frequency
def adjust_frequency(Frequency, Re, Im):
    Fid_unshifted = np.array(Re + 1j * Im)

    # FFT
    FFT = np.fft.fftshift(np.fft.fft(Fid_unshifted))

    # Check the length of FFT and Frequency (it is always the same, this is just in case)
    if len(Frequency) != len(FFT):
        Frequency = np.linspace(Frequency[0], Frequency[-1], len(FFT))

    popt_fitting = []
    try:
        popt_fitting = fit_fft_with_voigh(Frequency, FFT)
        index_max = np.argmax(voigt(Frequency, *popt_fitting))
    except:
        print('The fitting of the FFT to adjust the frequency couldnt be done. Using simple adjustemnt')
        index_max = np.argmax(FFT)

    # Find index of zero (frequency)
    index_zero = find_nearest(Frequency, 0)

    # Find difference
    delta_index = index_max - index_zero
    # print(f'The max is {Frequency[index_max]}')

    # To avoid over correction
    if delta_index == 0:
        print('The frequency was not adjusted, it is already perfect')
        return Re, Im

    # Shift the spectra (amplitude) by the difference in indices
    FFT_shifted = np.concatenate((FFT[delta_index:], FFT[:delta_index]))

    # iFFT
    Fid_shifted = np.fft.ifft(np.fft.ifftshift(FFT_shifted))

    # Define Real, Imaginary and Amplitude
    Re_shifted = np.real(Fid_shifted)
    Im_shifted = np.imag(Fid_shifted)

    # plt.plot(Frequency, FFT, label="original")
    # plt.plot(Frequency, FFT_shifted, label="shifted")
    # plt.plot(Frequency, voigt(Frequency, *popt_fitting), 'r--', label="Gaussian fit")

    # plt.show()

    return Re_shifted, Im_shifted

# Correct spectra with baseline
def simple_baseline_correction(FFT):
    twentyperc = int(round(len(FFT) * 0.02))
    Baseline = np.mean(np.real(FFT[:twentyperc]))
    FFT_corrected = FFT - Baseline
    Re = np.real(FFT_corrected)
    return Re

# Apodize spectra
def apodization_fft(Real, Freq):
    # Find sigma at 0.1% from the max amplitude of the spectra
    Maximum = np.max(np.abs(Real))
    idx_max = np.argmax(np.abs(Real))
    ten_percent = Maximum * 0.03

    b = np.argmin(np.abs(Real[idx_max:] - ten_percent))
    sigma_ap = Freq[idx_max + b]

    apodization_function_s = np.exp(-(Freq / sigma_ap) ** 2)

    Real_apod = Real * apodization_function_s

    # plt.plot(Freq, Real, 'r', label = 'Original')
    # plt.plot(Freq, apodization_function_s, 'k--', label = 'Apodization')
    # plt.plot(Freq, Real_apod, 'b', label = 'Apodized')
    # plt.legend()
    # plt.show()

    return Real_apod

# Zero filling procedure
def add_zeros(Time, Real, Imaginary, number_of_points):
    length_diff = number_of_points - len(Time)
    amount_to_add = np.zeros(length_diff+1)

    Re_zero = np.concatenate((Real, amount_to_add))
    Im_zero = np.concatenate((np.zeros(len(Real)), amount_to_add))

    dt = Time[1] - Time[0]
    Time_to_add = Time[-1] + np.arange(1, length_diff + 1) * dt

    Time = np.concatenate((Time, Time_to_add))
    Fid = np.array(Re_zero + 1j * Im_zero)
    Fid = Fid[:-1]

    return Time, Fid

# Apodization of time Domain
def apodization(Time, Real, Imaginary, sigma):
    apodization_function = np.exp(-(Time / sigma) ** 4)
    Re_ap = Real * apodization_function
    Im_ap = Imaginary * apodization_function

    # plt.plot(Time, apodization_function, 'r--', label='apodization')
    # plt.plot(Time, Real, 'k', label='Re')
    # plt.plot(Time, Re_ap, 'k--', label='Re ap')
    # plt.xlim([-5,80])
    # plt.xlabel('Time, μs')
    # plt.ylabel('Amplitude, a.u.')
    # plt.legend()
    # plt.tight_layout()
    # plt.show()

    return Re_ap, Im_ap

# Cut the first part of the FID
def cut_beginning(Time, Data):
    Time_plot = Time[np.argmax(Data):]
    Data_plot = Data[np.argmax(Data):]
    return Time_plot, Data_plot

# Normalize the data to FID at long times
def normalize_to_fid(Fid, Data, Time_fid, Time_data, fr, to):
    Mean_amp_fid_long = np.mean(Fid[find_nearest(Time_fid, fr):find_nearest(Time_fid, to)])
    Mean_amp_dat_long = np.mean(Data[find_nearest(Time_data, fr):find_nearest(Time_data, to)])

    difference = Mean_amp_dat_long - Mean_amp_fid_long
    Normalized_amp = Data - difference

    # the division doesn't work very good
    # proportionality_coefficient = Mean_amp_fid_long /Mean_amp_dat_long
    # Normalized_amp = Data * proportionality_coefficient
    return Normalized_amp

# reference the long component
def reference_long_component(Time, Component_n, end):
    # 3. Cut the ranges for fitting
    minimum = find_nearest(Time, end)

    Time_range = Time[minimum:]
    Component_n_range = Component_n[minimum:]

    # Smooth data
    Smooth = savgol_filter(Component_n_range, 40, 0)

    p = [5, 30, 0.5]
    # 7. Fit data to exponential decay
    popt, _      = curve_fit(decaying_exponential, Time_range, Smooth, p0 =p)
    # TODO: If no covariance met, give the error window!

    # 9. Calculate the curves fitted to data within the desired range
    Component_f = decaying_exponential(Time, *popt)

    # 10. Subtract
    Component_sub = Component_n - Component_f

    # # For Debug
    # plt.plot(Time, Component_n, 'r', label='Original')
    # plt.plot(Time, Component_sub, 'b', label='Subtracted')
    # plt.plot(Time, Component_f, 'k--', label='Fitted')
    # plt.xlabel('Time, μs')
    # plt.ylabel('Amplitude, a.u.')
    # plt.legend()
    # plt.tight_layout()
    # plt.show()

    return Component_sub

# Calculate M2
def calculate_M2(FFT_real, Frequency):
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

# Functions that call other multiple functions :)
# Create spectra with option for corrections
def freq_domain_correction(Time, Real, Imaginary, Apodize, time_a, Adjust):
    number_of_points = 2**16

    if Apodize:
        sigma = time_a
        Real, Imaginary = apodization(Time, Real, Imaginary, sigma)

    if Adjust:
        Freq = calculate_frequency_scale(Time)
        Re_ph, Im_ph = adjust_phase(Real, Imaginary)
        Real_fft, Imaginary_fft = adjust_frequency(Freq, Re_ph, Im_ph)
    else:
        # Create complex FID
        FID = np.array(Real + 1j * Imaginary)
        # FFT
        FFT = np.fft.fftshift(np.fft.fft(FID))
        Real_fft = np.real(FFT)
        Imaginary_fft = np.imag(FFT)

    # 6. Add zeros
    Time_zero, Fid_zero = add_zeros(Time, Real_fft, Imaginary_fft, number_of_points)
    Frequency = calculate_frequency_scale(Time_zero)

    # 7. FFT
    FFT = np.fft.fftshift(np.fft.fft(Fid_zero))

    # 8. Simple baseline
    Re = simple_baseline_correction(FFT)

    # 9. Apodization
    Real_apod = apodization_fft(Re, Frequency)

    return Frequency, Real_apod

# Create NMR signal with option for corrections
def time_domain_correction(file_path):
    Time, Re_original, Im_original = read_data(file_path)

    R_phased, _ = adjust_phase(Re_original, Im_original)

    return Time, R_phased

# All the procedures with correction options
def nmr_signal_correction(file_path, file_path_fid, file_path_empty, file_path_empty_fid, Subtraction, Cut_beginning, Normalize_to_fid, n_start, n_end, Long_component, L_end, AdjustZero, zero, Smooth, SmoothOrder, SmoothWindow):
    Time_td, Re_td = time_domain_correction(file_path)
    Time_td_fid, Re_td_fid = time_domain_correction(file_path_fid)

    if Smooth:
        Re_td = smooth_noisy_signal(Re_td, SmoothOrder, SmoothWindow)
        Re_td_fid = smooth_noisy_signal(Re_td_fid, SmoothOrder, SmoothWindow)

    if Subtraction and file_path_empty!=[] and file_path_empty_fid!=[]:
        length_to_cut = len(Time_td)
        _, Re_td_empty, _ = read_data(file_path_empty)
        Re_td = Re_td - Re_td_empty[:length_to_cut]

        length_to_cut = len(Time_td_fid)
        _, Re_td_empty_fid, _ = read_data(file_path_empty_fid)
        Re_td_fid = Re_td_fid - Re_td_empty_fid[:length_to_cut]

    if AdjustZero:
        Time_td_fid = Time_td_fid + zero

    if Cut_beginning:
        _, Re_td          = cut_beginning(Time_td, Re_td)
        Time_td = np.linspace(0, 0+(0.5*len(Re_td)),len(Re_td),endpoint=False)
        Time_td_fid, Re_td_fid  = cut_beginning(Time_td_fid, Re_td_fid)

    if Normalize_to_fid:
        Re_td       = normalize_to_fid(Re_td_fid, Re_td, Time_td_fid, Time_td, n_start, n_end)

    if Long_component:
        Re_td       = reference_long_component(Time_td, Re_td, L_end)
        Re_td_fid   = reference_long_component(Time_td_fid, Re_td_fid, L_end)

    return Time_td, Re_td, Time_td_fid, Re_td_fid

def find_maximum_se(echo_time, maximum, start, end):
    echo_time = np.array(echo_time).astype(float)
    echo_time_to_fit = np.arange(0,echo_time[-1:], 0.01)
    if end == 0:
        end = None
    else:
        end = int(-end)
    echo_time = echo_time[int(start):end]
    maximum = maximum[int(start):end] 
    # Gaussian fit for SE maximum amplitude
    p = [10, 6, 1] # Initial guess
    popt, _ = curve_fit(gauss, echo_time, maximum, p0=p)
    Extrapolation = gauss(0, *popt)
    fitting_curve = gauss(echo_time_to_fit, *popt)

    return Extrapolation, fitting_curve, echo_time_to_fit

def build_up_fid(Time, Data, A, function_to_fit, start, finish):
    Time_cut    = Time[find_nearest(Time, start):find_nearest(Time, finish)]
    Data_cut    = Data[find_nearest(Time, start):find_nearest(Time, finish)]

    delta_time = Time_cut[1]-Time_cut[0]
    Time_build_from_zero = np.arange(0, start, delta_time)

    if function_to_fit == 'Polynom 4':
        # # Fit the small part of the FID with polynom (4 degree)
        popt, _ = curve_fit(polynom4_const_ampl(A), Time_cut, Data_cut, p0=[0.005, 0.005])
        Data_built = polynom4(Time, A, *popt)

        # Build-up the FID from time 0 to the first interception
        Data_build_from_zero = polynom4(Time_build_from_zero, A, *popt)

    elif function_to_fit == 'Polynom 6':
        # # Fit the small part of the FID with polynom (6 degree)
        popt, _ = curve_fit(polynom6_const_ampl(A), Time_cut, Data_cut, p0=[0.005, 0.005, 0.005])
        Data_built = polynom6(Time, A, *popt)

        # Build-up the FID from time 0 to the first interception
        Data_build_from_zero = polynom6(Time_build_from_zero, A, *popt)

    elif function_to_fit == 'Polynom 8':
        # # Fit the small part of the FID with polynom (4 degree)
        popt, _ = curve_fit(polynom8_const_ampl(A), Time_cut, Data_cut, p0=[0.005, 0.005, 0.005, 0.005])
        Data_built = polynom8(Time, A, *popt)

        # Build-up the FID from time 0 to the first interception
        Data_build_from_zero = polynom8(Time_build_from_zero, A, *popt)

    elif function_to_fit == 'Gaussian':
        # Fit the small part of the FID with gauss function with restricted amplitude
        popt, _ = curve_fit(gauss_const_ampl(A), Time_cut, Data_cut, p0=[8, 0])
        Data_built = gauss(Time, A-popt[1], *popt)

        # Build-up the FID from time 0 to the first interception
        Data_build_from_zero = gauss(Time_build_from_zero, A-popt[1], *popt)

    else:
        print('No function')
        return

    # Build the data from start to finish
    # Make an weighted average, where weight depends on X
    # So, in the beginning, no FID, all built ->0
    # In the end, only FID, no built ->1
    # Begin - where the start, end where is the finish

    start_idx = find_nearest(Time, start)
    finish_idx = find_nearest(Time, finish)
    Time_build_middle = Time[start_idx+1:finish_idx]
    length = len(Time_build_middle)
    data_fid = Data[start_idx+1:finish_idx]
    data_built = Data_built[start_idx+1:finish_idx]
    weight = np.linspace(0, 1, length)

    Data_build_middle = weight * data_fid + (1 - weight) * data_built

    # Build the data from 2d interception until the end
    Time_build_end  = Time[finish_idx+1:]
    Data_build_end   = Data[finish_idx+1:]

    Time_build_full = np.concatenate((Time_build_from_zero,Time_build_middle, Time_build_end))
    Data_build_full  = np.concatenate((Data_build_from_zero,Data_build_middle, Data_build_end))

    return Time_build_full, Data_build_full, Data_built

def for_the_sake_of_beauty(Time_td, Re_td, Time_td_fid, Re_td_fid, Apodize, time_a):

    if Apodize:
        sigma = time_a
        Re_td, _ = apodization(Time_td, Re_td, 0, sigma)
        Re_td_fid, _ = apodization(Time_td_fid, Re_td_fid, 0, sigma)

    return Time_td, Re_td, Time_td_fid, Re_td_fid


# ## Save build-up data for export
# df = pd.DataFrame({'Time' : Time_build_full_mse_r, 'Re': Re_build_full_mse})
# df.to_csv("MSE_build_up.dat", sep ='\t', index = 'none')

# df = pd.DataFrame({'Time' : Time_build_full_se_r, 'Re': Re_build_full_se})
# df.to_csv("SE_build_up.dat", sep ='\t', index = 'none')



