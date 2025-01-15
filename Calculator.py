# This is a distinct calculation part

# This Python file uses the following encoding: utf-8

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter



# Math procedures
def _find_nearest(array, value):
    # private api
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def read_data(file_path, header):
    data = np.loadtxt(file_path, skiprows=header)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    return x, y, z

def normalize_mq(DQ, MQ, state):
    maximum = np.max(MQ)
    DQ_norm = DQ/maximum
    MQ_norm = MQ/maximum

    if state == 'minus':
        End = MQ_norm - DQ_norm
    else:
        End = MQ_norm + DQ_norm

    return DQ_norm, MQ_norm, End

def dqmq(file_path, fit_from, fit_to, p, noise_level):

    def exponent(x, a, b, c):
        return a * np.exp(-power*(x/b)) + c

    Time, DQ, Ref = read_data(file_path, 1)

    power = p

    DQ_norm, Ref_norm, Diff = normalize_mq(DQ, Ref, 'minus')

    idx_min = _find_nearest(Time, fit_from)
    idx_max = _find_nearest(Time, fit_to)
    Time_cut = Time[idx_min:idx_max+1]
    Diff_cut = Diff[idx_min:idx_max+1]

    popt, _ = curve_fit(exponent, Time_cut, Diff_cut, p0=[(0, 10, 0)], maxfev=10000000)
    fitted_curve = exponent(Time, *popt)

    # Subtract difference function from Ref
    MQ = Ref_norm - fitted_curve

    DQ_normal, MQ_normal, _ = normalize_mq(DQ_norm, MQ, 'plus')

    additive_function = _exp_apodization(Time)
    nDQ = (DQ_normal+noise_level*additive_function)/(DQ_normal+MQ+2*noise_level*additive_function)

    nDQ = np.insert(nDQ, 0, 0)
    Time0 = np.insert(Time, 0, 0)

    return Time, DQ_norm, Ref_norm, Diff, DQ_normal, MQ_normal, Time0, nDQ, fitted_curve

def _cut_beginning(Time, Data, Data2):
    # private api
    Time_plot = Time[np.argmax(Data):]
    Data_plot = Data[np.argmax(Data):]
    Data_plot2 = Data2[np.argmax(Data):]
    return Time_plot, Data_plot, Data_plot2

def analysis_time_domain(file_path, file_empty, subtract):
    # 1. Read data
    Time, Real, Imag = read_data(file_path, 0)

    if subtract == True:
        Time_empty, Real_empty, Imag_empty = read_data(file_empty, 0)

        if len(Time) < len (Time_empty):
            Real = Real - Real_empty[:Real[-1:]]
            Imag = Imag - Imag_empty[:Imag[-1:]]
        elif len(Time)>len(Time_empty):
            Real = Real[:Real_empty[-1:]] - Real_empty
            Imag = Imag[:Imag_empty[-1:]] - Imag_empty
        elif len(Time) == len(Time_empty):
            Real = Real- Real_empty
            Imag = Imag - Imag_empty
        else:
            print("Couldn't subtract")
            return
    # 2. Crop time below zero
    T_cr, R_cr, I_cr = _crop_time_zero(Time, Real, Imag)

    # 2.5 Crop the ascending part of FID
    T_cr, R_cr, I_cr = _cut_beginning(T_cr, R_cr, I_cr)

    # 3. Phase the data
    R_ph, I_ph = _time_domain_phase(R_cr, I_cr)

    # 4. Adjust Frequency
    # 4.1 Calculate Freq
    Frequency = _calculate_frequency_scale(T_cr)
    # 4.2 Shift Freq
    R_sh, I_sh = _adjust_frequency(Frequency, R_ph, I_ph)

    # 5 Again Phase
    R_ph2, I_ph2 = _time_domain_phase(R_sh, I_sh)

    # Again frequency
    R_sh, I_sh = _adjust_frequency(Frequency, R_ph2, I_ph2)


    return T_cr, R_sh, I_sh

def _reference_long_component(Time, Component, Amplitude_gly, coeff):
    # private api
    # 2. Normalize (reference) components to Amplitude of the reference
    if Amplitude_gly is not None:
        Component_n = Component/Amplitude_gly
    else:
        Component_n = Component

    # 3. Cut the ranges for fitting
    minimum = _find_nearest(Time, 80)
    maximum = _find_nearest(Time, 200)
    # TODO: The user can adjust these numbers

    Time_range = Time[minimum:maximum]
    Component_n_range = Component_n[minimum:maximum]

    # 7. Fit data to exponential decay
    b = ([-max(Component_n)*10, 50, -30], [max(Component_n)*10, 800, 30])
    popt, _      = curve_fit(decaying_exponential, Time_range, Component_n_range, bounds = b, p0=coeff)

    # 8. Set the ranges for subtraction
    Time_cropped  = Time[0:maximum]
    Component_c   = Component_n[0:maximum]

    # 9. Calculate the curves fitted to data within the desired range
    Component_f = decaying_exponential(Time_cropped, *popt)

    # 10. Subtract
    Component_sub = Component_c - Component_f

    return Component_sub, Time_cropped

def subtract_long_component(Time, Re, Im):
    coeff_re = [max(Re), 500, 0.1]

    Real_smooth = savgol_filter(Re, 30, 1)

    Real_subtracted, Time_cropped   = _reference_long_component(Time, Real_smooth, None, coeff_re)
    Im_subtracted = np.zeros(len(Real_subtracted))
    return Time_cropped, Real_subtracted, Im_subtracted

def long_component(Time_s, Time_r, Re_s, Re_r, Im_s, Im_r):
    # r stands for reference, s stands for sample
    Amp_r = _calculate_amplitude(Re_r, Im_r)
    Amp_s = _calculate_amplitude(Re_s, Im_s)

    # 1. Crop the arrays together (they should be of the same length, but I know, I know...)
    if len(Time_s) > len(Time_r):
        Time  =   Time_s[:len(Time_r)]
        Amp_s   =   Amp_s[:len(Time_r)]
        Re_s    =   Re_s[:len(Time_r)]
        Im_s    =   Im_s[:len(Time_r)]
    else:
        Time  =   Time_r[:len(Time_s)]
        Amp_r   =   Amp_r[:len(Time_s)]

    coeff_re = [0.9, 400, 0.1]
    coeff_im = [1, 400, 0]

    Real_subtracted, Time_cropped   = _reference_long_component(Time, Re_s, Amp_r, coeff_re)
    Im_subtracted = np.zeros(len(Real_subtracted))

    # 11. Normalize
    Re_n, Im_n = _normalize(Real_subtracted, Im_subtracted)

    return Time_cropped, Re_n, Im_n

def final_analysis_time_domain(Time, Real, Imaginary, number_of_points):
    # 5. Apodize the time-domain
    Re_ap, Im_ap = _apodization(Time, Real, Imaginary)

    # 6. Add zeros
    Tim, Fid = _add_zeros(Time, Re_ap, Im_ap, number_of_points)

    #stophere
    return Tim, Fid

def frequency_domain_analysis(FFT, Frequency):

    # 8. Simple baseline
    _, Re, _ = _simple_baseline_correction(FFT)

    # 9. Apodization
    Real_apod = _calculate_apodization(Re, Frequency)

    # 10. M2 & T2
    M2, T2 = _calculate_M2(Real_apod, Frequency)

    return M2, T2

def read_data(file_path, header):
    data = np.loadtxt(file_path, skiprows=header)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    return x, y, z

def _crop_time_zero(Time, Real, Imaginary):
    # private api
    if Time[0] < 0:
        Time_start = 0
        Time_crop_idx = np.where(Time >= Time_start)[0][0]
        Time_cropped = Time[Time_crop_idx:]
        Real_cropped = Real[Time_crop_idx:]
        Imaginary_cropped = Imaginary[Time_crop_idx:]
        return Time_cropped, Real_cropped, Imaginary_cropped
    else:
        return Time, Real, Imaginary

def _time_domain_phase(Real, Imaginary):
    # private api
    delta = np.zeros(360)

    for phi in range(360):
        Re_phased = Real * np.cos(np.deg2rad(phi)) - Imaginary * np.sin(np.deg2rad(phi))
        Im_phased = Real * np.sin(np.deg2rad(phi)) + Imaginary * np.cos(np.deg2rad(phi))
        Magnitude_phased = _calculate_amplitude(Re_phased, Im_phased)
        
        Re_cut = Re_phased[:5]
        Ma_cut = Magnitude_phased[:5]
        
        delta[phi] = np.mean(Ma_cut - Re_cut)

    idx = np.argmin(delta)
    #print(idx)

    Re = Real * np.cos(np.deg2rad(idx)) - Imaginary * np.sin(np.deg2rad(idx))
    Im = Real * np.sin(np.deg2rad(idx)) + Imaginary * np.cos(np.deg2rad(idx))

    return Re, Im

def _adjust_frequency(Frequency, Re, Im):
    # private api
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
    index_zero = _find_nearest(Frequency, 0)

    # Find difference
    delta_index = index_max - index_zero

    # To avoid over correction
    if delta_index == 0:
        return Re, Im

    # Shift the spectra (amplitude) by the difference in indices
    FFT_shifted = np.concatenate((FFT[delta_index:], FFT[:delta_index]))

    # iFFT
    Fid_shifted = np.fft.ifft(np.fft.fftshift(FFT_shifted))

    # Define Real, Imaginary and Amplitude
    Re_shifted = np.real(Fid_shifted)
    Im_shifted = np.imag(Fid_shifted)

    return Re_shifted, Im_shifted

def _normalize(Real, Imaginary):
    # private api
    Amplitude = np.sqrt(Real ** 2 + Imaginary ** 2)
    Amplitude_max = np.max(Amplitude)
    Amp = Amplitude/Amplitude_max
    Re = Real/Amplitude_max
    Im = Imaginary/Amplitude_max
    return Re, Im

def _apodization(Time, Real, Imaginary):
    # private api
    Amplitude = _calculate_amplitude(Real, Imaginary)
    coeffs = np.polyfit(Time, Amplitude, 1)  # Fit an exponential decay function
    c = np.polyval(coeffs, Time)
    d = np.argmin(np.abs(c - 1e-5))
    sigma = Time[d]
    if sigma == 0:
        sigma = 1000
    apodization_function = np.exp(-(Time / sigma) ** 4)
    Re_ap = Real * apodization_function
    Im_ap = Imaginary * apodization_function
    return Re_ap, Im_ap

def _exp_apodization(Time):
    # private api
    Time_end = Time[-1]+10
    noise = Time_end*0.3
    additive_function =  0.5*np.exp(((Time - Time_end) / noise) ** 3)
    return additive_function

def _add_zeros(Time, Real, Imaginary, number_of_points):
    # private api
    length_diff = number_of_points - len(Time)
    amount_to_add = np.zeros(length_diff+1)

    Re_zero = np.concatenate((Real, amount_to_add))
    Im_zero = np.concatenate((Imaginary, amount_to_add))

    dt = Time[1] - Time[0]
    Time_to_add = Time[-1] + np.arange(1, length_diff + 1) * dt

    Time = np.concatenate((Time, Time_to_add))
    Fid = np.array(Re_zero + 1j * Im_zero)
    Fid = Fid[:-1]

    return Time, Fid

def _simple_baseline_correction(FFT):
    # private api
    twentyperc = int(round(len(FFT) * 0.02))
    Baseline = np.mean(np.real(FFT[:twentyperc]))
    FFT_corrected = FFT - Baseline
    Re = np.real(FFT_corrected)
    Im = np.imag(FFT_corrected)
    Amp = _calculate_amplitude(Re, Im)
    return Amp, Re, Im

def _calculate_apodization(Real, Freq):
    # private api
    # Find sigma at 2% from the max amplitude of the spectra
    Maximum = np.max(np.abs(Real))
    idx_max = np.argmax(np.abs(Real))
    ten_percent = Maximum * 0.02

    b = np.argmin(np.abs(Real[idx_max:] - ten_percent))
    sigma_ap = Freq[idx_max + b]

    apodization_function_s = np.exp(-(Freq / sigma_ap) ** 6)

    Real_apod = Real * apodization_function_s
    
    return Real_apod

def _calculate_amplitude(Real, Imaginary):
    # private api
    Amp = np.sqrt(Real ** 2 + Imaginary ** 2)
    return Amp

def _calculate_frequency_scale(Time):
    # private api
    numberp = len(Time)

    dt = Time[1] - Time[0]
    f_range = 1 / dt
    f_nyquist = f_range / 2
    df = 2 * (f_nyquist / numberp)
    Freq = np.arange(-f_nyquist, f_nyquist + df, df)
    Freq = Freq[:-1]

    return Freq

def _find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def _calculate_M2(FFT_real, Frequency):
    # private api
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

def calculate_SC(Amplitude):
    S = np.mean(Amplitude[2:20])
    L = np.mean(Amplitude[120:160])
    # solid_content = (S-L)/S
    solid_content = L
    return solid_content

def calculate_DQ_intensity(Time, Amplitude):
    idx_time = np.argmin(np.abs(Time - 4))
    DQ = np.mean(Amplitude[:idx_time])
    return DQ

def decaying_exponential(x, a, b, c):
    return a * np.exp(-x/b) + c

def decaying_2exponential(x, a1, b1, a2, b2, c):
    return a1 * np.exp(-x/b1) + a2 * np.exp(-x/b2) + c

def decaying_3exponential(x, a1, b1, a2, b2, a3, b3, c):
    return a1 * np.exp(-x/b1) + a2 * np.exp(-x/b2) + a3 * np.exp(-x/b3) + c

def fit_exponent(Time, Signal, order):
    Time_fit = np.linspace(start=min(Time), stop=max(Time), num=60)

    if order == 1:
        p = [-10, 100, 15]
        b = ([-np.inf, 0, -np.inf], [np.inf, 50000, np.inf])
        popt_, _ = curve_fit(decaying_exponential, Time, Signal, p0 = p, bounds = b, maxfev=10000000)
        fitted_curve = decaying_exponential(Time_fit, *popt_)
        popt = np.round(popt_, 3)

        tau1 = popt[1]
        A1   = popt[0]

        tau2 = 0
        tau3 = 0
        A2   = 0
        A3   = 0
        decrease_order = False

        r2_curve = decaying_exponential(Time, *popt_)
        R2 = round(calculate_r_squared(Signal, r2_curve),4)

    elif order == 2:
        b=([-np.inf, 0, -np.inf, 0, -np.inf], [np.inf, 50000, np.inf, 50000, np.inf])
        try:
            p1 = [-10, 100, 15]
            b1 = ([-np.inf, 0, -np.inf], [np.inf, 50000, np.inf])
            popt1, _ = curve_fit(decaying_exponential, Time, Signal, p0 = p1, bounds = b1, maxfev=10000000)
            p = [popt1[0], popt1[1], 10, 10, popt1[2]]
            popt_, _ = curve_fit(decaying_2exponential, Time, Signal, bounds = b, maxfev=10000000, p0=p)
        except:
            print('No covariance for exp fitting of the 2d order')
            return
        fitted_curve = decaying_2exponential(Time_fit, *popt_)
        popt = np.round(popt_, 3)

        A1_, tau1_, A2_, tau2_ = popt[0:4]
        tau_A_pairs = sorted([(tau1_, A1_), (tau2_, A2_)])
        tau1, A1 = tau_A_pairs[0]
        tau2, A2 = tau_A_pairs[1]

        tau3 = 0
        A3   = 0

        decrease_order = check_tau_values(tau1, tau2, tau3)

        r2_curve = decaying_2exponential(Time, *popt_)
        R2 = round(calculate_r_squared(Signal, r2_curve),4)

    elif order ==3:
        # p = [-10, 10, -10, 50, -10, 100, 15]
        b=([-np.inf, 0, -np.inf, 0, -np.inf, 0, -np.inf], [np.inf, 50000, np.inf, 50000, np.inf, 50000, np.inf])
        # popt, pcov = curve_fit(decaying_3exponential, Time, Signal, p0 = p, bounds = b, maxfev=10000000)
        popt_, _ = curve_fit(decaying_3exponential, Time, Signal, bounds = b, maxfev=10000000)

        fitted_curve = decaying_3exponential(Time_fit, *popt_)

        popt = np.round(popt_, 3)
        A1_, tau1_, A2_, tau2_ , A3_, tau3_ = popt[0:6]
        tau_A_pairs = sorted([(tau1_, A1_), (tau2_, A2_), (tau3_, A3_)])
        tau1, A1 = tau_A_pairs[0]
        tau2, A2 = tau_A_pairs[1]
        tau3, A3 = tau_A_pairs[2]

        decrease_order = check_tau_values(tau1, tau2, tau3)

        r2_curve = decaying_3exponential(Time, *popt_)
        R2 = round(calculate_r_squared(Signal, r2_curve),4)

    return Time_fit, fitted_curve, tau1, tau2, tau3, R2, A1, A2, A3, decrease_order

def check_tau_values(tau1, tau2, tau3):
    decrease_order = False

    if tau2 == 0:
        return
    else:
        tau_rounded = round(tau1, 1)
        tau2_rounded = round(tau2, 1)
        tau3_rounded = round(tau3, 1)
        
        # Check if the rounded values are equal
        if (tau2_rounded == tau3_rounded) or (tau_rounded == tau2_rounded):
            decrease_order = True
        else:
            decrease_order = False
    
    return decrease_order

def calculate_r_squared(y_true, y_pred):

    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    # Calculate the mean of the observed values
    y_mean = np.mean(y_true)

    # Calculate the total sum of squares (TSS)
    ss_tot = np.sum((y_true - y_mean) ** 2)

    # Calculate the residual sum of squares (RSS)
    ss_res = np.sum((y_true - y_pred) ** 2)

    # Calculate R-squared
    r_squared = 1 - (ss_res / ss_tot)

    return r_squared

def gaussian(x, amp, cen, wid, y0):
    # NoClass
    return amp * np.exp(-(x - cen)**2 / (2 * wid**2)) + y0

def lorenz(x, amp, cen, wid, y0):
    # NoClass
    return (amp * (wid**2)) / ((x - cen)**2 + (wid**2)) + y0

def voigt(x, amp, cen, wid, frac, y0):
    lorentzian = amp *(2 * wid) / (np.pi * (4 * (x - cen)**2 + wid**2))
    gaussian = amp *(np.exp((-4 * np.log(2) * (x - cen)**2) / wid**2)) / (wid * np.sqrt(np.pi / (4 * np.log(2))))
    return  (frac * lorentzian + (1 - frac) * gaussian) + y0

def twod_model(x, CDD, tauc, A, Ctrans, tautrans, taures):
    term1 = tauc / (1 + (np.pi * 2 * x * tauc) ** 2)
    term2 = 4 * tauc / (1 + (2 * np.pi * 2 * x * tauc) ** 2)
    log_term1 = np.log((1 + (np.pi * 2 * x * tautrans) ** 2) / ((tautrans / taures) ** 2 + (np.pi * 2 * x * tautrans) ** 2))
    log_term2 = 4 * np.log((1 + (np.pi * 2 * 2 * x * tautrans) ** 2) / ((tautrans / taures) ** 2 + (np.pi * 2 * 2 * x * tautrans) ** 2))
    
    result = CDD * (term1 + term2) + A + Ctrans * tautrans * (log_term1 + log_term2)
    
    return result

def fit_model(Omega, Rate, fixed_CDD, initial_parameters):
    Omega_fit = np.arange(min(Omega), max(Omega) + 0.1, 0.05)

    bs = ([0, 0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])

    if fixed_CDD is not None and initial_parameters is not None:
        initial_parameter = initial_parameters[1:]
        b = bs[0][1:],bs[1][1:]
        def twod_model_fixed(x, tauc, A, Ctrans, tautrans, taures):
            return twod_model(x, fixed_CDD, tauc, A, Ctrans, tautrans, taures)
        popt, _ = curve_fit(twod_model_fixed, Omega, Rate, p0 = initial_parameter, bounds = b, maxfev=100000080)
        popt = [fixed_CDD] + popt.tolist()
    else:
        popt, _ = curve_fit(twod_model, Omega, Rate, bounds = bs, maxfev=100000000)

        
    fitted_curve = twod_model(Omega_fit, *popt)

    r2_curve = twod_model(Omega, *popt)
    R2 = round(calculate_r_squared(Rate, r2_curve),4)

    return Omega_fit, fitted_curve, popt, R2

def repeating_part(t, x):
    return t / (1 + (np.pi * 2 * x * t)**2) + 4 * t / (1 + (np.pi * 4 * x * t)**2)

def simplified_expression1(x, C1, C2, C3, A, t1, t2, t3):
    return C1 * repeating_part(t1, x) + C2 * repeating_part(t2, x) + C3 * repeating_part(t3, x) + A

def simplified_expression(x, C1, C2, C3, A, t1, t2, t3, alpha):
    return C1 * repeating_part(t1, x) + C2 * repeating_part(t2, x) + C3 * repeating_part(t3, x) + A * (x)**(-alpha)

def simulation(Omega, Rate, state):
    if state == 'One':
        initial = [1, 0.01, 1, 0.5]
        bounds = ([0, 0, 0, -2], [np.inf, np.inf, np.inf, 2])
        popt, _ = curve_fit(lambda x, C1, A, t1, alpha: simplified_expression(x, C1, 0, 0, A, t1, 0, 0, alpha), Omega, Rate, p0=initial, bounds=bounds)
        C1, A, t1, alpha = popt[0], popt[1], popt[2], popt[3]
        fitting = simplified_expression(Omega, C1, 0, 0, A, t1, 0, 0, alpha)
        short, middle, long = 0, 0, 0
        R2 = round(calculate_r_squared(Rate, fitting), 4)

    elif state == 'Two':
        initial = [1, 0.01, 1, 0.01, 0.1, 0.5]
        bounds = ([0, 0, 0, 0, 0, -2], [np.inf, np.inf, np.inf, np.inf, np.inf, 2])
        popt, _ = curve_fit(lambda x, C1, C2, A, t1, t2, alpha: simplified_expression(x, C1, C2, 0, A, t1, t2, 0, alpha), Omega, Rate, p0=initial, bounds=bounds)
        C1, C2, A, t1, t2, alpha = popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]
        fitting = simplified_expression(Omega, C1, C2, 0, A, t1, t2, 0, alpha)
        short = simplified_expression(Omega, C1, 0, 0, A, t1, 0, 0, alpha)
        middle = simplified_expression(Omega, 0, C2, 0, A, 0, t2, 0, alpha)
        long = 0
        R2 = round(calculate_r_squared(Rate, fitting), 4)

    elif state == 'Three':
        initial = [1, 0.1, 0.01, 1, 0.01, 0.1, 1, 0.5]
        bounds = ([0, 0, 0, 0, 0, 0, 0, -2], [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 2])
        popt, _ = curve_fit(simplified_expression, Omega, Rate, p0=initial, bounds=bounds)
        C1, C2, C3, A, t1, t2, t3, alpha = popt
        fitting = simplified_expression(Omega, C1, C2, C3, A, t1, t2, t3, alpha)
        short = simplified_expression(Omega, C1, 0, 0, A, t1, 0, 0, alpha)
        middle = simplified_expression(Omega, 0, C2, 0, A, 0, t2, 0, alpha)
        long = simplified_expression(Omega, 0, 0, C3, A, 0, 0, t3, alpha)
        R2 = round(calculate_r_squared(Rate, fitting), 4)

    else:
        return  # Invalid state

    return fitting, short, middle, long, popt, R2
