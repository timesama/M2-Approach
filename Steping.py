# This Python file uses the following encoding: utf-8
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Step by step time-domain and frequency domain FID processing
# No functions!
def decaying_exponential(x, a, b, c):
    return a * np.exp(-x/b) + c
# I lied

file_path = r'C:\Mega\NMR\003_Temperature\2024_03_01_Chocolates\2024_03_07_Chocolate_85per_MSE_SE\SE_Ch85_30_c.dat'
file_path_glycerol = r'C:\Mega\NMR\003_Temperature\2024_03_08_SE_Temperature_Glycerol\SE_Glycerol_30_c.dat'

# 1. Read the data:
data = np.loadtxt(file_path)
x, y, z = data[:, 0], data[:, 1], data[:, 2]

Time = x
Real = y
Imag = z
Amp = np.sqrt(Real ** 2 + Imag ** 2)

# # Plot the initial data
# plt.figure()
# plt.plot(Time, Real, 'r', label='Initial Re')
# plt.plot(Time, Imag, 'b', label='Initial Im')
# plt.plot(Time, Amp, 'k', label='Initial Amp')
# plt.legend()

# 2. Crop the data from Time = 0:
if Time[0] < 0:
    Time_start = 0
    Time_crop_idx = np.where(Time >= Time_start)[0][0]
    Time = Time[Time_crop_idx:]
    Real = Real[Time_crop_idx:]
    Imag = Imag[Time_crop_idx:]
else:
    pass

# 3. Phase the data:
delta = np.zeros(360)

for phi in range(360):
    Re = Real * np.cos(np.deg2rad(phi)) - Imag * np.sin(np.deg2rad(phi))
    Im = Real * np.sin(np.deg2rad(phi)) + Imag * np.cos(np.deg2rad(phi))
    Ma = np.sqrt(Real ** 2 + Imag ** 2)
    
    # Take the short time range up to ~50 microsec
    Re_cut = Re[:100]
    Ma_cut = Ma[:100]
    
    delta[phi] = np.mean(Ma_cut - Re_cut)

idx = np.argmin(delta)

Re_phased = Real * np.cos(np.deg2rad(idx)) - Imag * np.sin(np.deg2rad(idx))
Im_phased = Real * np.sin(np.deg2rad(idx)) + Imag * np.cos(np.deg2rad(idx))
Am_phased = np.sqrt(Re_phased ** 2 + Im_phased ** 2)

# # Plot the phased data
# plt.figure()
# plt.plot(Time, Re_phased, 'r', label='Phased Re')
# plt.plot(Time, Im_phased, 'b', label='Phased Im')
# plt.plot(Time, Am_phased, 'k', label='Phased Amp')
# plt.legend()


# 4. Adjust the frequency
# Calculate Frequency
numberp = len(Time)

dt = Time[1] - Time[0]
f_range = 1 / dt
f_nyquist = f_range / 2
df = 2 * (f_nyquist / numberp)
Frequency = np.arange(-f_nyquist, f_nyquist + df, df)

# Create complex FID
Fid_unshifted = np.array(Re_phased + 1j * Im_phased)

# FFT
FFT = np.fft.fftshift(np.fft.fft(Fid_unshifted))

# Check the length of FFT and Frequency (it is always the same, this is just in case)
if len(Frequency) != len(FFT):
    Frequency = np.linspace(Frequency[0], Frequency[-1], len(FFT))

# Find index of max spectrum (amplitude)
index_max = np.argmax(FFT)

# Find index of zero (frequency)
array = np.asarray(Frequency)
index_zero = (np.abs(array - 0)).argmin()

# Find difference
delta_index = index_max - index_zero

# Shift the spectra (amplitude) by the difference in indices
FFT_shifted = np.concatenate((FFT[delta_index:], FFT[:delta_index]))

# # Plot the spectra
# plt.figure()
# plt.plot(Frequency, FFT, 'r', label='Original')
# plt.plot(Frequency, FFT_shifted, 'b', label='Adjusted')
# plt.legend()

# iFFT
Fid_shifted = np.fft.ifft(np.fft.fftshift(FFT_shifted))

# Define Real, Imaginary and Amplitude
Re_s = np.real(Fid_shifted)
Im_s = np.imag(Fid_shifted)
Amp_s = np.sqrt(Re_s ** 2 + Im_s ** 2)
Time_s = Time

# Plot the data with adjusted frequency
plt.figure()
plt.plot(Time_s, Re_s, 'r', label='Adjusted Re')
plt.plot(Time_s, Im_s, 'b', label='Adjusted Im')
plt.plot(Time_s, Amp_s, 'k', label='Adjusted Amp')
plt.legend()

# 5. Perform just the same thing for Glycerol
# 1. Read the data:
data_gly = np.loadtxt(file_path_glycerol)
x_gly, y_gly, z_gly = data_gly[:, 0], data_gly[:, 1], data_gly[:, 2]

Time_gly = x_gly
Real_gly = y_gly
Imag_gly = z_gly
Amp_gly = np.sqrt(Real_gly ** 2 + Imag_gly ** 2)

# 2. Crop the data from Time = 0:
if Time_gly[0] < 0:
    Time_start = 0
    Time_crop_idx = np.where(Time_gly >= Time_start)[0][0]
    Time_gly = Time_gly[Time_crop_idx:]
    Real_gly = Real_gly[Time_crop_idx:]
    Imag_gly = Imag_gly[Time_crop_idx:]


# 3. Phase the data:
delta_gly = np.zeros(360)

for phi in range(360):
    Re_gly = Real_gly * np.cos(np.deg2rad(phi)) - Imag_gly * np.sin(np.deg2rad(phi))
    Im_gly = Real_gly * np.sin(np.deg2rad(phi)) + Imag_gly * np.cos(np.deg2rad(phi))
    Ma_gly = np.sqrt(Real_gly ** 2 + Imag_gly ** 2)
    
    # Take the short time range up to ~50 microsec
    Re_cut_gly = Re_gly[:100]
    Ma_cut_gly = Ma_gly[:100]
    
    delta_gly[phi] = np.mean(Ma_cut_gly - Re_cut_gly)

idx_gly = np.argmin(delta_gly)

Re_phased_gly = Real_gly * np.cos(np.deg2rad(idx_gly)) - Imag_gly * np.sin(np.deg2rad(idx_gly))
Im_phased_gly = Real_gly * np.sin(np.deg2rad(idx_gly)) + Imag_gly * np.cos(np.deg2rad(idx_gly))
Am_phased_gly = np.sqrt(Re_phased_gly ** 2 + Im_phased_gly ** 2)

# 4. Adjust the frequency
# Calculate Frequency
numberp_gly = len(Time_gly)

dt = Time_gly[1] - Time_gly[0]
f_range = 1 / dt
f_nyquist = f_range / 2
df = 2 * (f_nyquist / numberp)
Frequency_gly = np.arange(-f_nyquist, f_nyquist + df, df)

# Create complex FID
Fid_unshifted_gly = np.array(Re_phased_gly + 1j * Im_phased_gly)

# FFT
FFT_gly = np.fft.fftshift(np.fft.fft(Fid_unshifted_gly))

# Check the length of FFT and Frequency (it is always the same, this is just in case)
if len(Frequency_gly) != len(FFT_gly):
    Frequency_gly = np.linspace(Frequency_gly[0], Frequency_gly[-1], len(FFT))

# Find index of max spectrum (amplitude)
index_max_gly = np.argmax(FFT_gly)

# Find index of zero (frequency)
array_gly = np.asarray(Frequency_gly)
index_zero_gly = (np.abs(array_gly - 0)).argmin()

# Find difference
delta_index_gly = index_max_gly - index_zero_gly

# Shift the spectra (amplitude) by the difference in indices
FFT_shifted_gly = np.concatenate((FFT_gly[delta_index_gly:], FFT_gly[:delta_index_gly]))

# iFFT
Fid_shifted_gly = np.fft.ifft(np.fft.fftshift(FFT_shifted_gly))

# Define Real, Imaginary and Amplitude
Re_r = np.real(Fid_shifted_gly)
Im_r = np.imag(Fid_shifted_gly)
Amp_r = np.sqrt(Re_r ** 2 + Im_r ** 2)
Time_r = Time_gly

# Plot the data with adjusted frequency
plt.figure()
plt.plot(Time_r, Re_r, 'r', label='Adjusted Re Glycerol')
plt.plot(Time_r, Im_r, 'b', label='Adjusted Im Glycerol')
plt.plot(Time_r, Amp_r, 'k', label='Adjusted Amp Glycerol')
plt.legend()


# 5. Crop the arrays together (they should be of the same length, but I know, I know...)
if len(Time_s) > len(Time_r):
    Time_s=Time_s[:len(Time_r)]
    Amp_s=Amp_s[:len(Time_r)]
    Re_s=Re_s[:len(Time_r)]
    Im_s=Im_s[:len(Time_r)]
else:  
    Time_r= Time_r[:len(Time_s)]
    Amp_r = Amp_r[:len(Time_s)]
    Re_r  = Re_r[:len(Time_s)]
    Im_r  = Im_r[:len(Time_s)]

Re_gly_norm = Re_s/Amp_r
Amp_gly_norm = Amp_s/Amp_r
Im_gly_norm = Im_s/Amp_r

# Plot the normalized to glycerol FID
plt.figure()
plt.plot(Time_s, Re_gly_norm, 'r-', label='Re normalized to Gly')
plt.plot(Time_s, Im_gly_norm, 'b-', label='Im normalized to Gly')
plt.plot(Time_s, Amp_gly_norm, 'k-', label='Amp normalized to Gly')
plt.legend()

# 6. Cut the ranges for fitting
array_max = np.asarray(Time_s)
maximum = (np.abs(array_max - 250)).argmin()
array_min = np.asarray(Time_s)
minimum = (np.abs(array_min - 50)).argmin()

Time_range = Time_s[minimum:maximum]
Re_gly_norm_range = Re_gly_norm[minimum:maximum]
Amp_gly_norm_range = Amp_gly_norm[minimum:maximum]
Im_gly_norm_range = Im_gly_norm[minimum:maximum]

# 7. Fit data to exponential decay
coeff = [0.9, 400, 0.1]
popt, pcov = curve_fit(decaying_exponential, Time_range, Re_gly_norm_range, p0=coeff, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
popt__, pcov__ = curve_fit(decaying_exponential, Time_range, Amp_gly_norm_range,p0=coeff, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
popt_, pcov_ = curve_fit(decaying_exponential, Time_range, Im_gly_norm_range, p0=[1, 400, 0])

# 8. Set the ranges for subtraction
Time_cropped = Time[0:maximum]
Real_cropped = Re_gly_norm[0:maximum]
Amp_cropped = Amp_gly_norm[0:maximum]
Im_cropped = Im_gly_norm[0:maximum]

# 9. Calculate the curves fitted to data within the desired range
Real_fitted = decaying_exponential(Time_cropped, *popt)
Amp_fitted = decaying_exponential(Time_cropped, *popt__)
Im_fitted = decaying_exponential(Time_cropped, *popt_)

# Plot the fitted curves and FID
plt.figure()
plt.plot(Time_cropped, Real_fitted, 'r--', label='Re Fitted curve')
plt.plot(Time_cropped, Real_cropped, 'r-', label='Im Fitted curve')
plt.plot(Time_cropped, Amp_fitted, 'k--', label='Amp Fitted curve')
plt.plot(Time_cropped, Amp_cropped, 'k-', label='Re')
plt.plot(Time_cropped, Im_fitted, 'b--', label='Im')
plt.plot(Time_cropped, Im_cropped, 'b-', label='Amp')


# 10. Subtract
Real_subtracted = Real_cropped - Real_fitted
Im_subtracted = Im_cropped - Im_fitted
Amp_subtracted = np.sqrt(Real_subtracted ** 2 + Im_subtracted ** 2)
#Im_subtracted = np.zeros((np.shape(Real_subtracted)))

# 11. Normalize
Amplitude_max_n = np.max(Amp_subtracted)
Amp_normalized = Amp_subtracted/Amplitude_max_n
Re_normalized = Real_subtracted/Amplitude_max_n
Im_normalized = Im_subtracted/Amplitude_max_n

plt.figure()
plt.plot(Time_cropped, Im_normalized, 'b-', label='Im')
plt.plot(Time_cropped, Re_normalized, 'r-', label='Re')
plt.plot(Time_cropped, Amp_normalized, 'k-', label='Am')
plt.legend()

# 12. Apodize the time-domain
coeffs = np.polyfit(Time_cropped, Amp_normalized, 1)
c = np.polyval(coeffs, Time_cropped)
d = np.argmin(np.abs(c - 1e-5))
sigma = Time_cropped[d]
apodization_function = np.exp(-(Time_cropped / sigma) ** 4)
Re_ap = Re_normalized * apodization_function
Im_ap = Im_normalized * apodization_function

# 13. Add zeros
length_diff = 16383 - len(Time_cropped)
amount_to_add = np.zeros(length_diff+1)

Re_zero = np.concatenate((Re_ap, amount_to_add))
Im_zero = np.concatenate((Im_ap, amount_to_add))

dt = Time_cropped[1] - Time_cropped[0]
Time_to_add = Time_cropped[-1] + np.arange(1, length_diff + 1) * dt

Time_FID = np.concatenate((Time_cropped, Time_to_add))
Fid_final = np.array(Re_zero + 1j * Im_zero)

# 14. Perform FFT
FFT_final = np.fft.fftshift(np.fft.fft(Fid_final))

# Frequency scale
numberp = len(Time_FID+1)
dt = Time_FID[1] - Time_FID[0]
f_range = 1 / dt
f_nyquist = f_range / 2
df = 2 * (f_nyquist / numberp)
Frequency_final = np.arange(-f_nyquist, f_nyquist + df, df)

# 15. Simple baseline
twentyperc = int(round(len(FFT_final) * 0.02))
Baseline = np.mean(np.real(FFT_final[:twentyperc]))
FFT_corrected = FFT_final - Baseline
Re_spectra = np.real(FFT_corrected)
Im_spectra = np.imag(FFT_corrected)
Amp_spectra = np.sqrt(Re_spectra ** 2 + Im_spectra ** 2)

# Plot the spectra
plt.figure()
plt.plot(Frequency_final, Re_spectra, 'r', label = 'Re')
plt.plot(Frequency_final, Im_spectra, 'b', label = 'Im')
plt.plot(Frequency_final, Amp_spectra, 'k', label = 'Amp')
plt.legend()

# 16. Perform Apodization
Maximum = np.max(np.abs(Re_spectra))
idx_max = np.argmax(np.abs(Re_spectra))
two_percent = Maximum * 0.02

b = np.argmin(np.abs(Re_spectra[idx_max:] - two_percent))
Amplitudes = np.interp(two_percent, Re_spectra, Re_spectra)
sigma_ap = Frequency_final[idx_max + b]
apodization_function_s = np.exp(-(Frequency_final / sigma_ap) ** 6)
Real_apod = Re_spectra * apodization_function_s

# 17. Calculate M2 and T2

Integral = np.trapz(np.real(Real_apod))

# Normalize FFT to the Integral value
Fur_normalized = np.real(Real_apod) / Integral

# Calculate the integral of normalized FFT to receive 1
Integral_one = np.trapz(Fur_normalized)

# Multiplication (the power ^n will give the nth moment (here it is n=2)
Multiplication = (Frequency_final ** 2) * Fur_normalized

# Plot the spectra
plt.figure()
plt.plot(Frequency_final, Multiplication, 'r', label = 'Second moment')
plt.legend()

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

print(f'M2 is: {M2}')
print(f'T2 is: {T2}')

plt.show()