"""Numerical helpers for reconstructing FIDs from build-up FID/RecFID data.

This module contains no Qt/UI code.  It is a cleaned-up, application-local
version of the old standalone BuildUpFID ``Math.py`` pipeline.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AnalysisOptions:
    subtract_empty: bool = True
    cut_beginning: bool = True
    normalize_to_fid: bool = True
    normalize_from: float = 70.0
    normalize_to: float = 90.0
    long_component: bool = False
    long_component_from: float = 55.0
    apodize_time_domain: bool = True
    apodization_time: float = 100.0
    adjust_frequency_phase: bool = True
    adjust_fid_zero: bool = False
    fid_zero_shift: float = 0.0
    smooth: bool = False
    smooth_order: int = 1
    smooth_window: int = 5


@dataclass(frozen=True)
class SignalAnalysisResult:
    time_data: np.ndarray
    signal_data: np.ndarray
    time_fid: np.ndarray
    signal_fid: np.ndarray
    frequency_data: np.ndarray
    spectrum_data: np.ndarray
    frequency_fid: np.ndarray
    spectrum_fid: np.ndarray
    m2_data: float
    t2_data: float
    m2_fid: float
    t2_fid: float


def find_nearest(array: Sequence[float], value: float) -> int:
    values = np.asarray(array, dtype=float)
    if values.size == 0:
        raise ValueError("Cannot find a nearest value in an empty array.")
    return int(np.abs(values - value).argmin())


def gauss(x, amplitude, sigma, y0):
    return amplitude * np.exp(-(np.asarray(x) ** 2) / (2 * sigma**2)) + y0


def gauss_const_ampl(amplitude):
    return lambda x, sigma, y0: gauss(x, amplitude - y0, sigma, y0)


def polynom4(x, amplitude, c, g):
    x = np.asarray(x)
    return amplitude + c * x**2 + g * x**4


def polynom4_const_ampl(amplitude):
    return lambda x, c, g: polynom4(x, amplitude, c, g)


def polynom6(x, amplitude, c, g, h):
    x = np.asarray(x)
    return amplitude + c * x**2 + g * x**4 + h * x**6


def polynom6_const_ampl(amplitude):
    return lambda x, c, g, h: polynom6(x, amplitude, c, g, h)


def polynom8(x, amplitude, c, g, h, j):
    x = np.asarray(x)
    return amplitude + c * x**2 + g * x**4 + h * x**6 + j * x**8


def polynom8_const_ampl(amplitude):
    return lambda x, c, g, h, j: polynom8(x, amplitude, c, g, h, j)


def decaying_exponential(x, a, b, c):
    return a * np.exp(-np.asarray(x) / b) + c


def read_data(file_path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(file_path)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"{file_path} must contain at least three numeric columns.")
    return data[:, 0], data[:, 1], data[:, 2]


def calculate_frequency_scale(time: Sequence[float]) -> np.ndarray:
    time = np.asarray(time, dtype=float)
    if time.size < 2:
        raise ValueError("At least two time points are required to calculate frequency scale.")
    numberp = len(time)
    dt = time[1] - time[0]
    if dt == 0:
        raise ValueError("Time step cannot be zero.")
    f_range = 1 / dt
    f_nyquist = f_range / 2
    df = 2 * (f_nyquist / numberp)
    return np.arange(-f_nyquist, f_nyquist + df, df)[:-1]


def calculate_amplitude(real, imaginary):
    return np.sqrt(np.asarray(real) ** 2 + np.asarray(imaginary) ** 2)


def smooth_noisy_signal(signal: Sequence[float], smooth_order: int, smooth_window: int) -> np.ndarray:
    signal = np.asarray(signal, dtype=float)
    if signal.size == 0:
        return signal
    window = int(smooth_window)
    order = int(smooth_order)
    window = max(window, order + 2)
    if window % 2 == 0:
        window += 1
    if window > signal.size:
        window = signal.size if signal.size % 2 == 1 else signal.size - 1
    if window <= order or window < 3:
        logger.warning("Skipping smoothing because window/order are invalid for signal length.")
        return signal
    return savgol_filter(signal, window, order)


def adjust_phase(real, imaginary):
    real = np.asarray(real, dtype=float)
    imaginary = np.asarray(imaginary, dtype=float)
    delta = np.zeros(360)
    for phi in range(360):
        re_phased = real * np.cos(np.deg2rad(phi)) - imaginary * np.sin(np.deg2rad(phi))
        im_phased = real * np.sin(np.deg2rad(phi)) + imaginary * np.cos(np.deg2rad(phi))
        magnitude_phased = calculate_amplitude(re_phased, im_phased)
        delta[phi] = np.mean(magnitude_phased[:5] - re_phased[:5])
    idx = int(np.argmin(delta))
    re = real * np.cos(np.deg2rad(idx)) - imaginary * np.sin(np.deg2rad(idx))
    im = real * np.sin(np.deg2rad(idx)) + imaginary * np.cos(np.deg2rad(idx))
    return re, im


def voigt(x, amp, cen, wid, frac, y0):
    x = np.asarray(x)
    lorentzian = amp * (2 * wid) / (np.pi * (4 * (x - cen) ** 2 + wid**2))
    gaussian = amp * np.exp((-4 * np.log(2) * (x - cen) ** 2) / wid**2) / (
        wid * np.sqrt(np.pi / (4 * np.log(2)))
    )
    return frac * lorentzian + (1 - frac) * gaussian + y0


def fit_fft_with_voigh(frequency, fft):
    frequency = np.asarray(frequency, dtype=float)
    spectrum = np.abs(fft)
    idx_peak = int(np.argmax(spectrum))
    mu0 = frequency[idx_peak]
    window = (frequency[-1] - frequency[0]) * 0.1
    mask = (frequency >= (mu0 - window / 2)) & (frequency <= (mu0 + window / 2))
    x = frequency[mask]
    y = spectrum[mask]
    if x.size < 5:
        raise ValueError("Not enough points to fit FFT with Voigt profile.")
    p0 = [np.max(spectrum), frequency[idx_peak], (frequency[-1] - frequency[0]) / 50, 0.5, 0]
    lower = [0, frequency[0], frequency[1] - frequency[0], 0, -10]
    upper = [np.inf, frequency[-1], frequency[-1] - frequency[0], 1, 10]
    popt, _ = curve_fit(voigt, x, y, bounds=(lower, upper), p0=p0, maxfev=20000)
    return popt


def adjust_frequency(frequency, real, imaginary):
    fid_unshifted = np.asarray(real) + 1j * np.asarray(imaginary)
    fft = np.fft.fftshift(np.fft.fft(fid_unshifted))
    frequency = np.asarray(frequency, dtype=float)
    if len(frequency) != len(fft):
        frequency = np.linspace(frequency[0], frequency[-1], len(fft))
    try:
        popt_fitting = fit_fft_with_voigh(frequency, fft)
        index_max = int(np.argmax(voigt(frequency, *popt_fitting)))
    except Exception:
        logger.exception("FFT Voigt fitting failed; using simple maximum for frequency adjustment.")
        index_max = int(np.argmax(np.abs(fft)))
    index_zero = find_nearest(frequency, 0)
    delta_index = index_max - index_zero
    if delta_index == 0:
        logger.info("Frequency was not adjusted; maximum is already at zero.")
        return np.asarray(real), np.asarray(imaginary)
    fft_shifted = np.concatenate((fft[delta_index:], fft[:delta_index]))
    fid_shifted = np.fft.ifft(np.fft.ifftshift(fft_shifted))
    return np.real(fid_shifted), np.imag(fid_shifted)


def simple_baseline_correction(fft):
    fft = np.asarray(fft)
    n_baseline = max(1, int(round(len(fft) * 0.02)))
    baseline = np.mean(np.real(fft[:n_baseline]))
    return np.real(fft - baseline)


def apodization_fft(real, freq):
    real = np.asarray(real, dtype=float)
    freq = np.asarray(freq, dtype=float)
    maximum = np.max(np.abs(real))
    if maximum == 0:
        return real
    idx_max = int(np.argmax(np.abs(real)))
    target = maximum * 0.03
    relative_idx = int(np.argmin(np.abs(real[idx_max:] - target))) if idx_max < len(real) else 0
    sigma_ap = abs(freq[min(idx_max + relative_idx, len(freq) - 1)])
    if sigma_ap == 0:
        logger.warning("Skipping FFT apodization because sigma is zero.")
        return real
    return real * np.exp(-(freq / sigma_ap) ** 2)


def add_zeros(time, real, imaginary, number_of_points):
    time = np.asarray(time, dtype=float)
    real = np.asarray(real, dtype=float)
    imaginary = np.asarray(imaginary, dtype=float)
    length_diff = int(number_of_points) - len(time)
    if length_diff <= 0:
        return time, real + 1j * imaginary
    amount_to_add = np.zeros(length_diff + 1)
    re_zero = np.concatenate((real, amount_to_add))
    im_zero = np.concatenate((np.zeros(len(real)), amount_to_add))
    dt = time[1] - time[0]
    time_to_add = time[-1] + np.arange(1, length_diff + 1) * dt
    time_zero = np.concatenate((time, time_to_add))
    return time_zero, (re_zero + 1j * im_zero)[:-1]


def apodization(time, real, imaginary, sigma):
    time = np.asarray(time, dtype=float)
    if sigma == 0:
        return np.asarray(real), np.asarray(imaginary)
    apodization_function = np.exp(-(time / sigma) ** 4)
    return np.asarray(real) * apodization_function, np.asarray(imaginary) * apodization_function


def cut_beginning(time, data):
    data = np.asarray(data, dtype=float)
    time = np.asarray(time, dtype=float)
    if data.size == 0:
        return time, data
    idx = int(np.argmax(data))
    return time[idx:], data[idx:]


def normalize_to_fid(fid, data, time_fid, time_data, start, stop):
    fid = np.asarray(fid, dtype=float)
    data = np.asarray(data, dtype=float)
    mean_fid = np.mean(fid[find_nearest(time_fid, start) : find_nearest(time_fid, stop)])
    mean_data = np.mean(data[find_nearest(time_data, start) : find_nearest(time_data, stop)])
    return data - (mean_data - mean_fid)


def reference_long_component(time, component, end):
    time = np.asarray(time, dtype=float)
    component = np.asarray(component, dtype=float)
    minimum = find_nearest(time, end)
    time_range = time[minimum:]
    component_range = component[minimum:]
    if len(time_range) < 5:
        raise ValueError("Not enough long-component points to fit exponential reference.")
    window = min(41, len(component_range) if len(component_range) % 2 == 1 else len(component_range) - 1)
    smooth = savgol_filter(component_range, max(window, 3), 0) if window >= 3 else component_range
    popt, _ = curve_fit(decaying_exponential, time_range, smooth, p0=[5, 30, 0.5], maxfev=20000)
    return component - decaying_exponential(time, *popt)


def calculate_M2(fft_real, frequency):
    fft_real = np.asarray(fft_real, dtype=float)
    frequency = np.asarray(frequency, dtype=float)
    integral = np.trapezoid(np.real(fft_real))
    if integral == 0:
        raise ValueError("Cannot calculate M2 because FFT integral is zero.")
    fur_normalized = np.real(fft_real) / integral
    multiplication = (frequency**2) * fur_normalized
    m2 = np.trapezoid(multiplication) * 4 * np.pi**2
    if np.abs(np.mean(multiplication[:10])) > 10**-6:
        logger.warning("FFT edge moment is high; apodization may be insufficient.")
    if m2 < 0:
        return 0.0, 0.0
    return float(m2), float(np.sqrt(2 / m2)) if m2 else 0.0


def freq_domain_correction(time, real, imaginary=0, apodize=True, time_a=100.0, adjust=True):
    number_of_points = 2**16
    imaginary_array = np.zeros_like(np.asarray(real, dtype=float)) if np.isscalar(imaginary) else np.asarray(imaginary)
    real_array = np.asarray(real, dtype=float)
    if apodize:
        real_array, imaginary_array = apodization(time, real_array, imaginary_array, time_a)
    if adjust:
        freq = calculate_frequency_scale(time)
        re_ph, im_ph = adjust_phase(real_array, imaginary_array)
        real_fft, imaginary_fft = adjust_frequency(freq, re_ph, im_ph)
    else:
        fft = np.fft.fftshift(np.fft.fft(real_array + 1j * imaginary_array))
        real_fft = np.real(fft)
        imaginary_fft = np.imag(fft)
    time_zero, fid_zero = add_zeros(time, real_fft, imaginary_fft, number_of_points)
    frequency = calculate_frequency_scale(time_zero)
    fft = np.fft.fftshift(np.fft.fft(fid_zero))
    return frequency, apodization_fft(simple_baseline_correction(fft), frequency)


def time_domain_correction(file_path):
    time, re_original, im_original = read_data(file_path)
    re_phased, _ = adjust_phase(re_original, im_original)
    return time, re_phased


def nmr_signal_correction(
    file_path,
    file_path_fid,
    file_path_empty=None,
    file_path_empty_fid=None,
    options: AnalysisOptions | None = None,
):
    options = options or AnalysisOptions()
    time_td, re_td = time_domain_correction(file_path)
    time_td_fid, re_td_fid = time_domain_correction(file_path_fid)
    if options.smooth:
        re_td = smooth_noisy_signal(re_td, options.smooth_order, options.smooth_window)
        re_td_fid = smooth_noisy_signal(re_td_fid, options.smooth_order, options.smooth_window)
    if options.subtract_empty and file_path_empty and file_path_empty_fid:
        _, re_empty, _ = read_data(file_path_empty)
        re_td = re_td - re_empty[: len(time_td)]
        _, re_empty_fid, _ = read_data(file_path_empty_fid)
        re_td_fid = re_td_fid - re_empty_fid[: len(time_td_fid)]
    if options.adjust_fid_zero:
        time_td_fid = time_td_fid + options.fid_zero_shift
    if options.cut_beginning:
        _, re_td = cut_beginning(time_td, re_td)
        time_td = np.linspace(0, 0 + (0.5 * len(re_td)), len(re_td), endpoint=False)
        time_td_fid, re_td_fid = cut_beginning(time_td_fid, re_td_fid)
    if options.normalize_to_fid:
        re_td = normalize_to_fid(
            re_td_fid, re_td, time_td_fid, time_td, options.normalize_from, options.normalize_to
        )
    if options.long_component:
        re_td = reference_long_component(time_td, re_td, options.long_component_from)
        re_td_fid = reference_long_component(time_td_fid, re_td_fid, options.long_component_from)
    return time_td, re_td, time_td_fid, re_td_fid


def for_the_sake_of_beauty(time_td, re_td, time_td_fid, re_td_fid, apodize=True, time_a=100.0):
    if apodize:
        re_td, _ = apodization(time_td, re_td, 0, time_a)
        re_td_fid, _ = apodization(time_td_fid, re_td_fid, 0, time_a)
    return time_td, re_td, time_td_fid, re_td_fid


def analyze_signal(data_file, fid_file, data_empty=None, fid_empty=None, options: AnalysisOptions | None = None):
    options = options or AnalysisOptions()
    time_td, re_td, time_fid, re_fid = nmr_signal_correction(
        data_file, fid_file, data_empty, fid_empty, options
    )
    time_td, re_td, time_fid, re_fid = for_the_sake_of_beauty(
        time_td, re_td, time_fid, re_fid, options.apodize_time_domain, options.apodization_time
    )
    freq_data, spectrum_data = freq_domain_correction(
        time_td, re_td, 0, options.apodize_time_domain, options.apodization_time, options.adjust_frequency_phase
    )
    freq_fid, spectrum_fid = freq_domain_correction(
        time_fid, re_fid, 0, options.apodize_time_domain, options.apodization_time, options.adjust_frequency_phase
    )
    m2_data, t2_data = calculate_M2(spectrum_data, freq_data)
    m2_fid, t2_fid = calculate_M2(spectrum_fid, freq_fid)
    return SignalAnalysisResult(
        time_td,
        re_td,
        time_fid,
        re_fid,
        freq_data,
        spectrum_data,
        freq_fid,
        spectrum_fid,
        round(m2_data, 5),
        round(t2_data, 5),
        round(m2_fid, 5),
        round(t2_fid, 5),
    )


def extract_echo_time(file_path: str | Path) -> float:
    match = re.search(r".*_\s*(\d+(?:\.\d+)?)_c\.dat$", str(file_path))
    if not match:
        raise ValueError("Echo time could not be parsed from file name.")
    return float(match.group(1))


def find_maximum_se(echo_time, maximum, start=0, end=0):
    echo_time = np.asarray(echo_time, dtype=float)
    maximum = np.asarray(maximum, dtype=float)
    if echo_time.size != maximum.size or echo_time.size < 3:
        raise ValueError("At least three echo-time/maximum pairs are required.")
    echo_time_to_fit = np.arange(0, echo_time[-1], 0.01)
    end_idx = None if int(end) == 0 else int(-end)
    echo_time_cut = echo_time[int(start) : end_idx]
    maximum_cut = maximum[int(start) : end_idx]
    popt, _ = curve_fit(gauss, echo_time_cut, maximum_cut, p0=[10, 6, 1], maxfev=20000)
    return float(gauss(0, *popt)), gauss(echo_time_to_fit, *popt), echo_time_to_fit


def build_up_fid(time, data, amplitude, function_to_fit, start, finish):
    time = np.asarray(time, dtype=float)
    data = np.asarray(data, dtype=float)
    start_idx = find_nearest(time, start)
    finish_idx = find_nearest(time, finish)
    if finish_idx <= start_idx + 1:
        raise ValueError("Finish time must be after begin time and include at least two data points.")
    time_cut = time[start_idx:finish_idx]
    data_cut = data[start_idx:finish_idx]
    delta_time = time_cut[1] - time_cut[0]
    time_build_from_zero = np.arange(0, start, delta_time)
    if function_to_fit == "Polynom 4":
        popt, _ = curve_fit(polynom4_const_ampl(amplitude), time_cut, data_cut, p0=[0.005, 0.005])
        data_built = polynom4(time, amplitude, *popt)
        data_build_from_zero = polynom4(time_build_from_zero, amplitude, *popt)
    elif function_to_fit == "Polynom 6":
        popt, _ = curve_fit(polynom6_const_ampl(amplitude), time_cut, data_cut, p0=[0.005] * 3)
        data_built = polynom6(time, amplitude, *popt)
        data_build_from_zero = polynom6(time_build_from_zero, amplitude, *popt)
    elif function_to_fit == "Polynom 8":
        popt, _ = curve_fit(polynom8_const_ampl(amplitude), time_cut, data_cut, p0=[0.005] * 4)
        data_built = polynom8(time, amplitude, *popt)
        data_build_from_zero = polynom8(time_build_from_zero, amplitude, *popt)
    elif function_to_fit == "Gaussian":
        popt, _ = curve_fit(gauss_const_ampl(amplitude), time_cut, data_cut, p0=[8, 0])
        data_built = gauss(time, amplitude - popt[1], *popt)
        data_build_from_zero = gauss(time_build_from_zero, amplitude - popt[1], *popt)
    else:
        raise ValueError(f"Unsupported build-up function: {function_to_fit}")
    time_build_middle = time[start_idx + 1 : finish_idx]
    length = len(time_build_middle)
    weight = np.linspace(0, 1, length)
    data_build_middle = weight * data[start_idx + 1 : finish_idx] + (1 - weight) * data_built[start_idx + 1 : finish_idx]
    time_build_end = time[finish_idx + 1 :]
    data_build_end = data[finish_idx + 1 :]
    return (
        np.concatenate((time_build_from_zero, time_build_middle, time_build_end)),
        np.concatenate((data_build_from_zero, data_build_middle, data_build_end)),
        data_built,
    )


def analyze_build_up(time, data, extrapolation, function_to_fit, begin, finish, apodization_time=100.0):
    time_build, data_build, data_fit = build_up_fid(time, data, extrapolation, function_to_fit, begin, finish)
    freq_build, spectrum_build = freq_domain_correction(
        time_build, data_build, 0, True, apodization_time, True
    )
    m2_build, t2_build = calculate_M2(spectrum_build, freq_build)
    return time_build, data_build, data_fit, freq_build, spectrum_build, m2_build, t2_build


def time_range_grid(time: Sequence[float], maximum: float):
    minimum = min(time)
    start_range = np.arange(minimum, maximum, 0.5)
    finish_range = np.arange(minimum + 5, maximum + 5, 0.5)
    return start_range, finish_range
