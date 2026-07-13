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
# from scipy.integrate import trapezoid
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

import Calculator as Cal

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


def read_data(file_path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(file_path)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"{file_path} must contain at least three numeric columns.")
    return data[:, 0], data[:, 1], data[:, 2]

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

def voigt(x, amp, cen, wid, frac, y0):
    x = np.asarray(x)
    lorentzian = amp * (2 * wid) / (np.pi * (4 * (x - cen) ** 2 + wid**2))
    gaussian = amp * np.exp((-4 * np.log(2) * (x - cen) ** 2) / wid**2) / (
        wid * np.sqrt(np.pi / (4 * np.log(2)))
    )
    return frac * lorentzian + (1 - frac) * gaussian + y0

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

    start_data  = Cal._find_nearest(time_data, start)
    end_data    = Cal._find_nearest(time_data, stop)

    start_fid   = Cal._find_nearest(time_fid,start)
    end_fid     = Cal._find_nearest(time_fid,stop)

    mean_fid = np.mean(fid[start_fid : end_fid])
    mean_data = np.mean(data[start_data : end_data])
    return data - (mean_data - mean_fid)


def reference_long_component(time, component, end):
    time = np.asarray(time, dtype=float)
    component = np.asarray(component, dtype=float)
    minimum = Cal._find_nearest(time, end)
    time_range = time[minimum:]
    component_range = component[minimum:]
    if len(time_range) < 5:
        raise ValueError("Not enough long-component points to fit exponential reference.")
    window = min(41, len(component_range) if len(component_range) % 2 == 1 else len(component_range) - 1)
    smooth = savgol_filter(component_range, max(window, 3), 0) if window >= 3 else component_range
    popt, _ = curve_fit(Cal.decaying_exponential, time_range, smooth, p0=[5, 30, 0.5], maxfev=20000)
    return component - Cal.decaying_exponential(time, *popt)

def freq_domain_correction(time, real, imaginary=0, apodize=True, time_a=100.0, adjust=True):
    number_of_points = 2**16
    imaginary_array = np.zeros_like(np.asarray(real, dtype=float)) if np.isscalar(imaginary) else np.asarray(imaginary)
    real_array = np.asarray(real, dtype=float)

    if adjust:
        freq = Cal._calculate_frequency_scale(time)
        real_shifted, imaginary_shifted = Cal._adjust_frequency(freq, real_array, imaginary_array)
        real_td_signal = real_shifted
        imag_td_signal = imaginary_shifted

    else:
        fft = np.fft.fftshift(np.fft.fft(real_array + 1j * imaginary_array))
        real_td_signal = np.real(fft)
        imag_td_signal = np.imag(fft)

    if apodize:
        real_array, imaginary_array = apodization(time, real_td_signal, imag_td_signal, time_a)

    time_zero, fid_zero = Cal._add_zeros(time, real_td_signal, imag_td_signal, number_of_points)

    frequency = Cal._calculate_frequency_scale(time_zero)
    fft = np.fft.fftshift(np.fft.fft(fid_zero))

    _, real_baseline, _ = Cal._simple_baseline_correction(fft)

    real_apod = Cal._calculate_apodization(real_baseline, frequency)

    return frequency, real_apod


def time_domain_correction(file_path):

    time, re_original, im_original = read_data(file_path)
    re_phased, _ = Cal._time_domain_phase(re_original, im_original)
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

    m2_data, t2_data = Cal._calculate_M2(spectrum_data, freq_data)
    m2_fid, t2_fid = Cal._calculate_M2(spectrum_fid, freq_fid)

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

    start_idx = Cal._find_nearest(time, start)
    finish_idx = Cal._find_nearest(time, finish)

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

    m2_build, t2_build = Cal._calculate_M2(spectrum_build, freq_build)

    return time_build, data_build, data_fit, freq_build, spectrum_build, m2_build, t2_build


def time_range_grid(time: Sequence[float], maximum: float):
    minimum = min(time)
    start_range = np.arange(minimum, maximum, 0.5)
    finish_range = np.arange(minimum + 5, maximum + 5, 0.5)
    return start_range, finish_range
