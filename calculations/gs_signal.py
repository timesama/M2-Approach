import os
import re

import numpy as np
from scipy.optimize import curve_fit

import Calculator as Cal


X_AXIS_PATTERN = r"_([0-9]+)\.dat$"


def clean_tabbed_line(line):
    """Collapse repeated tabs from Spin Diffusion text rows."""
    while "\t\t" in line:
        line = line.replace("\t\t", "\t")
    return line.strip()


def x_axis_from_filename(file_path, fallback):
    """Extract the x-axis value from a GS filename or use a fallback row."""
    current_file = os.path.basename(file_path)
    match = re.search(X_AXIS_PATTERN, current_file)
    if match:
        return match.group(1)

    return fallback


def read_spin_diffusion_file(file_path):
    """Read time and short/medium/long signals from a GS data file."""
    sqrt_time = []
    short_signal = []
    medium_signal = []
    long_signal = []

    with open(file_path, "r") as data:
        lines = [clean_tabbed_line(line.rstrip("\n")) for line in data if line.strip()]

    for line in lines[1:]:
        parts = line.split("\t")
        sqrt_time.append(float(parts[0]))
        short_signal.append(float(parts[4]))
        medium_signal.append(float(parts[5]))
        long_signal.append(float(parts[6]))

    return sqrt_time, short_signal, medium_signal, long_signal


def transformed_time(time_values, use_sqrt_time):
    """Return the raw or square-root transformed GS time axis."""
    time_array = np.array(time_values).flatten()
    if use_sqrt_time:
        return np.sqrt(time_array)

    return time_array


def signal_arrays(dictionary_entry):
    """Convert stored GS signal lists into NumPy arrays by source."""
    return {
        "short": np.array(dictionary_entry["short"]).flatten(),
        "medium": np.array(dictionary_entry["medium"]).flatten(),
        "long": np.array(dictionary_entry["long"]).flatten(),
    }


def selected_signal(signals, source):
    """Return the signal array selected by the GS source radio buttons."""
    return signals[source]


def fit_range(time_values, signal_values, fit_from, fit_to):
    """Slice GS data to the nearest requested fit bounds."""
    start_index = (np.abs(time_values - fit_from)).argmin()
    end_index = (np.abs(time_values - fit_to)).argmin() + 1
    fit_time = time_values[start_index:end_index]
    fit_signal = signal_values[start_index:end_index]
    return fit_time, fit_signal


def fit_spin_diffusion(time_values, signal_values, time_plateau, signal_plateau, beta, r2, m2):
    """Fit the GS line and calculate the domain size from it."""
    fit_time_curve, fit_signal_curve, fit_time_plateau, fit_signal_plateau, sqrt_time, r2_value = linear_fit_GS(time_values, signal_values, time_plateau, signal_plateau)
    diffusion_distance = calculate_domain_size(sqrt_time, beta, r2, m2)
    return fit_time_curve, fit_signal_curve, fit_time_plateau, fit_signal_plateau, sqrt_time, r2_value, diffusion_distance


def is_valid_fit_range(time_values, fit_time, fit_signal):
    """Check that the selected GS fit range has enough points."""
    return len(time_values) > 0 and len(fit_time) >= 2 and len(fit_signal) >= 2


def range_is_valid(selected_range, allowed_minimum, data_max, time_values):
    selected_from, selected_to = selected_range

    if not (
        allowed_minimum <= selected_from <= data_max
        and allowed_minimum <= selected_to <= data_max
    ):
        return False

    # The nearest indices must be ordered and distinct.
    return Cal._find_nearest(time_values, selected_from) < Cal._find_nearest(time_values, selected_to)

def calculate_domain_size(t05, beta, r2, M2):
    pi = np.pi
    Dsd = np.sqrt(pi * M2) * ((r2*10**(-10))**2)/6
    d = (2 * beta * t05 * np.sqrt(Dsd) / np.sqrt(pi))*10**(9)
    d = np.round(d, 3)
    return d

def linear_fit_GS(Time, Signal, Time_plateau, Signal_plateau):
    # Fit using curve_fit
    Time = np.array(Time).flatten()
    Signal = np.array(Signal).flatten()
    popt, _ = curve_fit(Cal.linear_model, Time, Signal)
    a_1, b_1 = popt

    # Fit plateau
    Time_plateau = np.array(Time_plateau).flatten()
    Signal_plateau = np.array(Signal_plateau).flatten()
    popt_plateau, _ = curve_fit(Cal.linear_model, Time_plateau, Signal_plateau)
    a_2, b_2 = popt_plateau

    # Extrapolate to find Time (sqrtT) where Signal = 0: a * sqrtT + b = 0 -> sqrtT = -b / a
    # sqrtT = -b / a if a != 0 else np.nan
    diff_a = a_1 - a_2
    diff_b = b_2 - b_1

    if diff_a == 0:
        sqrtT = np.nan
    else:
        sqrtT = diff_b/diff_a

    # Generate Time_fit and fitted_curve

    Time_fit = np.linspace(min(Time), sqrtT + 50, 300)
    fitted_curve = Cal.linear_model(Time_fit, a_1, b_1)

    Time_fit_plateau = np.linspace(sqrtT-100, max(Time_plateau), 300)
    fitted_curve_plateau = Cal.linear_model(Time_fit_plateau, a_2, b_2)

    # Compute R^2 manually
    residuals = Signal - Cal.linear_model(Time, a_1, b_1)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((Signal - np.mean(Signal))**2)
    R2 = 1 - (ss_res / ss_tot)

    R2 = np.round(R2, 4)

    return Time_fit, fitted_curve, Time_fit_plateau, fitted_curve_plateau, sqrtT, R2