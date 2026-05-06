import os
import re

import numpy as np

import Calculator as Cal


X_AXIS_PATTERN = r"_([0-9]+)\.dat$"


def clean_tabbed_line(line):
    while "\t\t" in line:
        line = line.replace("\t\t", "\t")
    return line.strip()


def x_axis_from_filename(file_path, fallback):
    current_file = os.path.basename(file_path)
    match = re.search(X_AXIS_PATTERN, current_file)
    if match:
        return match.group(1)

    return fallback


def read_spin_diffusion_file(file_path):
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
    time_array = np.array(time_values).flatten()
    if use_sqrt_time:
        return np.sqrt(time_array)

    return time_array


def signal_arrays(dictionary_entry):
    return {
        "short": np.array(dictionary_entry["short"]).flatten(),
        "medium": np.array(dictionary_entry["medium"]).flatten(),
        "long": np.array(dictionary_entry["long"]).flatten(),
    }


def selected_signal(signals, source):
    return signals[source]


def fit_range(time_values, signal_values, fit_from, fit_to):
    start_index = (np.abs(time_values - fit_from)).argmin()
    end_index = (np.abs(time_values - fit_to)).argmin() + 1
    fit_time = time_values[start_index:end_index]
    fit_signal = signal_values[start_index:end_index]
    return fit_time, fit_signal


def fit_spin_diffusion(time_values, signal_values, beta, r2, m2):
    fit_time_curve, fit_signal_curve, sqrt_time, r2_value = Cal.linear_fit_GS(time_values, signal_values)
    diffusion_distance = Cal.calculate_domain_size(sqrt_time, beta, r2, m2)
    return fit_time_curve, fit_signal_curve, sqrt_time, r2_value, diffusion_distance


def is_valid_fit_range(time_values, fit_time, fit_signal):
    return len(time_values) > 0 and len(fit_time) >= 2 and len(fit_signal) >= 2
