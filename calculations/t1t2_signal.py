import numpy as np

import Calculator as Cal


def clean_tabbed_line(line):
    while "\t\t" in line:
        line = line.replace("\t\t", "\t")
    return line.strip()


def add_curve(dictionary, file_path, suffix, x_axis, time_values, signal_values):
    key = file_path + suffix
    dictionary[key]["X Axis"].append(x_axis)
    dictionary[key]["Time"].extend(time_values)
    dictionary[key]["Signal"].extend(signal_values)
    return dictionary


def fit_range(time_values, signal_values, start, end, denominator):
    time_array = np.array(time_values) / denominator
    signal_array = np.array(signal_values)
    fit_time = time_array[start:end]
    fit_signal = signal_array[start:end]
    return time_array, signal_array, fit_time, fit_signal


def initial_parameters(order, signal_values, tau1, tau2, tau3):
    if len(signal_values) == 0:
        raise ValueError("selected range contains no signal points")

    first_signal = signal_values[0]
    if order == 1:
        return [first_signal, tau1, 1]
    if order == 2:
        return [first_signal, tau1, first_signal, tau2, 1]

    return [first_signal, tau1, first_signal, tau2, first_signal, tau3, 1]


def fit_relaxation(time_values, signal_values, order, initial_params):
    return Cal.fit_exponent(time_values, signal_values, order, initial_params)
