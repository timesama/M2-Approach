import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit

import Calculator as Cal


DEFAULT_FILTER_FROM = 0
DEFAULT_FILTER_TO = 20


def calculate_linearization(dq_time, t2, dq, time_min, time_max):
    """Build DQ linearization arrays for the selected DQ time window."""
    selected_mask = (dq_time >= time_min) & (dq_time <= time_max)
    selected_time = dq_time[selected_mask]

    if len(selected_time) < 3:
        time_min = DEFAULT_FILTER_FROM
        time_max = DEFAULT_FILTER_TO
        selected_mask = (dq_time >= time_min) & (dq_time <= time_max)
        selected_time = dq_time[selected_mask]

    selected_t2 = t2[selected_mask]
    if len(selected_time) <= 1 or len(selected_t2) <= 1:
        return None

    coefficients = np.polyfit(selected_time, selected_t2, 1)
    integral = trapezoid(dq)
    dq_norm = dq / integral
    t2_linearized = coefficients[0] * dq_time + coefficients[1]
    x_line = np.arange(0, 105.1, 0.1)
    y_line = np.polyval(coefficients, x_line)

    return {
        "time_min": time_min,
        "time_max": time_max,
        "t2_linearized": t2_linearized,
        "dq_norm": dq_norm,
        "x_line": x_line,
        "y_line": y_line,
    }


def prepare_t2_axis(t2_linearized, use_log_scale):
    """Return the T2 axis values and label for linear/log plotting."""
    if use_log_scale:
        return np.log10(t2_linearized), "log(T₂*)"

    return t2_linearized, "T₂*"


def fit_distribution(x_values, y_values, fit_function, use_log_scale):
    """Fit the selected DQ distribution model without touching the UI."""
    if use_log_scale:
        x = np.log10(x_values)
        initial_params = [1, 1, 1, 0]
        fit_bounds = ([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1, np.inf])
        axis_label = "log(T₂*)"
    else:
        x = x_values
        initial_params = [1, 10, 10, 0]
        fit_bounds = ([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1, np.inf])
        axis_label = "T₂*"

    try:
        x_fit = np.arange(0, np.max(x) + 0.001, 0.01)
    except Exception:
        x_fit = np.arange(0, 100 + 0.001, 0.01)

    gaussian_lorenz_bounds = ([0, 0, 0, -10], [np.inf, np.inf, np.inf, np.inf])
    if fit_function == "Gauss":
        params, _ = curve_fit(Cal.gaussian, x, y_values, p0=initial_params, bounds=gaussian_lorenz_bounds)
        y_fit = Cal.gaussian(x_fit, *params)
        y_r2 = Cal.gaussian(x, *params)
        center = params[1]
        fwhm = params[2]
        lorenz_fraction = 0
    elif fit_function == "Lorenz":
        params, _ = curve_fit(Cal.lorenz, x, y_values, p0=initial_params, bounds=gaussian_lorenz_bounds)
        y_fit = Cal.lorenz(x_fit, *params)
        y_r2 = Cal.lorenz(x, *params)
        center = params[1]
        fwhm = params[2]
        lorenz_fraction = 1
    elif fit_function == "Pseudo Voigt":
        params, _ = curve_fit(Cal.voigt, x, y_values, bounds=fit_bounds)
        y_fit = Cal.voigt(x_fit, *params)
        y_r2 = Cal.voigt(x, *params)
        center = params[1]
        fwhm = params[2]
        lorenz_fraction = params[3]
    else:
        return None

    return {
        "x_fit": x_fit,
        "y_fit": y_fit,
        "r_squared": Cal.calculate_r_squared(y_values, y_r2),
        "center": center,
        "fwhm": fwhm,
        "lorenz_fraction": lorenz_fraction,
        "axis_label": axis_label,
    }
