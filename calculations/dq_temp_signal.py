import numpy as np
from scipy.optimize import curve_fit

import Calculator as Cal

INITIAL_PARAMS = [1, 5, 5, 0]
VOIGT_BOUNDS = ([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, 1, np.inf])
GAUSS_LORENZ_BOUNDS = ([0, 0, 0, -10], [np.inf, np.inf, np.inf, np.inf])


def load_distribution_file(file_path):
    """Load saved DQ distribution columns needed by temperature comparison."""
    data = np.loadtxt(file_path, delimiter=",")
    if data.ndim != 2 or data.shape[1] < 6:
        raise ValueError("file must contain at least 6 columns")

    return data[:, [4, 5]]


def fit_distribution(data):
    """Fit Gauss/Lorenz/Voigt summaries for one DQ temperature file."""
    t2_linearized = data[:, 0]
    dq_normalized = data[:, 1]
    t2_fit = np.arange(0, np.max(t2_linearized) + 0.001, 0.01)

    gauss_params, _ = curve_fit(
        Cal.gaussian,
        t2_linearized,
        dq_normalized,
        p0=INITIAL_PARAMS,
        bounds=GAUSS_LORENZ_BOUNDS,
    )
    lorenz_params, _ = curve_fit(
        Cal.lorenz,
        t2_linearized,
        dq_normalized,
        p0=INITIAL_PARAMS,
        bounds=GAUSS_LORENZ_BOUNDS,
    )
    voigt_params, _ = curve_fit(
        Cal.voigt,
        t2_linearized,
        dq_normalized,
        bounds=VOIGT_BOUNDS,
    )

    voigt_fit = Cal.voigt(t2_fit, *voigt_params)
    derivative_x, derivative_y, derivative_center = Cal.derivative_peak_find(t2_linearized, dq_normalized)

    return {
        "t2_linearized": t2_linearized,
        "dq_normalized": dq_normalized,
        "t2_fit": t2_fit,
        "voigt_fit": voigt_fit,
        "derivative_x": derivative_x,
        "derivative_y": derivative_y,
        "center_gauss": gauss_params[1],
        "center_lorenz": lorenz_params[1],
        "center_voigt": voigt_params[1],
        "center_derivative": derivative_center,
        "fwhm_gauss": gauss_params[2],
        "fwhm_lorenz": lorenz_params[2],
        "fwhm_voigt": voigt_params[2],
    }


def finite_xy(x_values, y_values):
    """Return finite x/y pairs for robust summary plotting."""
    x_array = np.asarray(x_values, dtype=float)
    y_array = np.asarray(y_values, dtype=float)
    mask = np.isfinite(x_array) & np.isfinite(y_array)
    return x_array[mask], y_array[mask]
