"""Numerical Dres fitting helpers for the DQMQ tab.

Adapted from the standalone ``DQMQ.py`` reference script without any UI or
matplotlib dependencies.
"""

import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit

K = 0.4
D_GRID = np.linspace(0, 0.6, 20000)
VALID_KERNELS = ["gaussian", "abragam", "pake", "weibull", "a-l", "p-l"]


def normalize_distribution(p_values, d_values):
    area = trapezoid(p_values, d_values)
    if area <= 0 or not np.isfinite(area):
        return np.zeros_like(p_values)
    return p_values / area


def p_gaussian(d_values, mu, sigma):
    sigma = max(float(sigma), 1e-9)
    squared_offset = (d_values - mu) ** 2
    variance_scale = 2 * sigma**2
    p_values = np.exp(-squared_offset / variance_scale)
    return normalize_distribution(p_values, d_values)


def p_gaussian_2d(d_values, mu1, sigma1, mu2, sigma2, frac1):
    frac1 = np.clip(frac1, 0.0, 1.0)
    p1 = p_gaussian(d_values, mu1, sigma1)
    p2 = p_gaussian(d_values, mu2, sigma2)
    return frac1 * p1 + (1.0 - frac1) * p2


def dq_kernel(x_values, kernel, beta=2.0, k_value=K):
    kernel = kernel.lower()
    beta = max(float(beta), 1e-12)

    if kernel == "gaussian":
        return 1.0 - np.exp(-k_value * x_values**2)
    if kernel == "abragam":
        return 1.0 - np.exp(-k_value * x_values**2) * np.sinc(x_values / np.pi)
    if kernel == "pake":
        return 1.0 - np.exp(-k_value * x_values**2) * np.cos(x_values)
    if kernel == "weibull":
        return 1.0 - np.exp(-k_value * x_values**beta)
    if kernel == "a-l":
        return 1.0 - np.exp(-(k_value * x_values)**beta) * np.sinc(x_values/ np.pi)
    if kernel == "p-l":
        return 1.0 - np.exp(-(k_value * x_values)**beta) * np.cos(k_value * x_values)

    raise ValueError(f"Unknown kernel: {kernel}. Use one of {VALID_KERNELS}")


def _integrated_ndq_response(time_values, p_values, kernel, beta, k_value=K):
    integrated_values = []
    for tau_i in np.asarray(time_values, dtype=float):
        scaled_dres = D_GRID * tau_i
        kernel_response = dq_kernel(scaled_dres, kernel, beta, k_value)
        integrated_response = trapezoid(p_values * kernel_response, D_GRID)
        integrated_values.append(0.5 * integrated_response)

    return np.array(integrated_values)


def ndq_1d(time_values, mu, sigma, beta, kernel, k_value=K):
    p_values = p_gaussian(D_GRID, mu, sigma)
    return _integrated_ndq_response(time_values, p_values, kernel, beta, k_value)


def ndq_2d(time_values, mu1, sigma1, mu2, sigma2, frac1, beta, kernel, k_value=K):
    p_values = p_gaussian_2d(D_GRID, mu1, sigma1, mu2, sigma2, frac1)
    return _integrated_ndq_response(time_values, p_values, kernel, beta, k_value)


def ndq_single(time_values, dres, k_value=K):
    time_array = np.asarray(time_values, dtype=float)
    exponent = -k_value * dres**2 * time_array**2
    return 0.5 * (1.0 - np.exp(exponent))


def make_fit_model(kernel, n_components, k_value=K):
    kernel = kernel.lower()
    if kernel not in VALID_KERNELS:
        raise ValueError(f"kernel must be one of {VALID_KERNELS}")

    if n_components == 1:

        def model(time_values, mu, sigma, beta):
            return ndq_1d(time_values, mu, sigma, beta, kernel, k_value)

        return model

    if n_components == 2:

        def model(time_values, mu1, sigma1, mu2, sigma2, frac1, beta):
            return ndq_2d(
                time_values, mu1, sigma1, mu2, sigma2, frac1, beta, kernel, k_value
            )

        return model

    raise ValueError("n_components must be 1 or 2")


def fit_selected_model(
    fullTimearray, time0, ndq0, kernel="gaussian", n_components=1, p0=None, k_value=K
):
    kernel = kernel.lower()
    if kernel not in VALID_KERNELS:
        raise ValueError(f"kernel must be one of {VALID_KERNELS}")

    k_value = float(k_value)
    if not np.isfinite(k_value) or k_value <= 0:
        raise ValueError("Dres K value must be a positive finite number.")

    model = make_fit_model(kernel, n_components, k_value)
    if n_components == 1:
        default_p0 = [0.25, 1e-3, 2.0]
        bounds_min = [0, 0, 0]
        bounds_max = [0.500, 1.0, 6.0]
        param_names = ["mu", "sigma", "beta"]
    elif n_components == 2:
        default_p0 = [0.25, 0.001, 0.05, 0.001, 0.5, 2.0]
        bounds_min = [0, 0, 0, 0, 0.0, 0]
        bounds_max = [0.5000, 1.0, 0.5000, 1.0, 1.0, 6.0]
        param_names = ["mu1", "sigma1", "mu2", "sigma2", "frac1", "beta"]
    else:
        raise ValueError("n_components must be 1 or 2")

    if p0 is None:
        p0 = default_p0
    if len(p0) != len(default_p0):
        raise ValueError(
            f"Expected {len(default_p0)} initial parameters for {n_components} Dres component(s)."
        )

    p0 = np.asarray(p0, dtype=float)
    if not np.all(np.isfinite(p0)):
        raise ValueError("Dres initial parameters must be finite numbers.")

    bounds_min = np.asarray(bounds_min, dtype=float)
    bounds_max = np.asarray(bounds_max, dtype=float)
    below_bounds = p0 < bounds_min
    above_bounds = p0 > bounds_max
    if np.any(below_bounds) or np.any(above_bounds):
        raise ValueError(
            "Dres initial parameters are outside the allowed fitting bounds."
        )

    popt, pcov = curve_fit(
        model,
        time0,
        ndq0,
        p0=p0,
        bounds=(bounds_min, bounds_max),
        maxfev=50000,
    )

    fit_x = np.arange(0, fullTimearray[-1], 0.1)
    fit_y = model(fit_x, *popt)

    return {
        "kernel": kernel,
        "n_components": n_components,
        "popt": popt,
        "pcov": pcov,
        "fit_x": fit_x,
        "fit_y": fit_y,
        "param_names": param_names,
        "k_value": k_value,
    }


def build_distribution(fit_result):
    n_components = fit_result["n_components"]
    popt = fit_result["popt"]

    if n_components == 1:
        p_values = p_gaussian(D_GRID, popt[0], popt[1])
    elif n_components == 2:
        p_values = p_gaussian_2d(D_GRID, popt[0], popt[1], popt[2], popt[3], popt[4])
    else:
        raise ValueError("n_components must be 1 or 2")

    d_plot = D_GRID / (2 * np.pi) * 1000.0
    return d_plot, p_values


def build_singular_distribution(center, sigma):
    p_values = p_gaussian(D_GRID, center, sigma)

    d_plot = D_GRID / (2 * np.pi) * 1000.0
    return d_plot, p_values
