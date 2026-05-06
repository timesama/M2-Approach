"""Numerical Dres fitting helpers for the DQMQ tab.

Adapted from the standalone ``DQMQ.py`` reference script without any UI or
matplotlib dependencies.
"""

import numpy as np
from scipy.optimize import curve_fit

K = 0.4
D_GRID = np.linspace(0, 0.94, 60000)
VALID_KERNELS = ["gaussian", "abragam", "pake", "weibull", "a-l"]


def normalize_distribution(p_values, d_values):
    area = np.trapz(p_values, d_values)
    if area <= 0 or not np.isfinite(area):
        return np.zeros_like(p_values)
    return p_values / area


def p_gaussian(d_values, mu, sigma):
    sigma = max(float(sigma), 1e-9)
    p_values = np.exp(-((d_values - mu) ** 2) / (2 * sigma**2))
    return normalize_distribution(p_values, d_values)


def p_gaussian_2d(d_values, mu1, sigma1, mu2, sigma2, frac1):
    frac1 = np.clip(frac1, 0.0, 1.0)
    p1 = p_gaussian(d_values, mu1, sigma1)
    p2 = p_gaussian(d_values, mu2, sigma2)
    return frac1 * p1 + (1.0 - frac1) * p2


def dq_kernel(x_values, kernel, beta=2.0):
    kernel = kernel.lower()
    beta = max(float(beta), 1e-12)

    if kernel == "gaussian":
        return 1.0 - np.exp(-K * x_values**2)
    if kernel == "abragam":
        return 1.0 - np.exp(-K * x_values**2) * np.sinc(x_values / np.pi)
    if kernel == "pake":
        return 1.0 - np.exp(-K * x_values**2) * np.cos(x_values)
    if kernel == "weibull":
        return 1.0 - np.exp(-K * x_values**beta)
    if kernel == "a-l":
        return 1.0 - np.exp(-K * x_values**beta) * np.cos(x_values)

    raise ValueError(f"Unknown kernel: {kernel}. Use one of {VALID_KERNELS}")


def ndq_1d(time_values, mu, sigma, beta, kernel):
    p_values = p_gaussian(D_GRID, mu, sigma)
    return np.array([
        0.5 * np.trapz(p_values * dq_kernel(D_GRID * tau_i, kernel, beta), D_GRID)
        for tau_i in np.asarray(time_values, dtype=float)
    ])


def ndq_2d(time_values, mu1, sigma1, mu2, sigma2, frac1, beta, kernel):
    p_values = p_gaussian_2d(D_GRID, mu1, sigma1, mu2, sigma2, frac1)
    return np.array([
        0.5 * np.trapz(p_values * dq_kernel(D_GRID * tau_i, kernel, beta), D_GRID)
        for tau_i in np.asarray(time_values, dtype=float)
    ])


def ndq_single(time_values, dres):
    return 0.5 * (1.0 - np.exp(-K * dres**2 * np.asarray(time_values, dtype=float) ** 2))


def make_fit_model(kernel, n_components):
    kernel = kernel.lower()
    if kernel not in VALID_KERNELS:
        raise ValueError(f"kernel must be one of {VALID_KERNELS}")

    if n_components == 1:
        def model(time_values, mu, sigma, beta):
            return ndq_1d(time_values, mu, sigma, beta, kernel)
        return model

    if n_components == 2:
        def model(time_values, mu1, sigma1, mu2, sigma2, frac1, beta):
            return ndq_2d(time_values, mu1, sigma1, mu2, sigma2, frac1, beta, kernel)
        return model

    raise ValueError("n_components must be 1 or 2")


def fit_selected_model(time0, ndq0, kernel="gaussian", n_components=1):
    kernel = kernel.lower()
    if kernel not in VALID_KERNELS:
        raise ValueError(f"kernel must be one of {VALID_KERNELS}")

    model = make_fit_model(kernel, n_components)
    if n_components == 1:
        p0 = [0.25, 1e-3, 2.0]
        bounds_min = [1e-7, 1e-7, 1e-7]
        bounds_max = [1.000, 0.8, 6.0]
        param_names = ["mu", "sigma", "beta"]
    elif n_components == 2:
        p0 = [0.061, 0.05, 0.313, 0.033, 0.359, 0.96]
        bounds_min = [1e-3, 1e-4, 1e-5, 1e-15, 0.0, 1e-7]
        bounds_max = [1.0000, 1.5, 1.9000, 1.1, 1.0, 6.0]
        param_names = ["mu1", "sigma1", "mu2", "sigma2", "frac1", "beta"]
    else:
        raise ValueError("n_components must be 1 or 2")

    popt, pcov = curve_fit(
        model,
        time0,
        ndq0,
        p0=p0,
        bounds=(bounds_min, bounds_max),
        maxfev=5000000,
    )
    fit = model(time0, *popt)

    return {
        "kernel": kernel,
        "n_components": n_components,
        "popt": popt,
        "pcov": pcov,
        "fit": fit,
        "param_names": param_names,
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