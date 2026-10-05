#fid fit calculator module

import numpy as np
from scipy.optimize import curve_fit

def weibull(t, A, T2, n):
    return A * np.exp(-((t / T2) ** n))


def lorentz(t, A, T2):
    return A * np.exp(-(t / T2))


def gaussian(t, A, T2):
    return A * np.exp(-((t / T2) ** 2))


def abragam(t, A, T2, B):
    gaussian_part = np.exp(- (t /(2* T2)) ** 2)

    # np.sinc(x) = sin(pi*x)/(pi*x)
    # therefore np.sinc(B*t/pi) = sin(B*t)/(B*t)
    oscillating_part = np.sinc((B * t) / np.pi)

    return A * gaussian_part * oscillating_part


def pake(t, A, T2, B):
    return (A * np.exp(-(t / (2*T2)) ** 2) * np.cos(B * t / 2))


MODEL_SPECS = {
    "Weibull": {
        "function": weibull,
        "equation": "A · exp[-(t / T₂)ⁿ]",
        "parameters": [
            {
                "name": "A",
                "label": "A",
                "default": 1.0,
                "minimum": 0,
                "maximum": 1e6,
            },
            {
                "name": "T2",
                "label": "T₂",
                "default": 50.0,
                "minimum": 1e-9,
                "maximum": 1e9,
            },
            {
                "name": "n",
                "label": "n",
                "default": 2.0,
                "minimum": 0.01,
                "maximum": 10.0,
            },
        ],
    },

    "Lorentz": {
        "function": lorentz,
        "equation": "A · exp(-t / T₂)",
        "parameters": [
            {
                "name": "A",
                "label": "A",
                "default": 1.0,
                "minimum": 0,
                "maximum": 1e6,
            },
            {
                "name": "T2",
                "label": "T₂",
                "default": 100.0,
                "minimum": 1e-9,
                "maximum": 1e9,
            },
        ],
    },

    "Gaussian": {
        "function": gaussian,
        "equation": "A · exp[-(t / T₂)²]",
        "parameters": [
            {
                "name": "A",
                "label": "A",
                "default": 1.0,
                "minimum": 0,
                "maximum": 1e6,
            },
            {
                "name": "T2",
                "label": "T₂",
                "default": 30.0,
                "minimum": 1e-9,
                "maximum": 1e9,
            },
        ],
    },

    "Abragam": {
        "function": abragam,
        "equation": "A · exp[-(t / 2 T₂)²] · sin(bt)/(bt)",
        "parameters": [
            {
                "name": "A",
                "label": "A",
                "default": 1.0,
                "minimum": 0,
                "maximum": 1e6,
            },
            {
                "name": "T2",
                "label": "T₂",
                "default": 10.0,
                "minimum": 1e-9,
                "maximum": 1e9,
            },
            {
                "name": "B",
                "label": "B",
                "default": 0.2,
                "minimum": 0.0,
                "maximum": 1e6,
            },
        ],
    },

    "Pake": {
        "function": pake,
        "equation": "A · exp[-(t / 2 T₂)²] · cos(bt/2)",
        "parameters": [
            {
                "name": "A",
                "label": "A",
                "default": 1.0,
                "minimum": 0,
                "maximum": 1e6,
            },
            {
                "name": "T2",
                "label": "T₂",
                "default": 10.0,
                "minimum": 1e-9,
                "maximum": 1e9,
            },
            {
                "name": "B",
                "label": "B",
                "default": 0.2,
                "minimum": 0.0,
                "maximum": 1e6,
            },
        ],
    },
}


def fit_components(time, signal, components):
    free_parameters, p0, lower, upper = _prepare_fit(components)
    total_function, free_lookup = _build_total_function(components, free_parameters)

    if free_parameters:
        popt, pcov = curve_fit(
            total_function,
            time,
            signal,
            p0=p0,
            bounds=(lower, upper),
            maxfev=50000,
        )
    else:
        popt = np.array([])
        pcov = np.empty((0, 0))

    return _build_fit_result(
        time,
        signal,
        components,
        free_parameters,
        free_lookup,
        popt,
        pcov,
    )


def _build_fit_result(time, signal, components, free_parameters, free_lookup, popt, pcov):
    parameter_results = []
    contributions = {}

    # Reconstruct final parameters and individual component curves
    for component in components:
        final_parameters = {}

        for parameter in component["parameters"]:
            key = (component["row"], parameter["name"])

            if parameter["fixed"]:
                value = parameter["value"]
                error = np.nan
                free_index = None
            else:
                free_index = free_lookup[key]
                value = popt[free_index]

                if pcov.size:
                    variance = pcov[free_index, free_index]
                    error = np.sqrt(variance) if variance >= 0 else np.nan
                else:
                    error = np.nan

            final_parameters[parameter["name"]] = value

            parameter_results.append({
                "row": component["row"],
                "model": component["model"],
                "name": parameter["name"],
                "label": parameter["label"],
                "value": float(value),
                "error": float(error) if np.isfinite(error) else np.nan,
                "fixed": parameter["fixed"],
                "free_index": free_index,
            })

        amplitudes = {row: 0.0 for row in range(1, 5)}

        for parameter in parameter_results:
            if parameter["name"] == "A":
                amplitudes[parameter["row"]] = parameter["value"]

        total_amplitude = sum(amplitudes.values())

        if total_amplitude != 0:
            amplitude_ratios = {row: amplitudes[row] / total_amplitude for row in range(1, 5)}
        else:
            amplitude_ratios = {row: 0.0 for row in range(1, 5)}

        contributions[component["row"]] = component["function"](time, **final_parameters)

    # Cumulative fitted curve
    cumulative = np.zeros_like(time, dtype=float)

    for contribution in contributions.values():
        cumulative += contribution

    # Residuals and statistics
    residuals = signal - cumulative
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((signal - np.mean(signal)) ** 2)

    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    n_free = len(free_parameters)
    degrees_of_freedom = len(signal) - n_free
    reduced_residual_variance = ss_res / degrees_of_freedom if degrees_of_freedom > 0 else np.nan

    # Covariance -> correlation
    correlation = covariance_to_correlation(pcov)

    # Strongest parameter correlation for display in the table
    dependencies = _parameter_dependencies(parameter_results, correlation)

    return {
        "parameters": parameter_results,
        "contributions": contributions,
        "cumulative": cumulative,
        "residuals": residuals,
        "covariance": pcov,
        "correlation": correlation,
        "dependency": dependencies,
        "r_squared": r_squared,
        "reduced_residual_variance": reduced_residual_variance,
        "amplitude_ratios": amplitude_ratios,
    }

def _prepare_fit(components):
    free_parameters = []
    p0, lower, upper = [], [], []

    for component in components:
        for parameter in component["parameters"]:
            if parameter["fixed"]:
                continue

            free_parameters.append((component, parameter))
            p0.append(parameter["value"])
            lower.append(parameter["minimum"])
            upper.append(parameter["maximum"])

    return (
        free_parameters,
        np.asarray(p0, dtype=float),
        np.asarray(lower, dtype=float),
        np.asarray(upper, dtype=float),
    )

def _build_total_function(components, free_parameters):
    free_lookup = {
        (component["row"], parameter["name"]): index
        for index, (component, parameter) in enumerate(free_parameters)
    }

    def total_function(t, *free_values):
        result = np.zeros_like(t, dtype=float)

        for component in components:
            values = {}

            for parameter in component["parameters"]:
                key = (component["row"], parameter["name"])
                values[parameter["name"]] = parameter["value"] if parameter["fixed"] else free_values[free_lookup[key]]

            result += component["function"](t, **values)

        return result

    return total_function, free_lookup

def covariance_to_correlation(pcov):
    if pcov.size == 0:
        return np.empty((0, 0))

    std = np.sqrt(np.clip(np.diag(pcov), 0, None))
    denominator = np.outer(std, std)

    return np.divide(
        pcov,
        denominator,
        out=np.full_like(pcov, np.nan, dtype=float),
        where=denominator != 0,
    )

def _parameter_dependencies(parameter_results, correlation):
    dependencies = {}

    free_results = [result for result in parameter_results if not result["fixed"]]

    for result in parameter_results:
        key = (result["row"], result["name"])

        if result["fixed"]:
            dependencies[key] = "Fixed"
            continue

        i = result["free_index"]

        if correlation.size == 0 or len(free_results) <= 1:
            dependencies[key] = "—"
            continue

        candidates = []

        for other in free_results:
            j = other["free_index"]

            if i == j:
                continue

            rho = correlation[i, j]

            if np.isfinite(rho):
                candidates.append((abs(rho), rho, other))

        if not candidates:
            dependencies[key] = "—"
            continue

        _, rho, other = max(candidates, key=lambda item: item[0])

        dependencies[key] = f"{rho:.3f} ({other['row']}:{other['label']})"

    return dependencies