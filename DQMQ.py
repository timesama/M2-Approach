import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

mpl.rc("font", size=14)

# ============================================================
# GLOBAL SETTINGS
# ============================================================
parent_directory = os.getcwd()

k = 0.4
D_grid = np.linspace(0, 0.94, 60000)

VALID_KERNELS = ["gaussian", "abragam", "pake", "weibull", "a-l"]


# ============================================================
# DATA IO
# ============================================================
def read_dqmq_file(material, filename):
    full_path = os.path.join(parent_directory, material, filename)

    return pd.read_csv(
        full_path,
        header=None,
        names=["tau", "DQ", "Ref", "nDQ"],
        sep=",",
        engine="python"
    )


def prepare_dqmq_arrays(df):
    ref_max = np.max(df["Ref"])

    tau = df["tau"].to_numpy(dtype=float)
    DQ = df["DQ"].to_numpy(dtype=float) / ref_max
    Ref = df["Ref"].to_numpy(dtype=float) / ref_max
    nDQ_original = df["nDQ"].to_numpy(dtype=float)

    nDQ = savgol_filter(nDQ_original, 3, 1)
    tau = tau + 1

    MQ_raw = df["DQ"].to_numpy(dtype=float) + df["Ref"].to_numpy(dtype=float)
    MQ = MQ_raw / np.max(MQ_raw)

    Time0 = np.insert(tau, 0, 0.0)
    DQ0 = np.insert(DQ, 0, 0.0)
    nDQ0 = np.insert(nDQ, 0, 0.0)

    return {
        "tau": tau,
        "DQ": DQ,
        "Ref": Ref,
        "nDQ": nDQ,
        "MQ": MQ,
        "Time0": Time0,
        "DQ0": DQ0,
        "nDQ0": nDQ0,
    }


# ============================================================
# DISTRIBUTION: ONLY GAUSSIAN P(Dres)
# ============================================================
def normalize(P, D):
    area = np.trapz(P, D)
    if area <= 0 or not np.isfinite(area):
        return np.zeros_like(P)
    return P / area


def P_gaussian(D, mu, sigma):
    sigma = max(float(sigma), 1e-9)
    P = np.exp(-(D - mu) ** 2 / (2 * sigma**2))
    return normalize(P, D)


def P_gaussian_2D(D, mu1, sigma1, mu2, sigma2, frac1):
    frac1 = np.clip(frac1, 0.0, 1.0)
    P1 = P_gaussian(D, mu1, sigma1)
    P2 = P_gaussian(D, mu2, sigma2)
    return frac1 * P1 + (1.0 - frac1) * P2


# ============================================================
# KERNELS
# ============================================================

def dq_kernel(x, kernel, beta=2.0):
    kernel = kernel.lower()
    beta = max(float(beta), 1e-12)

    if kernel == "gaussian":
        return 1.0 - np.exp(-k * x**2)

    elif kernel == "abragam":
        # np.sinc(x / np.pi) = sin(x) / x
        return 1.0 - np.exp(-k * x**2) * np.sinc(x / np.pi)

    elif kernel == "pake":
        return 1.0 - np.exp(-k * x**2) * np.cos(x)

    elif kernel == "weibull":
        return 1.0 - np.exp(-k * x**beta)

    elif kernel == "a-l":
        return 1.0 - np.exp(-k * x**beta) * np.cos(x)

    else:
        raise ValueError(f"Unknown kernel: {kernel}. Use one of {VALID_KERNELS}")


# ============================================================
# nDQ MODELS
# ============================================================
def ndq_1D(t, mu, sigma, beta, kernel):
    P = P_gaussian(D_grid, mu, sigma)

    y = []
    for tau_i in np.asarray(t, dtype=float):
        x = D_grid * tau_i
        K = dq_kernel(x, kernel, beta)
        y.append(0.5 * np.trapz(P * K, D_grid))

    return np.array(y)


def ndq_2D(t, mu1, sigma1, mu2, sigma2, frac1, beta, kernel):
    P = P_gaussian_2D(D_grid, mu1, sigma1, mu2, sigma2, frac1)

    y = []
    for tau_i in np.asarray(t, dtype=float):
        x = D_grid * tau_i
        K = dq_kernel(x, kernel, beta)
        y.append(0.5 * np.trapz(P * K, D_grid))

    return np.array(y)


def ndq_single(t, Dres):
    return 0.5 * (1.0 - np.exp(-k * Dres**2 * np.asarray(t)**2))


def make_fit_model(kernel, n_components):
    if n_components == 1:
        def model(t, mu, sigma, beta):
            return ndq_1D(t, mu, sigma, beta, kernel)
        return model

    elif n_components == 2:
        def model(t, mu1, sigma1, mu2, sigma2, frac1, beta):
            return ndq_2D(t, mu1, sigma1, mu2, sigma2, frac1, beta, kernel)
        return model

    else:
        raise ValueError("n_components must be 1 or 2")


# ============================================================
# FITTING
# ============================================================
def fit_selected_model(Time0, nDQ0, kernel="gaussian", n_components=1):
    kernel = kernel.lower()

    if kernel not in VALID_KERNELS:
        raise ValueError(f"kernel must be one of {VALID_KERNELS}")

    model = make_fit_model(kernel, n_components)

    if n_components == 1:
        p0 = [0.25, 1e-3, 2.0]

        bounds_min = [1e-7, 1e-7, 1e-7]
        bounds_max = [1.000, 0.8, 6.0]

        param_names = ["mu", "sigma", "beta"]

    else:
        p0 = [0.061, 0.05, 0.313, 0.033, 0.359, 0.96]

        bounds_min = [1e-3, 1e-4, 1e-5, 1e-15, 0.0, 0.1]
        bounds_max = [1.0000, 1.5, 1.9000, 1.1, 1.0, 4.0]

        param_names = ["mu1", "sigma1", "mu2", "sigma2", "frac1", "beta"]

    popt, pcov = curve_fit(
        model,
        # Time0[20:44],
        # nDQ0[20:44],
        Time0,
        nDQ0,
        p0=p0,
        bounds=(bounds_min, bounds_max),
        maxfev=5000000
    )

    fit = model(Time0, *popt)

    return {
        "kernel": kernel,
        "n_components": n_components,
        "popt": popt,
        "pcov": pcov,
        "fit": fit,
        "param_names": param_names,
    }


def fit_single_reference(Time0, nDQ0):
    popt, pcov = curve_fit(
        ndq_single,
        Time0,
        nDQ0,
        p0=[0.25],
        bounds=([0.001], [1.0]),
        maxfev=20000
    )

    return {
        "popt": popt,
        "pcov": pcov,
        "fit": ndq_single(Time0, *popt),
    }


# ============================================================
# PLOTTING
# ============================================================
def plot_single_fit(name, arrays, fit_result, single_result=None):
    tau = arrays["tau"]
    Time0 = arrays["Time0"]

    plt.figure(figsize=(7, 5))

    plt.hlines(0.5, 0, 300, "red", "dashed", linewidth=1)

    plt.plot(Time0, arrays["DQ0"], "k", linewidth=2, label=r"$I_{DQ}$")
    plt.plot(tau, arrays["Ref"], "r", linewidth=2, label=r"$I_{ref}$")
    plt.plot(Time0, arrays["nDQ0"], "b", linewidth=2, label=r"$I_{nDQ}$")
    plt.plot(tau, arrays["MQ"], "g", linewidth=2, label=r"$I_{MQ}$")

    label = f"{fit_result['kernel']} {fit_result['n_components']}D"
    plt.plot(Time0, fit_result["fit"], "m--", linewidth=3, label=label)

    # if single_result is not None:
    #     plt.plot(Time0, single_result["fit"], "c:", linewidth=2, label="single Dres")

    plt.xlabel(r"$\tau_{DQ}$, µs")
    plt.ylabel("Norm. amplitude, a.u.")
    plt.xlim([0, 100])
    plt.ylim([0, 1.1])
    plt.title(r"$\omega_{cr}$ = " + str(name))
    plt.legend()
    plt.tight_layout()


def plot_fitted_distribution(name, fit_result):
    D_plot = D_grid / (2 * np.pi) * 1000

    if fit_result["n_components"] == 1:
        mu, sigma, beta = fit_result["popt"]
        P = P_gaussian(D_grid, mu, sigma)

    else:
        mu1, sigma1, mu2, sigma2, frac1, beta = fit_result["popt"]
        P = P_gaussian_2D(D_grid, mu1, sigma1, mu2, sigma2, frac1)

    plt.figure(figsize=(7, 5))
    plt.plot(D_plot, P, linewidth=2)

    plt.xlabel(r"$D_{res}/2\pi$, kHz")
    plt.ylabel(r"$P(D_{res})$")
    plt.title(
        f"{name}: Gaussian P(Dres), "
        f"{fit_result['n_components']}D, kernel={fit_result['kernel']}"
    )
    plt.tight_layout()


def print_fit_result(name, fit_result, single_result=None):
    print("\n" + "=" * 70)
    print(f"Sample: {name}")
    print(f"Kernel: {fit_result['kernel']}")
    print(f"Dres components: {fit_result['n_components']}")

    for p_name, p_val in zip(fit_result["param_names"], fit_result["popt"]):
        print(f"{p_name:>8s} = {p_val:.8g}")

    if single_result is not None:
        print(f"{'Dres_single':>8s} = {single_result['popt'][0]:.8g}")

    print("=" * 70)


# ============================================================
# PROCESSING
# ============================================================
def process_one(material, filename, kernel="gaussian", n_components=1, do_plot=True):
    pattern = re.compile(r"DQMQ_data_(.*).csv")
    match = pattern.match(filename)

    if match is None:
        raise ValueError(f"Filename does not match DQMQ_data_*.csv: {filename}")

    name = match.group(1)

    df = read_dqmq_file(material, filename)
    arrays = prepare_dqmq_arrays(df)

    fit_result = fit_selected_model(
        arrays["Time0"],
        arrays["nDQ0"],
        kernel=kernel,
        n_components=n_components
    )

    single_result = fit_single_reference(arrays["Time0"], arrays["nDQ0"])

    print_fit_result(name, fit_result, single_result)

    if do_plot:
        plot_single_fit(name, arrays, fit_result, single_result)
        plot_fitted_distribution(name, fit_result)

    return name, {
        "arrays": arrays,
        "fit": fit_result,
        "single": single_result,
    }


def process_multiple(material, kernel="gaussian", n_components=1, do_plot=True):
    dir_path = os.path.join(parent_directory, material)

    results = {}

    for filename in sorted(os.listdir(dir_path)):
        if not filename.endswith(".csv"):
            continue

        if not filename.startswith("DQMQ_data_"):
            continue

        name, result = process_one(
            material,
            filename,
            kernel=kernel,
            n_components=n_components,
            do_plot=do_plot
        )

        results[name] = result

    return results


def process_selected_files(material, filenames, kernel="gaussian", n_components=1, do_plot=True):
    results = {}

    for filename in filenames:
        name, result = process_one(
            material,
            filename,
            kernel=kernel,
            n_components=n_components,
            do_plot=do_plot
        )

        results[name] = result

    return results


# ============================================================
# SUMMARY PLOTS
# ============================================================
def plot_parameter_vs_name(results, parameter="mu"):
    x = []
    y = []

    for name, result in results.items():
        fit = result["fit"]
        params = dict(zip(fit["param_names"], fit["popt"]))

        if parameter not in params:
            continue

        x.append(float(name))
        y.append(params[parameter])

    x = np.array(x)
    y = np.array(y)

    order = np.argsort(x)

    plt.figure(figsize=(6, 5))
    plt.plot(x[order], y[order], "o-", markersize=8)
    plt.xlabel(r"$\omega_{cr}$")
    plt.ylabel(parameter)
    plt.title(f"{parameter} vs crystallinity")
    plt.tight_layout()


def plot_all_distributions(results):
    plt.figure(figsize=(7, 5))

    D_plot = D_grid / (2 * np.pi) * 1000

    for name, result in results.items():
        fit = result["fit"]

        if fit["n_components"] == 1:
            mu, sigma, beta = fit["popt"]
            P = P_gaussian(D_grid, mu, sigma)

        else:
            mu1, sigma1, mu2, sigma2, frac1, beta = fit["popt"]
            P = P_gaussian_2D(D_grid, mu1, sigma1, mu2, sigma2, frac1)

        plt.plot(D_plot, P, linewidth=2, label=name)

    plt.xlabel(r"$D_{res}/2\pi$, kHz")
    plt.ylabel(r"$P(D_{res})$")
    plt.title("Fitted Gaussian Dres distributions")
    plt.legend()
    plt.tight_layout()


# ============================================================
# MAIN BODY
# ============================================================
if __name__ == "__main__":

    # Choose here
    MATERIAL = "PE"             # "PE", "PP", or "PS"
    KERNEL = "a-l"         # "gaussian", "abragam", "pake", "weibull", "a-l"
    N_COMPONENTS = 2            # 1 or 2

    # Either process all CSV files:
    results = process_multiple(
        MATERIAL,
        kernel=KERNEL,
        n_components=N_COMPONENTS,
        do_plot=True
    )

    # Or process only selected files:
    # results = process_selected_files(
    #     MATERIAL,
    #     filenames=[
    #         "DQMQ_data_0.353.csv",
    #         "DQMQ_data_0.398.csv",
    #     ],
    #     kernel=KERNEL,
    #     n_components=N_COMPONENTS,
    #     do_plot=True
    # )

    plot_all_distributions(results)


    plt.show()

    print("debug")