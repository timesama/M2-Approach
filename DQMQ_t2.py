# ### 1. Read the DQMQ_T2_distribution.csv file
# ### 1 column DQ time, 2: DQAmp, 3:M2, 4: t2*, 5:t2*lin, 6:DQNorm

# ### Take T2*lin (call it T2 for simplicity)
# ### Take DQAmp (call it ampDQ)

# ### in each point of T2
# ### ampDQ * exp (-(t)/T2)^2)
# ### where t is going from 0 to the end of DQ time
# ### take a sum of that. 

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt

# # =========================
# # Helper: Min-Max normalization
# # =========================
# def normalize(y, z):
#     ymin = np.min(z)
#     ymax = np.max(z)
#     if ymax == ymin:
#         return np.zeros_like(z)  # avoid division by zero
#     return (y - ymin) / (ymax - ymin)

# # =========================
# # PART 1: T2 distribution → summed signal
# # =========================
# data_T2 = pd.read_csv("DQMQ_T2_distribution.csv", header=None)
# data_T2.columns = ["DQ_time", "DQAmp", "M2", "T2star", "T2star_lin", "DQNorm"]

# ampDQ = data_T2["DQAmp"].values
# T2 = data_T2["T2star_lin"].values
# t_DQ = data_T2["DQ_time"].values

# t = np.linspace(0, t_DQ.max(), 500)

# dT2 = np.diff(T2)
# dT2 = np.append(dT2, dT2[-1])  # make same length as T2

# signal = np.zeros_like(t)
# for amp, t2, deltat2 in zip(ampDQ, T2, dT2):
#     signal += amp * np.exp(-(t / t2)**2) * deltat2

# # Normalize summed signal
# signal_norm = normalize(signal, signal)

# # =========================
# # PART 2: Dres distribution
# # =========================
# data_Dres = pd.read_csv("DQMQ_Dres_distribution.csv", header=None)

# tauDQ = data_Dres.iloc[:, 0].values
# IDQ = data_Dres.iloc[:, 1].values
# IRef = data_Dres.iloc[:, 2].values
# nDQ = data_Dres.iloc[:, 3].values

# # Normalize each curve independently
# DQ_normalized = normalize(IDQ, IRef)
# REF_normalized = normalize(IRef, IRef)
# IMQ = IDQ + IRef
# MQ_normalized = normalize(IMQ, IMQ)


# # =========================
# # PART 3: Plot
# # =========================
# plt.figure()

# plt.plot(tauDQ, DQ_normalized, label="DQ")
# plt.plot(tauDQ, REF_normalized, label="Ref")
# plt.plot(tauDQ, MQ_normalized, label="MQ", linewidth=3)
# plt.plot(tauDQ, nDQ, label="nDQ")

# plt.plot(t+9, signal_norm, label="Summed signal (norm)", linewidth=3)

# plt.xlabel("tauDQ / t")
# plt.ylabel("Normalized Signal")
# plt.legend()

# plt.show()

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# SETTINGS
# ============================================================
parent_directory = Path.cwd()

dqmq_pattern = re.compile(r"DQMQ_data_(.*)\.csv")
dqt2_pattern = re.compile(r"DQT2_data_(.*)\.csv")

output_directory = parent_directory / "figures_T2_DQMQ"
output_directory.mkdir(exist_ok=True)

SHOW_FIGURES = False
SAVE_FIGURES = True


# ============================================================
# Helper: Min-Max normalization
# ============================================================
def normalize(y, reference=None):
    """
    Min-max normalize y.

    If reference is given, use min/max of reference.
    If reference is None, use min/max of y.
    """
    y = np.asarray(y, dtype=float)

    if reference is None:
        reference = y
    else:
        reference = np.asarray(reference, dtype=float)

    ymin = np.min(reference)
    ymax = np.max(reference)

    if ymax == ymin:
        return np.zeros_like(y)

    return (y - ymin) / (ymax - ymin)


# ============================================================
# Read files
# ============================================================
def read_dqt2_file(filename):
    """
    DQT2 file columns:
    1: DQ_time
    2: DQAmp
    3: M2
    4: T2star
    5: T2star_lin
    6: DQNorm
    """
    df = pd.read_csv(filename, header=None)
    df.columns = ["DQ_time", "DQAmp", "M2", "T2star", "T2star_lin", "DQNorm"]
    return df


def read_dqmq_file(filename):
    """
    DQMQ file columns:
    1: tauDQ
    2: IDQ
    3: IRef
    4: nDQ
    """
    df = pd.read_csv(filename, header=None)
    df.columns = ["tauDQ", "IDQ", "IRef", "nDQ"]
    return df


# ============================================================
# T2 distribution -> summed signal
# ============================================================
def calculate_t2_summed_signal(df_t2, n_points=500):
    ampDQ = df_t2["DQAmp"].to_numpy(dtype=float)
    T2 = df_t2["T2star_lin"].to_numpy(dtype=float)
    t_DQ = df_t2["DQ_time"].to_numpy(dtype=float)

    valid = np.isfinite(ampDQ) & np.isfinite(T2) & np.isfinite(t_DQ) & (T2 > 0)
    ampDQ = ampDQ[valid]
    T2 = T2[valid]
    t_DQ = t_DQ[valid]

    order = np.argsort(T2)
    ampDQ = ampDQ[order]
    T2 = T2[order]

    t = np.linspace(0, np.max(t_DQ), n_points)

    if len(T2) > 1:
        dT2 = np.gradient(T2)
    else:
        dT2 = np.ones_like(T2)

    signal = np.zeros_like(t)

    for amp, t2, delta_t2 in zip(ampDQ, T2, dT2):
        signal += amp * np.exp(-(t / t2) ** 2) * delta_t2

    signal_norm = normalize(signal)

    return t, signal, signal_norm


# ============================================================
# DQMQ curves
# ============================================================
def calculate_dqmq_curves(df_dqmq):
    tauDQ = df_dqmq["tauDQ"].to_numpy(dtype=float)
    IDQ = df_dqmq["IDQ"].to_numpy(dtype=float)
    IRef = df_dqmq["IRef"].to_numpy(dtype=float)
    nDQ = df_dqmq["nDQ"].to_numpy(dtype=float)

    DQ_normalized = normalize(IDQ, IRef)
    REF_normalized = normalize(IRef, IRef)

    IMQ = IDQ + IRef
    MQ_normalized = normalize(IMQ)

    return {
        "tauDQ": tauDQ,
        "DQ_normalized": DQ_normalized,
        "REF_normalized": REF_normalized,
        "MQ_normalized": MQ_normalized,
        "nDQ": nDQ,
    }


# ============================================================
# Plot one paired dataset
# ============================================================
def plot_one_pair(name, dqmq_file, dqt2_file, shift=9.0):
    df_dqmq = read_dqmq_file(dqmq_file)
    df_t2 = read_dqt2_file(dqt2_file)

    curves = calculate_dqmq_curves(df_dqmq)
    t, signal, signal_norm = calculate_t2_summed_signal(df_t2)

    plt.figure(figsize=(7, 5))

    plt.plot(curves["tauDQ"], curves["DQ_normalized"], label="DQ")
    plt.plot(curves["tauDQ"], curves["REF_normalized"], label="Ref")
    plt.plot(curves["tauDQ"], curves["MQ_normalized"], label="MQ", linewidth=3)
    plt.plot(curves["tauDQ"], curves["nDQ"], label="nDQ")

    plt.plot(
        t + shift,
        signal_norm,
        label=f"T2 summed signal, shifted by {shift:g}",
        linewidth=3
    )

    plt.xlabel(r"$\tau_{DQ}$ / t")
    plt.ylabel("Normalized signal")
    plt.title(f"DQMQ and T2 summed signal: {name}")
    plt.legend()
    plt.tight_layout()

    if SAVE_FIGURES:
        output_file = output_directory / f"DQMQ_T2_overlay_{name}.png"
        plt.savefig(output_file, dpi=300)
        print(f"Saved: {output_file}")

    if SHOW_FIGURES:
        plt.show()
    else:
        plt.close()


# ============================================================
# File pairing
# ============================================================
def find_paired_files(directory):
    dqmq_files = {}
    dqt2_files = {}

    for file in directory.glob("*.csv"):
        m_dqmq = dqmq_pattern.match(file.name)
        m_dqt2 = dqt2_pattern.match(file.name)

        if m_dqmq:
            name = m_dqmq.group(1)
            dqmq_files[name] = file

        if m_dqt2:
            name = m_dqt2.group(1)
            dqt2_files[name] = file

    common_names = sorted(set(dqmq_files.keys()) & set(dqt2_files.keys()), key=float)

    missing_dqt2 = sorted(set(dqmq_files.keys()) - set(dqt2_files.keys()))
    missing_dqmq = sorted(set(dqt2_files.keys()) - set(dqmq_files.keys()))

    if missing_dqt2:
        print("DQMQ files without matching DQT2 files:")
        for name in missing_dqt2:
            print(f"  DQMQ_data_{name}.csv")

    if missing_dqmq:
        print("DQT2 files without matching DQMQ files:")
        for name in missing_dqmq:
            print(f"  DQT2_data_{name}.csv")

    pairs = [
        (name, dqmq_files[name], dqt2_files[name])
        for name in common_names
    ]

    return pairs


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    pairs = find_paired_files(parent_directory)

    print(f"Found {len(pairs)} matched file pairs.")

    for name, dqmq_file, dqt2_file in pairs:
        print(f"Processing {name}:")
        print(f"  DQMQ: {dqmq_file.name}")
        print(f"  DQT2: {dqt2_file.name}")

        plot_one_pair(
            name=name,
            dqmq_file=dqmq_file,
            dqt2_file=dqt2_file,
            shift=9.0
        )

    print("Done.")