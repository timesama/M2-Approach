"""Pure signal-analysis helpers for the DQMQ tab."""

import logging

import numpy as np
from scipy.optimize import curve_fit

logger = logging.getLogger(__name__)


def normalize_to_reference_max(dq_raw, ref_raw):
    dq_raw = np.asarray(dq_raw, dtype=float)
    ref_raw = np.asarray(ref_raw, dtype=float)
    finite_ref = ref_raw[np.isfinite(ref_raw)]
    if finite_ref.size == 0:
        raise ValueError("Ref data does not contain finite values for normalization.")

    ref_max = np.max(finite_ref)
    if not np.isfinite(ref_max) or np.isclose(ref_max, 0.0):
        raise ValueError(
            "Ref maximum is zero or non-finite; cannot normalize DQMQ data."
        )

    dq_norm = dq_raw / ref_max
    ref_norm = ref_raw / ref_max
    finite_dq_norm = dq_norm[np.isfinite(dq_norm)]
    if finite_dq_norm.size > 0:
        dq_norm_max = np.max(finite_dq_norm)
        if dq_norm_max > 1.05:
            logger.warning(
                "DQMQ normalized DQ exceeds 1.05: max DQ_norm=%s raw DQ max=%s raw Ref max=%s",
                dq_norm_max,
                np.nanmax(dq_raw),
                ref_max,
            )

    return dq_norm, ref_norm, ref_max


def normalize_to_own_max(signal):
    signal = np.asarray(signal, dtype=float)
    finite_signal = signal[np.isfinite(signal)]
    if finite_signal.size == 0:
        raise ValueError("Signal does not contain finite values for normalization.")

    signal_max = np.max(finite_signal)
    if not np.isfinite(signal_max) or np.isclose(signal_max, 0.0):
        raise ValueError(
            "Signal maximum is zero or non-finite; cannot normalize signal."
        )

    return signal / signal_max, signal_max


def fit_tail(time, dq_norm, ref_norm, fit_from, fit_to, exponent):
    time = np.asarray(time, dtype=float)
    dq_norm = np.asarray(dq_norm, dtype=float)
    ref_norm = np.asarray(ref_norm, dtype=float)
    tail_target = ref_norm - dq_norm
    idx_min = _nearest_index(time, fit_from)
    idx_max = _nearest_index(time, fit_to)
    start_idx = min(idx_min, idx_max)
    stop_idx = max(idx_min, idx_max) + 1
    time_cut = time[start_idx:stop_idx]
    tail_cut = tail_target[start_idx:stop_idx]
    valid_fit_points = np.isfinite(time_cut) & np.isfinite(tail_cut)
    time_cut = time_cut[valid_fit_points]
    tail_cut = tail_cut[valid_fit_points]

    if len(time_cut) < 3:
        logger.warning("DQMQ tail fitting skipped: fewer than 3 valid fit points.")
        return np.zeros_like(time, dtype=float)

    def exponent_model(x_values, amplitude, tau, offset):
        safe_tau = tau + 1e-12
        scaled_time = x_values / safe_tau
        exponent_values = -(scaled_time**exponent)
        return amplitude * np.exp(exponent_values) + offset

    try:
        amplitude0 = np.max(tail_cut) - np.min(tail_cut)
        tau0 = max(1e-6, (time_cut[-1] - time_cut[0]) / 4.0)
        offset0 = np.median(tail_cut[-5:])
        initial_params = [amplitude0, tau0, offset0]
        fitted_params, _ = curve_fit(
            exponent_model,
            time_cut,
            tail_cut,
            p0=initial_params,
            maxfev=10000000,
        )
        return exponent_model(time, *fitted_params)
    except Exception:
        logger.exception("DQMQ tail fitting failed; using zero tail fit.")
        return np.zeros_like(time, dtype=float)


def calculate_dqmq_analysis(
    raw_time,
    dq_raw,
    ref_raw,
    fit_from,
    fit_to,
    fitting_exponent,
    noise_level,
    time_shift,
    smoothing=None,
):
    raw_time = np.asarray(raw_time, dtype=float)
    time = raw_time + time_shift
    dq_norm, ref_norm, ref_max = normalize_to_reference_max(dq_raw, ref_raw)
    mq_raw_norm = dq_norm + ref_norm
    mq_norm, mq_max = normalize_to_own_max(mq_raw_norm)
    tail_fit = fit_tail(time, dq_norm, ref_norm, fit_from, fit_to, fitting_exponent)
    additive = exp_apodization(time, fit_from, noise_level)
    denominator = mq_raw_norm - tail_fit + 2 * noise_level * additive
    numerator = dq_norm + noise_level * additive
    n_dq = calculate_safe_ndq(numerator, denominator)
    n_dq = smooth_ndq_if_requested(time, n_dq, smoothing)
    time0 = np.insert(time, 0, 0.0)
    n_dq0 = np.insert(n_dq, 0, 0.0)

    return {
        "time": time,
        "dq_norm": dq_norm,
        "ref_norm": ref_norm,
        "mq_raw_norm": mq_raw_norm,
        "mq_norm": mq_norm,
        "tail_fit": tail_fit,
        "additive": additive,
        "denominator": denominator,
        "time0": time0,
        "nDQ": n_dq0,
        "fit_from": fit_from,
        "fit_to": fit_to,
        "power": fitting_exponent,
        "noise": noise_level,
        "time_shift": time_shift,
        "smoothing": smoothing,
        "ref_max": ref_max,
        "mq_max": mq_max,
    }


def exp_apodization(time, fit_from, noise_level):
    time = np.asarray(time, dtype=float)
    noise_tau = fit_from * 0.19 if fit_from != 0 else 1e-9
    exp_function = np.empty_like(time)
    before_fit = time < fit_from
    before_scaled = (time[before_fit] - fit_from) / noise_tau
    after_scaled = (fit_from - time[~before_fit]) / noise_tau
    exp_function[before_fit] = np.exp(before_scaled**3)
    exp_function[~before_fit] = 1.0 - np.exp(after_scaled**3)
    return 0.5 * noise_level * exp_function


def calculate_safe_ndq(numerator, denominator):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    n_dq = np.full_like(numerator, np.nan, dtype=float)
    valid_denominator = np.isfinite(denominator) & (denominator > 0)
    valid_numerator = np.isfinite(numerator)
    valid_points = valid_denominator & valid_numerator
    n_dq[valid_points] = numerator[valid_points] / denominator[valid_points]
    invalid_count = len(n_dq) - np.count_nonzero(valid_points)
    if invalid_count:
        logger.warning(
            "DQMQ nDQ calculation skipped %d point(s) with invalid denominator or numerator.",
            invalid_count,
        )
    return n_dq


def smooth_ndq_if_requested(time, n_dq, smoothing):
    if smoothing is None:
        return n_dq

    smooth_from, smooth_to, smooth_window = smoothing
    smoothing_enabled = (
        smooth_window > 1
        and smooth_to > smooth_from
        and np.isfinite(smooth_from)
        and np.isfinite(smooth_to)
    )
    if not smoothing_enabled:
        return n_dq

    smoothed_n_dq = np.array(n_dq, copy=True)
    start_idx = _nearest_index(time, smooth_from)
    stop_idx = _nearest_index(time, smooth_to)
    start_idx, stop_idx = sorted((start_idx, stop_idx))
    stop_idx += 1
    segment = smoothed_n_dq[start_idx:stop_idx]
    finite_segment = np.isfinite(segment)
    if np.count_nonzero(finite_segment) < smooth_window:
        return smoothed_n_dq

    filled_segment = np.array(segment, copy=True)
    if not np.all(finite_segment):
        valid_indices = np.flatnonzero(finite_segment)
        all_indices = np.arange(len(segment))
        filled_segment = np.interp(all_indices, valid_indices, segment[finite_segment])

    averaged_segment = moving_average(filled_segment, smooth_window)
    averaged_segment[~finite_segment] = np.nan
    smoothed_n_dq[start_idx:stop_idx] = averaged_segment
    return smoothed_n_dq


def moving_average(values, window_size):
    kernel = np.ones(window_size) / window_size
    pad_width = window_size // 2
    padded = np.pad(values, pad_width, mode="edge")
    smoothed = np.convolve(padded, kernel, mode="valid")
    return smoothed[: len(values)]


def _nearest_index(values, target):
    values = np.asarray(values, dtype=float)
    return int(np.nanargmin(np.abs(values - target)))
