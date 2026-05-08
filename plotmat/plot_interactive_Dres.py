import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from scipy.integrate import trapezoid as trapz

# -----------------------------
# Constants / grid
# -----------------------------
k = 0.4
D_grid = np.linspace(0, 0.60, 6000)
tau = np.linspace(0, 100, 500)


# -----------------------------
# Probability distributions
# -----------------------------
def normalize(P, D):
    area = trapz(P, D)
    if area <= 0:
        return np.zeros_like(P)
    return P / area


def P_gaussian(D, mu, sigma):
    sigma = max(sigma, 1e-6)
    P = np.exp(-(D - mu) ** 2 / (2 * sigma ** 2))
    return normalize(P, D)


def P_abragam(D, D0, w):
    w = max(w, 1e-6)
    P = (D / w**2) * np.exp(-(D - D0) ** 2 / (2 * w**2))
    P[P < 0] = 0
    return normalize(P, D)


def P_weibull(D, lam, beta):
    lam = max(lam, 1e-12)
    beta = max(beta, 1e-12)

    P = (beta / lam) * (D / lam) ** (beta - 1) * np.exp(-(D / lam) ** beta)
    P[D <= 0] = 0

    return normalize(P, D)


def P_pake(D, D0, w):
    # Pake-like singular distribution, regularized by small epsilon
    w = max(w, 1e-6)
    eps = 1e-12

    x = (D - D0) / w
    P = np.zeros_like(D)

    mask = np.abs(x) < 1
    P[mask] = 1.0 / np.sqrt(1.0 - x[mask] ** 2 + eps)

    return normalize(P, D)


def P_gaussian_mix2(D, mu1, sigma1, mu2, sigma2, frac1):
    return frac1 * P_gaussian(D, mu1, sigma1) + (1 - frac1) * P_gaussian(D, mu2, sigma2)


def P_abragam_mix2(D, D01, w1, D02, w2, frac1):
    return frac1 * P_abragam(D, D01, w1) + (1 - frac1) * P_abragam(D, D02, w2)


def P_pake_mix2(D, D01, w1, D02, w2, frac1):
    return frac1 * P_pake(D, D01, w1) + (1 - frac1) * P_pake(D, D02, w2)


def P_weibull_mix2(D, lam1, beta1, lam2, beta2, frac1):
    return frac1 * P_weibull(D, lam1, beta1) + (1 - frac1) * P_weibull(D, lam2, beta2)


# -----------------------------
# nDQ models
# -----------------------------
def ndq_from_distribution(t, P, kernel_model, beta=2.0):
    result = []
    beta = max(beta, 1e-12)

    for tau_i in t:
        x = D_grid * tau_i

        if kernel_model == "Gaussian":
            kernel = 1 - np.exp(-k * x**2)

        elif kernel_model == "Abragam":
            kernel = 1 - np.exp(-k * x**2) * np.sinc(x)

        elif kernel_model == "Pake":
            kernel = 1 - np.exp(-k * x**2) * np.cos(x)

        elif kernel_model == "Weibull":
            kernel = 1 - np.exp(-k * x**beta)

        elif kernel_model == "A-L":
            kernel = 1 - (np.exp(-k * x**beta) * np.cos(x))

        else:
            raise ValueError("Unknown kernel model")

        result.append(0.5 * trapz(P * kernel, D_grid))

    return np.array(result)


def open_interactive_dres(parent=None):
    """Open the interactive Dres matplotlib tool and return its figure."""
    # -----------------------------
    # Figure layout
    # -----------------------------
    fig, (ax_curve, ax_dist) = plt.subplots(1, 2, figsize=(12, 7))
    if fig.canvas.manager is not None:
        fig.canvas.manager.set_window_title("Interactive Dres")
    fig.subplots_adjust(left=0.28, bottom=0.4)

    ndq0 = np.zeros_like(tau)
    P0 = np.zeros_like(D_grid)

    (line_curve,) = ax_curve.plot(tau, ndq0, lw=2)
    (line_dist,) = ax_dist.plot(D_grid / (2 * np.pi) * 1000, P0, lw=2)

    ax_curve.set_xlabel(r"$\tau_{DQ}$, $\mu$s")
    ax_curve.set_ylabel(r"$I_{nDQ}$")
    ax_curve.set_xlim(0, 100)
    ax_curve.set_ylim(0, 0.55)
    ax_curve.grid(True)

    ax_dist.set_xlabel(r"$D_{res}/2\pi$, kHz")
    ax_dist.set_ylabel(r"$P(D_{res})$")
    ax_dist.set_xlim(0, 100)
    ax_dist.grid(True)

    # -----------------------------
    # Sliders
    # -----------------------------
    ax_s1 = fig.add_axes([0.28, 0.28, 0.55, 0.03])
    ax_s2 = fig.add_axes([0.28, 0.24, 0.55, 0.03])
    ax_s3 = fig.add_axes([0.28, 0.20, 0.55, 0.03])
    ax_s4 = fig.add_axes([0.28, 0.16, 0.55, 0.03])
    ax_s5 = fig.add_axes([0.28, 0.12, 0.55, 0.03])
    ax_s6 = fig.add_axes([0.28, 0.08, 0.55, 0.03])

    s1 = Slider(ax_s1, "center 1", 0.001, 0.5, valinit=0.15, valstep=0.001)
    s2 = Slider(ax_s2, "width 1", 0.000000001, 0.1, valinit=0.001, valstep=0.000000001)
    s3 = Slider(ax_s3, "center 2", 0.1, 0.5, valinit=0.35, valstep=0.001)
    s4 = Slider(ax_s4, "width 2", 0.000000001, 0.1, valinit=0.001, valstep=0.000000001)
    s5 = Slider(ax_s5, "fraction 1", 0.0, 1.0, valinit=0.5, valstep=0.001)
    s6 = Slider(ax_s6, "Weibull beta", 0.1, 5.0, valinit=2.0, valstep=0.01)

    # -----------------------------
    # Radio buttons
    # -----------------------------
    ax_radio_model = fig.add_axes([0.04, 0.62, 0.16, 0.16])
    radio_model = RadioButtons(
        ax_radio_model,
        ["Gaussian", "Abragam", "Pake", "Weibull", "A-L"],
    )

    ax_radio_order = fig.add_axes([0.04, 0.42, 0.16, 0.14])
    radio_order = RadioButtons(ax_radio_order, ["1D", "2D"])

    def compute_curve():
        kernel_model = radio_model.value_selected
        order = radio_order.value_selected

        mu1 = s1.val
        sigma1 = s2.val
        mu2 = s3.val
        sigma2 = s4.val
        frac1 = s5.val
        beta = s6.val

        if order == "1D":
            P = P_gaussian(D_grid, mu1, sigma1)

        elif order == "2D":
            P = P_gaussian_mix2(D_grid, mu1, sigma1, mu2, sigma2, frac1)

        else:
            raise ValueError("Unknown distribution order")

        ndq = ndq_from_distribution(tau, P, kernel_model, beta=beta)
        return ndq, P

    # -----------------------------
    # Update function
    # -----------------------------
    def update(_=None):
        ndq, P = compute_curve()

        line_curve.set_ydata(ndq)
        line_dist.set_ydata(P)

        max_p = max(P)
        ax_dist.set_ylim(0, max_p * 1.1 if max_p > 0 else 1)

        model = radio_model.value_selected
        order = radio_order.value_selected

        ax_curve.set_title(f"{model} {order}")
        ax_dist.set_title(r"$P(D_{res})$")

        fig.canvas.draw_idle()

    slider_callback_ids = []
    for slider in [s1, s2, s3, s4, s5, s6]:
        slider_callback_ids.append((slider, slider.on_changed(update)))

    radio_callback_ids = [
        (radio_model, radio_model.on_clicked(update)),
        (radio_order, radio_order.on_clicked(update)),
    ]

    fig._interactive_dres_widgets = {
        "sliders": [s1, s2, s3, s4, s5, s6],
        "radio_model": radio_model,
        "radio_order": radio_order,
        "slider_callback_ids": slider_callback_ids,
        "radio_callback_ids": radio_callback_ids,
    }
    fig._interactive_dres_closing = False

    def on_close(event):
        closing_figure = event.canvas.figure
        if getattr(closing_figure, "_interactive_dres_closing", False):
            return

        closing_figure._interactive_dres_closing = True
        widgets = getattr(closing_figure, "_interactive_dres_widgets", {}) or {}
        for widget, callback_id in widgets.get("slider_callback_ids", []):
            widget.disconnect(callback_id)
            widget.disconnect_events()
        for radio, callback_id in widgets.get("radio_callback_ids", []):
            radio.disconnect(callback_id)
            radio.disconnect_events()
        closing_figure._interactive_dres_widgets = None
        if plt.fignum_exists(closing_figure.number):
            plt.close(closing_figure)

    close_cid = fig.canvas.mpl_connect("close_event", on_close)
    fig._interactive_dres_close_cid = close_cid

    update()
    fig.show()
    plt.show(block=False)
    return fig


if __name__ == "__main__":
    open_interactive_dres()
