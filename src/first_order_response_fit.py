"""Fit a first-order (aperiodic) step response to pressure-vs-time data.

    p(t) = p_inf - A * exp(-t / T1)

``T1`` is recovered with ``scipy.optimize.curve_fit`` and the fit quality is
reported as R^2. The figure is saved to ``figures/``.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "figures")


def approx1(x, T1):
    return 600 - 100 * np.exp(-x / T1)


time = np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100])
pressure = np.array([500, 501, 502, 504, 508, 519, 540, 570, 589, 595, 598, 600, 600])

(T1,), _ = curve_fit(approx1, time, pressure, p0=[1])
print(f"T1 = {T1:.3f}")

y_pred = approx1(time, T1)
r, _ = pearsonr(pressure, y_pred)
print(f"R^2 (first-order model) = {r ** 2:.3f}")

t = np.arange(0, 100, 0.1)
plt.plot(t, approx1(t, T1), color="#FFD700", label="Первый порядок (аппроксимация)")
plt.plot(time, pressure, color="black")
plt.scatter(time, pressure, color="black", label="Исходные данные")
plt.xlabel("Время")
plt.ylabel("Давление")
plt.legend()
os.makedirs(FIG_DIR, exist_ok=True)
plt.savefig(os.path.join(FIG_DIR, "first_order_response_fit.png"), dpi=120, bbox_inches="tight")
print("saved figure -> figures/first_order_response_fit.png")
