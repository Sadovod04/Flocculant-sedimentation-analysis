"""Fit a two-exponential relaxation model to a batch-settling curve.

Interface height h(t) is modelled as a sum of two decaying exponentials; the
fast term (a1, T1) is fixed from the early slope and ``curve_fit`` recovers the
slow term (a2, T2), which characterises the consolidation stage.
"""
import numpy as np
from scipy.optimize import curve_fit


def func(t, a2, T2):
    a1 = 110
    T1 = 104.5
    return 1 + a1 * np.exp(-t / T1) + a2 * np.exp(-t / T2)


# Измеренная кривая осаждения: время (мин) и высота раздела (см)
t = np.array([0,1,2,2.42,3.15,4.15,5,5.55,6,7,8,9,10,
                  10.5,11,12,12.5,13,14,14.5,15,15.5,16,16.5,17,18,19,20,20.3,
                  21,22,23,24,25,26,27,28,29 ,30,31,32,33,34,35,36,37,38,39,40,42,44,46,48,50])
h = np.array([117.6,116.2,115.2,114.7,113.7,112.5,110.3,108.3,107.1,106,104.1,102.1,100,98.7,97.9,
              97.3,95.9,95.1,94.4,93.6,92.7,92.1,91.8,91.3,90.7,90.3,89.3,88.3,87.6,86.9,86.2,85.6,84.9,84.5,
              84,83.4,83,82.7,82.2,81.9,81.4,81.1,80.8,80.5,80.2,80,79.7,79.6,79.4,79,78.7,78.5,78.3,78.3])

# Используем curve_fit для подгонки кривой к данным
popt, pcov = curve_fit(func, t, h, maxfev=10000)

# Выводим результат
a2, T2 = popt
print("Коэффициенты: a2 = {:.2f}, T2 = {:.2f}".format(a2, T2))
