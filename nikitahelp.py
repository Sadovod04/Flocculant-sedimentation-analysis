import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

def approx1(x, T1):
    return 600-100* np.exp(-(x) / T1)

time = np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100])

pressure = np.array([500, 501, 502, 504, 508, 519, 540, 570, 589, 595, 598, 600, 600])

(T1), _ = curve_fit(approx1, time, pressure, p0=[1])
print(T1)

# Вычисление коэффициента детерминации для апериодического звена 2 порядка
y_pred_1 = approx1(time, T1)
r_1, p_value_1 = pearsonr(pressure, y_pred_1)
r_squared_1 = r_1**2
print(f'Коэффициент детерминации для апериодического звена 1 порядка (R^2): {r_squared_1:.3f}')

# Графическое представление кривых
t = np.arange(0, 100, 0.1)
plt.plot(t, approx1(t, T1), color="#FFD700", label="Апериодическое звено 1 порядка")
plt.plot(time, pressure, color="black")
plt.scatter(time, pressure, color="black", label="Изначальные данные")
plt.legend()
plt.show()
