"""Vertical solids-concentration profile in the thickener (Simpson variant).

Same 1-D settling/consolidation ODE as ``thickener_profile_rk4.py`` but
integrated with a composite Simpson rule. Plant constants come from
``config.py``; the mean floc diameter is read from the JSON produced by
``floc_population_balance.py`` (run that first). Figure -> figures/.
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import config

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, os.pardir, "figures")

coneHeight = config.CONE_HEIGHT
cylinderHeight = config.CYLINDER_HEIGHT
Qufeed = config.FEED_FLOW
Qunderfl = config.UNDERFLOW
Fifeed = config.FEED_SOLID_FRACTION
psolid = config.SOLID_DENSITY
pfluid = config.LIQUOR_DENSITY
muliqour = config.LIQUOR_VISCOSITY

_diam_file = os.path.join(HERE, os.pardir, config.MEAN_DIAMETER_FILE)
try:
    Mean_diameter = json.load(open(_diam_file, encoding="utf-8"))["mean_floc_diameter_m"]
except FileNotFoundError:
    sys.exit(f"{_diam_file} not found - run src/floc_population_balance.py first")

print("mean floc diameter:", Mean_diameter)
## расчет геометрии сгустителя

# расчет площади для  высоты
i = 0.01
height = []  # Текущая высота сгустителя
D = 30
height_step = 0.3
# coneHeight = 1.35 #м высота конуса
# cylinderHeight = 1 # высота цилиндра или высота ниже основания питающего колодца
lim = coneHeight + cylinderHeight
print(lim)
while i < lim:
    height.append(i)
    i = i + height_step

Diameter = np.zeros(len(height))  # Диаметр цилиндрической части
Square = np.zeros(len(height))

j = 0
while j <= 7:
    if height[j] <= cylinderHeight:
        Diameter[j] = D
    else:
        Diameter[j] = D - (height[j] - cylinderHeight) * D / coneHeight
    Square[j] = (3.14 * Diameter[j] ** 2) / 4
    j = j + 1

# Qufeed=350 #Расход питающего потока
Qinj = 20  # Расход разбавления
Q = Qufeed + Qinj  # Объемный расход всего, что поступает в сгуститель

# Qunderfl=90
QL = Q - Qunderfl  # Расход жидкого из сгустителя
QR = Qunderfl

# Fifeed=0.0159 #Объемная концентраци тв в питающей пульпе
cfeed = Fifeed * Qufeed / (Qufeed + Qinj)  # Объемная доля тв в питании
# Qfloc=2 #Расход флокулянта %
Rowater = 1020  # Плотность воды
Cfloc_w = 0.005
# Gfloc=Qfloc*Rowater*Cfloc_w #Массовый расход флокулянта
# psolid=3200 #Плотность тв
Gtv = Qufeed * psolid * Fifeed  # Массовый расход тв в питающем потоке
# R=Gfloc/Gtv
# pfluid=1240 #Плотность жидкого
cfeed_wt = cfeed * psolid / (psolid * cfeed + (1 - cfeed) * pfluid)
g = 9.81
ccr = 0.03
sigma0 = 2
k = 6.5
# muliqour=0.0021 #Вязкость раствора
# dfloc=(708+0.1133*R-112.6*100*cfeed_wt)*0.54*0.00001
dfloc = Mean_diameter
pfluidrazbab = (pfluid * Qufeed * (1 - Fifeed) + Rowater * Qinj) / (
            Qufeed * (1 - Fifeed) + Qinj)  # Плотность жидксоти (алюминат+вода)
c_out = Q / Qunderfl * cfeed  # Доля тв на выходе
v = (dfloc ** 2 * g * (psolid - pfluidrazbab)) / (18 * muliqour)  # Скорость осаждения Стокса

## расчет пространственных скоростей
qf = np.zeros(len(height))
qtv_from = np.zeros(len(height))
q = np.zeros(len(height))
qR = np.zeros(len(height))
qL = np.zeros(len(height))
c1 = np.zeros(len(height))
fd = np.zeros(len(height))
ff = np.zeros(len(height))
k = 0
for k in range(len(height)):
    qf[k] = Q / Square[k] / 3600  # Скорость всего, что поступает в сгуститель
    qtv_from[k] = Qunderfl / Square[k] / 3600  # Скорость потока тв из сгустителя
    qL[k] = QL / Square[k] / 3600
    c1[k] = (qtv_from[k] - qL[k]) * cfeed / (qtv_from[k] + v)
    fd[k] = qtv_from[k] * c_out
    ff[k] = qf[k] * Fifeed
    k = k + 1

##значения для функции
dp = psolid - pfluidrazbab  # Разница плотностей тв и жидкого
Cmax = 1  # Максимальная концентрация
n = 87  # Безразмерный индекс стесненного осаждения
eps1 = 10e-4
h = -0.02  # Шаг модели
Fc = c_out
qtvtek = 0.01
x = 2.35
fbk = 0


## функция для расчета концентрации
def fbk__fun(c, Cmax, v, n):
    if (c < Cmax):
        Fun1 = v * c * (1 - c) ** n
    else:
        Fun1 = 0
    return Fun1


def dsigma2__fun(c, ccr, k, sigma0):  # функция для dsigma2
    if c > ccr:
        Fun2 = sigma0 * k / ccr * ((c / ccr) ** (k - 1))
    else:
        Fun2 = 0.00001
    return Fun2


def a_fun(c, fbk, dsigma2, dp, g):  # функция для расчета а
    if c == 0:
        Fun3 = 0
    else:
        Fun3 = fbk * dsigma2 / (dp * g * Fc)
    return Fun3


def FCC_fun(qtv_from, c, c_out, fbk, a, eps1):  # функция для расчета концентрации
    f1 = qtv_from * (c - c_out)
    f11 = fbk
    f2 = a
    f3 = eps1
    Fun4 = (f1 + f11) / (f2 + f3)
    return Fun4


# -------------------------------------------
# методом Симпсона

def fbk_fun(c, Cmax, v, n):
    if c.any() < Cmax:
        return v * c * (1 - c) ** n
    else:
        return 0


def dsigma2_fun(c, ccr, k, sigma0):
    if c.any() > ccr:
        return sigma0 * k / ccr * ((c / ccr) ** (k - 1))
    else:
        return 0.00001


def a_fun(c, fbk, dsigma2, dp, g):
    if c.any() == 0:
        return 0
    else:
        return fbk * dsigma2 / (dp * g * Fc)


def FCC_fun(qtv_from, c, c_out, fbk, a, eps1):
    f1 = qtv_from * (c - c_out)
    f11 = fbk
    f2 = a
    f3 = eps1
    return (f1 + f11) / (f2 + f3)


# Параметры
h = 0.01  # Шаг интегрирования
x_start = -2  # Начальная точка
x_end = 2  # Конечная точка
n_steps = int((x_end - x_start) / h) + 1  # Количество шагов


x = np.zeros(n_steps)
c = np.zeros(n_steps)


x[0] = x_start
c[0] = 0.01  # Начальная концентрация


for i in range(1, n_steps):
    # начений функции в начальной и конечной
    c_start = c[i - 1]
    c_end = FCC_fun(qtv_from, c[i - 1], c_out, fbk_fun(c[i - 1], Cmax, v, n),
                    a_fun(c[i - 1], fbk_fun(c[i - 1], Cmax, v, n), dsigma2_fun(c[i - 1], ccr, k, sigma0), dp, g), eps1)

    # Вычисление промежуточных значений функции
    c_mid = 0.5 * (c_start + c_end)
    c_mid_plus_h2 = c_mid + h / 2
    c_mid_minus_h2 = c_mid - h / 2

    # Вычисление интеграла методом Симпсона
    integral = (h / 6) * (c_start + c_end + 4 * c_mid + 2 * (
                FCC_fun(qtv_from, c_mid_plus_h2, c_out, fbk_fun(c_mid_plus_h2, Cmax, v, n),
                        a_fun(c_mid_plus_h2, fbk_fun(c_mid_plus_h2, Cmax, v, n),
                              dsigma2_fun(c_mid_plus_h2, ccr, k, sigma0), dp, g), eps1) + FCC_fun(qtv_from,
                                                                                                  c_mid_minus_h2, c_out,
                                                                                                  fbk_fun(
                                                                                                      c_mid_minus_h2,
                                                                                                      Cmax, v, n),
                                                                                                  a_fun(c_mid_minus_h2,
                                                                                                        fbk_fun(
                                                                                                            c_mid_minus_h2,
                                                                                                            Cmax, v, n),
                                                                                                        dsigma2_fun(
                                                                                                            c_mid_minus_h2,
                                                                                                            ccr, k,
                                                                                                            sigma0), dp,
                                                                                                        g), eps1)))
    # Обновление значения концентрации
    c[i] = c[i - 1] + integral[0]

    # Обновление значения высоты
    x[i] = x[i - 1] + h

# Вывод результатов
print("concentration:", c)
print("height:", x)

plt.plot(c * 100, x)
plt.ylim(bottom=0)
plt.xlabel('Концентрация, %')
plt.ylabel('Высота, см')
plt.gca().invert_yaxis()
plt.title('Распределение концентрации по высоте')
plt.grid()
os.makedirs(FIG_DIR, exist_ok=True)
plt.savefig(os.path.join(FIG_DIR, "concentration_profile_simpson.png"), dpi=120, bbox_inches="tight")
print("saved figure -> figures/concentration_profile_simpson.png")

