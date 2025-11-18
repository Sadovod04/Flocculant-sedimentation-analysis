import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.stats as stats
import pandas as pd
import math
from scipy.interpolate import interp1d
import scipy.integrate
import sqlite3

# Устанавливаем соединение с базой данных
connection = sqlite3.connect('tasks.db')
cursor = connection.cursor()

# Создаем таблицу Tasks
cursor.execute('''
CREATE TABLE IF NOT EXISTS Tasks (
id INTEGER PRIMARY KEY,
Liquor_viscosity INTEGER ,
Particle_density INTEGER ,
surface_area INTEGER ,
Flocculant_dosage INTEGER ,
Feed_solid_concentration INTEGER ,
Liquor_density INTEGER ,
Pipe_flow_rate INTEGER ,
Well_inner_diameter INTEGER ,
Well_height INTEGER ,
Simulation_run_time INTEGER ,
Mean_diameter REAL 
)
''')

# Задаем данные распределения частиц по размерам (размеры в микрометрах и кумулятивные доли)
df = pd.read_csv("PSD_disc.csv", sep=";")
particle_sizes = df['mkm'].to_list()
percentages = df['%'].to_list()

def plot_cumulative_curve(particle_sizes, percentages):
    cumulative_percentages = np.cumsum(percentages)
    
    plt.plot(particle_sizes, cumulative_percentages, marker='o')
    plt.xlabel('Размер частиц')
    plt.ylabel('Кумулятивный процент')
    plt.title('Кумулятивная кривая размеров частиц')
    plt.grid(True)
    plt.show()
    return cumulative_percentages

plot_cumulative_curve(particle_sizes, percentages)


def cumulative_to_gaussian(particle_sizes, percentages):
    mean = np.average(particle_sizes, weights=percentages)
    variance = np.average((particle_sizes - mean)**2, weights=percentages)
    std_dev = np.sqrt(variance)

    gaussian_values = []
    for size, percentage in zip(particle_sizes, percentages):
        z_score = (size - mean) / std_dev
        gaussian_values.append((z_score, percentage))

    return gaussian_values, mean, std_dev


gaussian_values, mean, std_dev = cumulative_to_gaussian(particle_sizes, percentages)

# Вывод нормального гауссовского распределения
x = np.linspace(mean - 3*std_dev, mean + 3*std_dev, 20)
y_N = np.zeros_like(x)
for z_score, percentage in gaussian_values:
    y_N += percentage * np.exp(-0.5 * (x - mean)**2 / std_dev**2) / (std_dev * np.sqrt(2*np.pi))

plt.plot(x, y_N, 'r', label='Исходные данные')
plt.xlabel('Значение')
plt.ylabel('Плотность вероятности')
plt.title('Нормальное гауссовское распределение')
plt.show()
print(len(x))
print(len(y_N))


def odefunfik1(t, N):
    dNdt = np.zeros(20)

    file1 = np.zeros(20)
    file2 = np.zeros(20)
    file3 = np.zeros(20)
    file4 = np.zeros(20)
    file5 = np.zeros(20)
    file6 = np.zeros(20)

    Df = 2.4 #Fractal dimension (Df) NM
    Kp=1 #Packing factor (kp) NM
    V=20 #mean pipe flow velocity
    D=0.0254 #Pipe ID
    mu0 = 0.001 #Liquor viscosity
    k1=0.075
    k2=94.7
    k3=0.691
    k4=1220
    L=0.1
    mf=0.2 #mass floc
    As=2.971 #surface area of solid (floc)
    Teta=0.9#floc degradation
    ps=3200 #Particle density of the solid
    pl=1240 #density of the liquid
    fi=100
    fim=0.65
    k=2
    dp = 0.22 * 10 ** (-6)
    mu=mu0*(1-(fi/fim)*((100*10**(-6))/dp)**(3-Df))**(-k)
    pf=ps*fi+pl*(1-fi)
    Re=(D*V*pf)/mu
    f=0.0791/(Re**0.25)
    epsi=(2*f*V**3)/D
    G=(epsi*pf/mu)**0.5
    M=1 - (math.e) ** ((-k1 *f**0.5* L) / D)
    alpha=M
    Qf = mf / As * M * (1 - Teta)

    tzv = 0.1
    kb = 1
    G = 0.1
    nu0 = 10 ** (- 3)
    fit0 = 0.6
    q = 1.3
    ddd = 0.000015
    Df = 2.30
    Pi = 3.14
    a1 = 0
    a2 = 0
    a3 = 0
    a4 = 0

    fitot = 0.0
    for i in np.arange(0, 20, 1):
        dc = dp * ((2 ** i) ** (1 / Df))
    dc = np.zeros(20)
    for j in np.arange(0, 20, 1):
        dc[j] = dp * ((2 ** j) ** (1 / Df))
        fragkerns = k2 *epsi**k3*mu*(dp * Kp ** (-1 / 3) * (dc[j] / dp) ** (3 / Df))/(2*Qf)
        if j == 1:
            fragdistfun = 2
        else:
            fragdistfun = 0
        a4 = fragdistfun * fragkerns * N[j] + a4
    
    file6[1] = a4

    bet = np.zeros(20)
    for j in np.arange(0, 20, 1):
        bet[j] = 1.294 * G * alpha * ((2 ** (1 - 1)) ** (1 / Df) + (2 ** j) ** (1 / Df)) ** 3
        a3 = bet[j] * N[j] + a3
    file4[1] = a3 * N[1]
    dc = dp * ((2 ** (1 - 1)) ** (1 / Df))
    
    fragkerns = k2 *epsi**k3*mu*(dp * Kp ** (-1 / 3) * (dc/ dp) ** (3 / Df))/(2*Qf)
    A1=dp * Kp ** (-1 / 3) * (dc/ dp) ** (3 / Df)
    file5[1] = N[1] * fragkerns
    dNdt[1] = file6[1] - file4[1] - file5[1]
    
    file1[1] = 0
    file2[1] = 0
    file3[1] = 0
    a3 = 0
    a4 = 0
    bet = 1.294 * G * alpha * ((2 ** (2 - 1)) ** (1 / Df) + (2 ** (1 - 1)) ** (1 / Df)) ** 3
    a2 = 2 ** (1 - 2) * bet * N[1]
    bet = np.zeros(20)
    file3[2] = a2 * N[2]
    for j in np.arange(1, 20, 1):
        bet[j] =1.294 * G * alpha * ((2 ** (2 - 1)) ** (1 / Df) + (2 ** j) ** (1 / Df)) ** 3
        a3 = bet[j] * N[j] + a3
    file4[2] = N[2] * a3
    dc = np.zeros(20)
    fragkerns = np.zeros(20)
    for j in np.arange(1, 20, 1):
        dc[j] = dp * ((2 ** j) ** (1 / Df))
        fragkerns[j] = k2 *epsi**k3*mu*(dp * Kp ** (-1 / 3) * (dc[j]/ dp) ** (3 / Df))/(2*Qf)
        if j == 2:
            fragdistfun = 2
        else:
            fragdistfun = 0
        a4 = fragdistfun * fragkerns[j] * N[j] + a4
    
    file6[2] = a4
    dc = dp * ((2 ** (2 - 1)) ** (1 / Df))
    fragkerns =k2 *epsi**k3*mu*(dp * Kp ** (-1 / 3) * (dc/ dp) ** (3 / Df))/(2*Qf)
    A2=dp * Kp ** (-1 / 3) * (dc/ dp) ** (3 / Df)
    file5[2] = fragkerns * N[2]
    bet = 1.294 * G * alpha *  ((2 ** (1 - 1)) ** (1 / Df) + (2 ** (1 - 1)) ** (1 / Df)) ** 3
    file2[2] = (1 / 2) * bet * (N[1]) ** 2
    file1[2] = 0
    dNdt[2] = file2[2] - file3[2] - file4[2] - file5[2] + file6[2]
    a2 = 0
    a3 = 0
    a4 = 0
    bet = np.zeros(20)
    for i in np.arange(2, 20, 1):
        for j in np.arange(1, i - 2, 1):
            bet[j] = 1.294 * G * alpha *  ((2 ** (i - 1)) ** (1 / Df) + (2 ** j) ** (1 / Df)) ** 3
            if j - i < 0:
                a1 = 1 / (2 ** (i - j + 1)) * bet[j] * N[j] * N[i - 1] + a1
            else:
                a1 = 2 ** (j - i + 1) * bet[j] * N[j] * N[i - 1] + a1
        file1[i] = a1
        bet = 1.294 * G * alpha *  ((2 ** (i - 1)) ** (1 / Df) + (2 ** (i - 1)) ** (1 / Df)) ** 3
        file2[i] = (1 / 2) * bet * (N[i - 1]) ** 2
        bet = np.zeros(20)
        for j in np.arange(0, i+1, 1):
            bet[j] = 1.294 * G * alpha *  ((2 ** i) ** (1 / Df) + (2 ** j) ** (1 / Df)) ** 3
            if j - i < 0:
                a2 = 1 / (2 ** (i - j)) * bet[j] * N[j] + a2
            else:
                a2 = 2 ** (j - i) * bet[j] * N[j] + a2
        file3[i] = N[i] * a2
        for j in np.arange(i, 20, 1):
            bet[j] = 1.294 * G * alpha *  ((2 ** i) ** (1 / Df) + (2 ** j) ** (1 / Df)) ** 3
            a3 = bet[j] * N[j] + a3
        file4[i] = a3 * N[i]
        dc = np.zeros(20)
        fragkerns = np.zeros(20)
        for j in np.arange(i, 20, 1):
            dc[j] = dp * ((2 ** j) ** (1 / Df))

            fragkerns[j] = k2*epsi**k3*mu*(dp * Kp ** (-1 / 3) * (dc[j]/ dp) ** (3 / Df))/(2*Qf)
            if j == i + 1:
                fragdistfun = 2
            else:
                fragdistfun = 0
            a4 = fragdistfun * fragkerns[j] * N[j] + a4
            A3=dp * Kp ** (-1 / 3) * (dc[j]/ dp) ** (3 / Df)
        file6[i] = a4
        dc[i] = dp * ((2 ** (i - 1)) ** (1 / Df))

        fragkerns[i] = k2 *epsi**k3*mu*(dp * Kp ** (-1 / 3) * (dc[i]/ dp) ** (3 / Df))/(2*Qf)
        A4=dp * Kp ** (-1 / 3) * (dc[i]/ dp) ** (3 / Df)
        file5[i] = fragkerns[i] * N[i]
        a1 = 0
        a2 = 0
        a3 = 0
        a4 = 0
        dNdt[i] = file1[i] + file2[i] - file3[i] - file4[i] - file5[i] + file6[i]

    return dNdt

MM = y_N
N = scipy.integrate.odeint(odefunfik1, MM[:20], np.arange(0, 20), tfirst=True)
chisl=np.zeros(20)
znam=np.zeros(20)
mean_diameter=np.zeros(20)
for i in range(1,20):
    plt.plot(N[i], label='График ' + str(i))
    # Рассчитываем средний диаметр частиц по де Брукеру
    mean_diameter[i] = np.mean(N[i])
    chisl[i]=mean_diameter[i]**4
    znam[i]=mean_diameter[i]**3
    df=pd.DataFrame(N)
    df.to_excel('./ts.xlsx')
plt.title('Распределение частиц по размерам', fontsize=14)
plt.xlim(0)
plt.ylim(0)
plt.show()
mean=round(sum(chisl)/sum(znam),5)
print("Средний диаметр частиц:", mean)
cursor.execute('INSERT INTO Tasks (Mean_diameter) VALUES (?)', (mean,))
connection.commit()
# Закрываем соединение
connection.close()