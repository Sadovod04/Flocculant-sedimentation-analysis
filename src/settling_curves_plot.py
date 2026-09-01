"""Overlay the measured batch-settling curve with two model predictions.

The three series were sampled on different time grids; the model curves are
linearly resampled onto the measurement time axis before plotting. The figure
is written to ``figures/settling_curves.png``.
"""
x=[0, 1, 2, 2.42, 3.15, 4.15, 5, 5.55, 6, 7, 8, 9, 10,
              10.5, 11, 12, 12.5, 13, 14, 14.5, 15, 15.5, 16, 16.5, 17, 18, 19, 20, 20.3,
              21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 44, 46, 48, 50]
y = [117.6,116.2,115.2,114.7,113.7,112.5,110.3,108.3,107.1,106,104.1,102.1,100,98.7,97.9, 97.3,95.9,95.1,94.4,93.6,92.7,92.1,91.8,91.3,90.7,90.3,89.3,88.3,87.6,86.9,86.2,85.6,84.9,84.5, 84,83.4,83,82.7,82.2,81.9,81.4,81.1,80.8,80.5,80.2,80,79.7,79.6,79.4,79,78.7,78.5,78.3,78.3]
y_2=[120.47,
119.4213815,
116.0507537,
114.9047862,
113.4616506,
111.852938,
109.907461,
108.4403535,
107.5647478,
106.884907,
105.4733543,
104.1730185,
102.7599366,
101.8128944,
101.2605519,
100.5387812,
99.66984067,
99.1580491,
98.65399311,
97.66628729,
97.18145633,
96.70197628,
96.22742781,
95.75744466,
95.29170677,
94.82993427,
93.91733648,
93.01803483,
92.13078588,
91.25462887,
90.3888174,
89.53276755,
88.68601872,
87.84820383,
87.01902687,
86.19824579,
85.38565959,
84.58109857,
83.78441684,
82.99548677,
82.21419469,
81.44043768,
80.67412113,
79.91515685,
79.16346171,
78.41895654,
77.68156534,
76.95121467,
76.22783315,
74.80170045,
73.40262639,
72.0300879,
70.68357612,
69.36259436,
]
y_1=[
111.31,
110.9937747,
109.842018,
109.3950715,
108.7848046,
108.0344819,
107.0151136,
106.1562873,
105.6042871,
105.1548065,
104.1628635,
103.1803675,
102.042719,
101.2433574,
100.7648698,
101.272751,
99.34306689,
98.8736498,
98.40647337,
97.47879971,
97.01828126,
96.55996098,
96.10382838,
95.64987303,
95.19808451,
94.74845251,
93.85561689,
92.97128439,
92.09537405,
91.22780565,
90.36849974,
89.51737764,
88.67436139,
87.83937382,
87.01233844,
86.19317954,
85.38182208,
84.57819179,
83.78221505,
82.99381899,
82.21293141,
81.43948079,
80.67339632,
79.91460783,
79.16304584,
78.41864154,
77.68132674,
76.95103394,
76.22769625,
74.8016219,
73.40258132,
72.03006204,
70.68356129,
69.36258585]
import os
import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "figures")


def _against_time(series):
    """Model arrays were sampled on a finer grid than the measurements;
    interpolate them onto ``x`` so every curve shares one time axis."""
    idx = [i * (len(x) - 1) / (len(series) - 1) for i in range(len(series))]
    return [
        _interp(t, idx, series) for t in range(len(x))
    ]


def _interp(pos, idx, series):
    for j in range(1, len(idx)):
        if idx[j] >= pos:
            f = (pos - idx[j - 1]) / (idx[j] - idx[j - 1])
            return series[j - 1] + f * (series[j] - series[j - 1])
    return series[-1]


plt.plot(x, y, "ko-", label="Эксперимент", markersize=3)
plt.plot(x, _against_time(y_1), "b-", label="Модель 1")
plt.plot(x, _against_time(y_2), "r-", label="Модель 2")
plt.ylabel("Высота раздела, см")
plt.xlabel("Время, мин")
plt.legend()
os.makedirs(FIG_DIR, exist_ok=True)
plt.savefig(os.path.join(FIG_DIR, "settling_curves.png"), dpi=120, bbox_inches="tight")
print("saved figure -> figures/settling_curves.png")