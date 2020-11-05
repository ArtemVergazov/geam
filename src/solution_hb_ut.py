# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 15:22:39 2020

@author: kuzne
"""

import numpy as np
import matplotlib.pyplot as plt

lambd = 10
u0 = 1 / lambd * np.arcsinh((lambd - (lambd**2 - 4)**0.5) / 2)
umax = 1 / lambd * np.arcsinh((lambd + (lambd**2 - 4)**0.5) / 2)
tmax = 1 / lambd * np.log(np.tanh(lambd * umax / 2) / np.tanh(lambd * u0 / 2))

# u(t)
t = np.linspace(-0.2 * tmax, tmax, 500)
B = lambda t: np.exp(lambd * t) * np.tanh(lambd * u0 / 2)
u = lambda t: 1 / lambd * np.log((1 + B(t)) / (1 - B(t)))

fig, ax = plt.subplots()
ax.plot(t, u(t), '-k')
ax.set_xticks([])
ax.set_yticks([])
ax.spines['left'].set_position('zero')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
xlim = ax.get_xlim()[1]
ax.arrow(xlim, 0, .01 * xlim, 0)
ax.set_xlim(right = 1.2 * xlim)

ax.axvline(x=tmax, Color='k', LineWidth=1)
