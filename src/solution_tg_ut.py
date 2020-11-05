# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 15:22:39 2020

@author: kuzne
"""

import numpy as np
import matplotlib.pyplot as plt

lambd = 10
maxu = np.pi / 2 / lambd
u0 = maxu / 10
assert u0 < maxu, maxu
t_inf = -1 / lambd * np.log(np.sin(lambd * u0))

# u(t)
t = np.linspace(-0.25 * t_inf, t_inf, 150)
u = lambda t: \
    1 / lambd * np.arcsin(np.exp(lambd * t) * np.sin(lambd * u0))

ut = u(t)
plt.plot(t, ut)
plt.xlim([t[0], 1.1 * t_inf])
plt.ylim([0, 1.2 * u(0.99 * t_inf)])
plt.plot(t_inf, u(t_inf), 'o', MarkerSize=8)

# u(l) & t(l)
L_max = -1 / lambd * np.log(np.tan(lambd * u0 / 2))
L = np.linspace(-0.25 * L_max, L_max, 150)
u = lambda L: \
    1 / lambd * 2 * np.arctan(np.exp(lambd * L) * np.tan(lambd * u0 / 2))
t = lambda L: \
    1 / lambd * np.log(
        np.sin(2 * np.arctan(np.exp(lambd * L) *
                             np.tan(lambd * u0 / 2))) / np.sin(lambd * u0))
plt.figure()
plt.plot(L, u(L))
plt.plot(L[L >= 0], t(L[L >= 0]))
plt.plot(L_max, u(L_max), 'ok', MarkerSize=8)
plt.plot(L_max, t(L_max), 'ok', MarkerSize=8)
plt.vlines(L_max, 0, max(u(L_max), t(L_max)) * 1.2)
