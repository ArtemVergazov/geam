# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 15:22:39 2020

@author: kuzne
"""

import numpy as np
import matplotlib.pyplot as plt

lambd = 10
u0 = 5.6 / lambd * np.arcsinh((lambd - (lambd**2 - 4)**0.5) / 2)
t_star = -1 / lambd * np.log(np.tanh(lambd * u0 / 2))
tmax = .99 * t_star
Bmax = np.exp(lambd * tmax) * np.tanh(lambd * u0 / 2)
Lmax = 1 / lambd * np.log(
    np.sinh(np.log(1 + Bmax) - np.log(1 - Bmax)) /
    np.sinh(lambd * u0))
L = np.linspace(-0.2 * Lmax, Lmax, 500)

A = lambda L: np.exp(lambd * L) * np.sinh(lambd * u0)
u = lambda L: 1 / lambd * np.log(A(L) + np.sqrt(A(L)**2 + 1))
t = lambda L: 1 / lambd * np.log(
    np.tanh(1 / 2 * np.log(A(L) + np.sqrt(A(L)**2 + 1))) /
    np.tanh(lambd * u0 / 2))

fig, ax = plt.subplots()
for loc in ['bottom', 'left', 'top', 'right']:
    ax.spines[loc].set_visible(False)
for loc in ['bottom', 'left']:
    ax.spines[loc].set_position('zero')

ax.plot(L, u(L), '-k', LineWidth=1.5)
ind = t(L) >= -.15 * tmax
ax.plot(L[ind], t(L)[ind], '-k', LineWidth=1.5)
ax.plot(0, u0, 'ok', MarkerFaceColor='w')

ax.set_xticks([0])
ax.set_xticklabels(['0'])
ax.set_yticks([u0])
ax.set_yticklabels([r'$u_0$'], va='bottom')
ax.tick_params(length=0)

#ax.set_ylim(top=2.5 * tmax)

left, right = ax.get_xlim()
bottom, top = ax.get_ylim()
ax.arrow(left, 0, right - left, 0,
         length_includes_head=True,
         color='k',
         width=.0001,
         head_width=.01,
         head_length=.008,
         overhang=1)
ax.arrow(0, 0, 0, top,
         length_includes_head=True,
         color='k',
         width=.0001,
         head_width=.01,
         head_length=.008,
         overhang=1)

ax.axhline(t_star, color='k', LineWidth=.8)
k = 1
b = 1 / lambd * np.log(2 * np.sinh(lambd * u0))
y = k * L + b
ax.plot(L[y >= 0], y[y >= 0], '-k', LineWidth=.8)

ax.text(.95 * Lmax, .95 * u(Lmax), '$u(l)$',
        FontSize=14,
        ha='right',
        va='bottom')
ax.text(Lmax, .99 * t_star, '$t(l)$',
        FontSize=14,
        ha='right',
        va='top')
ax.text(1.05 * Lmax, t_star, '$t_*$',
        FontSize=14,
        ha='right',
        va='bottom')

ax.set_xlabel('$l$', FontSize=14, position=(1, 0))
