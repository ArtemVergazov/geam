# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 15:22:39 2020

@author: kuzne
"""

import numpy as np
import matplotlib.pyplot as plt

lambd = 10
umax = np.pi / 2 / lambd
u0 = umax / 4
tmax = -1 / lambd * np.log(np.sin(lambd * u0))
Lmax = 1 / lambd * np.log(
    np.tan(lambd * umax / 2) /
    np.tan(lambd * u0 / 2))
L = np.linspace(-0.2 * Lmax, Lmax, 500)

u = lambda L: 2 / lambd * np.arctan(np.exp(lambd * L) * np.tan(lambd * u0 / 2))
t = lambda L: 1 / lambd * np.log(
    np.sin(2 * np.arctan(np.exp(lambd * L) * np.tan(lambd * u0 / 2))) /
    np.sin(lambd * u0))

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

ax.plot(Lmax, tmax, 'ok')
ax.plot(Lmax, umax, 'ok')
ax.vlines(Lmax, 0, top, color='k', LineWidth=.8)

ax.text(.95 * Lmax, .95 * u(Lmax), '$u(l)$',
        FontSize=14,
        ha='right',
        va='bottom')
ax.text(Lmax, .99 * tmax, '$t(l)$',
        FontSize=14,
        ha='right',
        va='bottom')

ax.set_xlabel('$l$', FontSize=14, position=(1, 0))
