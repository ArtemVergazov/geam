# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 15:22:39 2020

@author: kuzne
"""

import numpy as np
import matplotlib.pyplot as plt

lambd = 10
u0 = 4 / lambd * np.arcsinh((lambd - (lambd**2 - 4)**0.5) / 2)
t_star = -1 / lambd * np.log(np.tanh(lambd * u0 / 2))
tmax = .99 * t_star

# u(t)
t = np.linspace(-0.2 * tmax, tmax, 500)
B = lambda t: np.exp(lambd * t) * np.tanh(lambd * u0 / 2)
u = lambda t: 1 / lambd * np.log((1 + B(t)) / (1 - B(t)))

fig, ax = plt.subplots()
for loc in ['bottom', 'left', 'top', 'right']:
    ax.spines[loc].set_visible(False)
for loc in ['bottom', 'left']:
    ax.spines[loc].set_position('zero')

ax.plot(t, u(t), '-k', LineWidth=1.5)
ax.plot(0, u0, 'ok', MarkerFaceColor='w')    

ax.set_xticks([0, t_star])
ax.set_xticklabels(['0', r'$t_*$'])
ax.set_yticks([u0])
ax.set_yticklabels([r'$u_0$'], va='bottom')
ax.tick_params(length=0)

ax.set_xlim(right=1.2 * tmax)

left, right = ax.get_xlim()
bottom, top = ax.get_ylim()
ax.arrow(left, 0, right - left, 0,
         length_includes_head=True,
         color='k',
         width=.0001,
         head_width=.012,
         head_length=.004,
         overhang=1)
ax.arrow(0, 0, 0, top,
         length_includes_head=True,
         color='k',
         width=.0001,
         head_width=.004,
         head_length=.01,
         overhang=1)

ax.vlines(t_star, 0, top, colors='k', LineWidth=.8)

ax.set_xlabel(r'$t$', FontSize=14, position=(1, 0))
ax.set_ylabel(r'$u$', FontSize=14, position=(3, .9), rotation='horizontal')
