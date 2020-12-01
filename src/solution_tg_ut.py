# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 15:22:39 2020

@author: kuzne
"""

import numpy as np
import matplotlib.pyplot as plt

lambd = 10
umax = np.pi / 2 / lambd
u0 = umax / 10
tmax = -1 / lambd * np.log(np.sin(lambd * u0))

# u(t)
t = np.linspace(-0.2 * tmax, tmax, 500)
u = lambda t: 1 / lambd * np.arcsin(np.exp(lambd * t) * np.sin(lambd * u0))

fig, ax = plt.subplots()
for loc in ['bottom', 'left', 'top', 'right']:
    ax.spines[loc].set_visible(False)
for loc in ['bottom', 'left']:
    ax.spines[loc].set_position('zero')

ax.plot(t, u(t), '-k', LineWidth=1.5)
ax.plot(0, u0, 'ok', MarkerFaceColor='w')
ax.plot(tmax, umax, 'ok') 

ax.set_xticks([0, tmax])
ax.set_xticklabels(['0', r'$t_*$'])
ax.set_yticks([u0, umax])
ax.set_yticklabels([r'$u_0$', r'$u_*$'], va='bottom')
ax.tick_params(length=0)

ax.set_xlim(right=1.2 * tmax)
ax.set_ylim(top=1.2 * umax, bottom=-.02 * umax)

left, right = ax.get_xlim()
bottom, top = ax.get_ylim()
ax.arrow(left, 0, right - left, 0,
         length_includes_head=True,
         color='k',
         width=.0001,
         head_width=.004,
         head_length=.004,
         overhang=1)
ax.arrow(0, 0, 0, top,
         length_includes_head=True,
         color='k',
         width=.0001,
         head_width=.004,
         head_length=.004,
         overhang=1)

ax.vlines(tmax, 0, top, colors='k', LineWidth=.8)

ax.set_xlabel(r'$t$', FontSize=14, position=(1, 0))
ax.set_ylabel(r'$u$', FontSize=14, position=(3, .95), rotation='horizontal')
