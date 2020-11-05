# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 02:52:18 2020

@author: kuzne
"""

from solver_data import SolverData
from GEAM import run
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Default solver params are set in SolverData constructor.

lambd = 1e5
data = {}
fig, ax = plt.subplots()

# Solve for 1-4.
case = SolverData(lambd=lambd)
case.approx_order_1 = 1
case.approx_order_2 = 4
run(case)
data['1-4'] = case.__dict__
grids, errors = case.log_data()
switch = case.switch

# Plot data 1-4.
ax.plot(grids[switch:], errors[switch:],
         '-ok',
         MarkerFaceColor='k')

ax.plot(grids[: switch + 1], errors[: switch + 1],
         '-ok',
         MarkerFaceColor='w')

# Annotations 1-4.
ax.text(grids[switch], errors[switch] + 0.3, 'ERK1', FontSize=12)

# Solve for 4-4.
case = SolverData(lambd=lambd)
case.approx_order_1 = 4
case.approx_order_2 = 4
run(case)
data['4-4'] = case.__dict__
grids, errors = case.log_data()
switch = case.switch

# Plot data 4-4.
ax.plot(grids[switch:], errors[switch:],
         '-ok',
         MarkerFaceColor='k')

ax.plot(grids[: switch + 1], errors[: switch + 1],
         '-ok',
         MarkerFaceColor='w')

# Annotations 4-4.
ax.text(grids[-2], errors[-2] + 0.3, 'ERK4', FontSize=12)

# Labels.
ax.set_xlabel('lg N',
              FontSize=14,
              FontStyle='italic')

ax.set_ylabel(r'lg $\Delta$',
              FontSize=14,
              FontStyle='italic')

# Ticks.
ax.xaxis.set_minor_locator(AutoMinorLocator(n=2))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=2))
ax.tick_params(which='both', direction='in')
ax.tick_params(which='major', length=6, pad=5, labelsize='large')
ax.tick_params(which='minor', length=3)
