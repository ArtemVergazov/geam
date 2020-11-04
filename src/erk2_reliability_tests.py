# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 22:28:27 2020

@author: kuzne
"""

from solver_data import SolverData
from GEAM import run
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

# Default solver params are set in SolverData constructor.

p = 2
max_nu = 7
data = {}

fig, ax = plt.subplots()
for nu in range(1, max_nu + 1):
    lambd = 10**nu
    case = SolverData(lambd=lambd)
    case.approx_order_1 = p
    case.approx_order_2 = p
    run(case)
    
    text = str(nu)
    data[nu] = case.__dict__
    
    grids, errors = case.log_data()
    switch = case.switch
    
    # Plot data.
    ax.plot(grids[switch:], errors[switch:],
             '-ok',
             MarkerFaceColor='k')
    
    ax.plot(grids[: switch + 1], errors[: switch + 1],
             '-ok',
             MarkerFaceColor='w')
    
    # Labels.
    ax.set_xlabel('lg N',
                  FontSize=14,
                  FontStyle='italic')
    
    ax.set_ylabel(r'lg $\Delta$',
                  FontSize=14,
                  FontStyle='italic')
    
    # Ticks.
    ax.xaxis.set_ticks(np.arange(1., 5.))
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=2))
    
    ax.yaxis.set_ticks(np.arange(-8., 1., 2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=2))
    
    ax.tick_params(which='both', direction='in')
    ax.tick_params(which='major', length=6, pad=5, labelsize='large')
    ax.tick_params(which='minor', length=3)
    
    # Annotations.
    ax.text(4.2, errors[-1], text,
            FontSize=8,
            FontStyle='oblique',
            VerticalAlignment='center')
