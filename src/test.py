# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 22:40:43 2020

@author: kuzne
"""

from solver_data import SolverData
from GEAM import run
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

# Default solver params are set in SolverData constructor.

fig, ax = plt.subplots()
ax.plot([1, 2, 3, 5, 8])
ax.set_xticks([3])
ax.set_xticklabels(['3'], FontSize=22)
