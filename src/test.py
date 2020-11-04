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

p = 2
lambd = 10**7
case = SolverData(lambd=lambd)
case.approx_order_1 = p
case.approx_order_2 = p
run(case)
data = {}
data[lambd] = case.__dict__
