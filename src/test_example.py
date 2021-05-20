# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 02:52:18 2020

@author: kuzne
"""

from geam.GEAM import run
import matplotlib.pyplot as plt
import numpy as np

# Solve for 1-4.
t, u = run(0, 1, [1])
plt.figure()
plt.plot(t, u, 'o', MarkerSize=2)
B = np.exp(t) * np.tanh(.5)
#plt.plot(t, np.log((1 + B) / (1 - B)))
plt.plot(t, np.exp(-10 * t))
plt.legend(['num', 'exact'])
