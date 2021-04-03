# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 04:06:15 2020

@author: kuzne
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


class SolverData:
    '''
    Storage for solver input and output.
    
    Settings:
        u0:
            vector of initial values for t and u
        T:
            final value for t
        crit1:
            stop criterion for stage 1
        crit2:
            stop criterion for stage 2
        switch:
            number of grid on which switch between stages took place
        approx_order_1:
            stage 1 approximation order
        approx_order_2:
            stage 2 approximation order
        length_0:
            initial arc length for first grid
        integral_0:
            initial value of integral for first grid

    Lists of ___ for each grid:
        grid_sizes:
            numbers of steps
        solutions:
             2D solutions
        lengths:
            1D L values
        steps:
            1D steps arrays
        rich_errors:
            Richardson errors
        integrals:
            numerical values of integrals used in step formula
        Nmin:
            parameters Nmin
        Nmax:
            parameters Nmax
            
    '''
    
    def __init__(self, t0, T, u0):
        '''
        Set solver settings before the run.

        Parameters
        ----------
        t0 : float
            initial value for t
        T : float
            final value for t
        u0 : float
            initial value for u
        '''
        self.u0 = np.array([t0, u0])
        self.T = T
        
        # Scalars known before the run. May be customized.
        self.length_0 = 1.
        self.integral_0 = 1.
        self.crit1 = 1e-2
        self.crit2 = 1e-4
        self.approx_order_1 = 1
        self.approx_order_2 = 4
        
        # Values dependent on grid. Set empty before the run.
        self.grid_sizes = []
        self.solutions = []
        self.lengths = []
        self.steps = []
        self.errors = []
        self.rich_errors = []
        self.integrals = []
        self.Nmin = [6]
        self.Nmax = [20]

    def log_data(self):
        grids = np.array(self.grid_sizes)[1:]
        errors = np.array(self.rich_errors)
        
        nan_grids = np.isnan(grids)
        nan_errors = np.isnan(errors)
        nan = np.logical_or(nan_grids, nan_errors)
        not_nan = np.logical_not(nan)
        
        grids = grids[not_nan]
        errors = errors[not_nan]
        
        assert (grids > 0.).all(), 'Error in grid sizes'
        assert (errors > 0.).all(), 'Errors are < 0'
        return np.log10(self.grid_sizes[1:]), np.log10(self.rich_errors)
    
    def plot_error(self):
        grids, errors = self.log_data()
        switch = self.switch - 1  # we do not include the first grid
        
        fig, ax = plt.subplots()
        ax.plot(grids[switch:], errors[switch:], 
                '-ok',
                MarkerFaceColor='k')
        ax.plot(grids[: switch + 1], errors[: switch + 1],
                '-ok',
                MarkerFaceColor='w')
        
        # Annotations.
        scheme1 = 'ERK' + str(self.approx_order_1)
        scheme2 = 'ERK' + str(self.approx_order_2)
        ax.text(grids[switch], errors[switch] + 0.3, scheme1, FontSize=12)
        ax.text(grids[-2], errors[-2] + 0.3, scheme2, FontSize=12)
        
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
