# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 04:06:15 2020

@author: kuzne
"""

import numpy as np


class SolverData:
    '''
    Solver input and output for single lambda.
    
    Scalar instance fields:
        lambd:
            stiffness of the problem
        u0:
            initial value for t and u
        pole:
            pole of the derivative
        length_extr:
            values of l with max kappa
        crit1:
            stop criterion for stage 1
        num2:
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
        errors:
            errors w.r.t. exact solution in L
        rich_errors:
            Richardson errors
        integrals:
            numerical values of integrals used in step formula
        Nmin:
            parameters Nmin
        Nmax:
            parameters Nmax
            
    '''
    
    def __init__(self, lambd=10):
        '''
        Set solver settings before the run.

        Parameters
        ----------
        lambd : float
            Stiffness of the problem.
        '''
        # Lambda and dependent properties.
        self.lambd = lambd
        self.u0 = np.array([0., 1 / lambd * np.arcsinh((lambd - (
            lambd**2 - 4)**0.5) / 2)])
        self.pole = -1 / lambd * np.log(np.tanh(lambd * self.u0[1] / 2))
        self.length_extr = 1 / lambd * \
            np.log(1 / np.sinh(lambd * self.u0[1]))
        
        # Scalars known before the run. May be customized by user.
        self.length_0 = 1.
        self.integral_0 = 1.
        self.crit1 = 1e-1
        self.num2 = 1e4
        self.approx_order_1 = 1
        self.approx_order_2 = 1
        
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
        grids = np.array(self.grid_sizes)
        errors = np.array(self.errors)
        
        nan_grids = np.isnan(grids)
        nan_errors = np.isnan(errors)
        nan = np.logical_or(nan_grids, nan_errors)
        not_nan = np.logical_not(nan)
        
        grids = grids[not_nan]
        errors = errors[not_nan]
        
        assert (grids > 0.).all(), 'Grid sizes'
        assert (errors > 0.).all(), ''
        return np.log10(self.grid_sizes), np.log10(self.errors)

