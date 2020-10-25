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
        lambda_:
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
    
    def __init__(self, lambda_):
        '''
        Set solver settings before run.

        Parameters
        ----------
        lambda_ : float
            Stiffness of the problem.
        '''
        # Lambda and dependent properties.
        self.lambda_ = lambda_
        self.u0 = np.array([0., 1 / lambda_ * np.arcsinh((lambda_ - (
            lambda_**2 - 4)**0.5) / 2)])
        self.pole = (-1 / lambda_) * np.log(np.sin(lambda_ * self.u0[1]))
        self.length_extr = 1 / lambda_ * \
            np.log(1 / np.sinh(lambda_ * self.u0[1]))
        
        # Scalars known before the run. May be customized by user.
        self.length_0 = 1.
        self.integral_0 = 1.
        self.crit1 = 1e-1
        self.num2 = 1e5
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

