# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 00:36:11 2021

@author: kuzne
"""

import numpy as np

def rhs(u):
    '''
    Right-hand side function of differential equation.
    
    Input
    -------
    u - vector, argument to the function
        
    Returns
    -------
    Vector rhs(u)

    '''
    lambda_ = 1e6
    return np.array([1., np.sinh(lambda_ * u[1])])
