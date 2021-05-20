# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 16:32:52 2019

@author: Artem
"""

from geam.rhs_func import rhs
from geam.solver_data import SolverData
import numpy as np


def erk1(u, h, F_n):
    return u + h * F_n


def erk2(u, h, F_n):
    a2 = 2 / 3
    b1, b2 = 1 / 4, 3 / 4
    w1 = F_n
    w2 = F(u + h * a2 * w1)
    return u + h * (b1 * w1 + b2 * w2)


def erk3(u, h, F_n):
    b1, b2, b3 = 2 / 9, 3 / 9, 4 / 9
    a21, a32 = 1 / 2, 3 / 4
    w1 = F_n
    w2 = F(u + h * a21 * w1)
    w3 = F(u + h * a32 * w2)
    return u + h * (b1 * w1 + b2 * w2 + b3 * w3)


def erk4(u, h, F_n):
    b1, b2, b3, b4 = 1 / 6, 1 / 3, 1 / 3, 1 / 6
    a21, a32, a43 = 1 / 2, 1 / 2, 1
    w1 = F_n
    w2 = F(u + h * a21 * w1)
    w3 = F(u + h * a32 * w2)
    w4 = F(u + h * a43 * w3)
    return u + h * (b1 * w1 + b2 * w2 + b3 * w3 + b4 * w4)


def erk(u, h, F_n, p):
    if p == 1:
        return erk1(u, h, F_n)
    elif p == 2:
        return erk2(u, h, F_n)
    elif p == 3:
        return erk3(u, h, F_n)
    elif p == 4:
        return erk4(u, h, F_n)
    else:
        raise ValueError('Invalid order of approximation p = {}!'.format(p))
        
        
def F(u):
    r = rhs(u)
    return r / np.linalg.norm(r)


def richardson(u1, u2, H, p):
    # H must be array of steps from previous grid.
    
    min_size = min(H.size, (len(u2) - 1) // 2)
    u1 = u1[:min_size]
    u2 = u2[::2][:min_size]
    H = H[:min_size]
    R = (u1 - u2) / (2**p - 1)
    
    # Left rectangles norm.
    R = np.linalg.norm(np.sqrt((R.T**2 * H).sum(1) / H.sum()))
    return R


def run_iterations(case):
    # Retrieve settings.
    u0 = case.u0
    u = u0
    Nmin = case.Nmin[-1]
    Nmax = case.Nmax[-1]
    if case.lengths == []:
        L_old = case.length_0
    else:
        L_old = case.lengths[-1][-1]
    if case.integrals == []:
        int_old = case.integral_0
    else:
        int_old = case.integrals[-1]
    p = case.approx_order_1

    # Temps.
    U = [u]
    F_last = F(u)
    kappa = 10.  # on the first step.
    kappas = []
    H = []  # steps
    L = [0]
    
    while True:
        # Extrapolate kappa from previous step.
        h = 1 / (Nmin / L_old + Nmax * kappa**.4 / int_old)
        H.append(h)
        L.append(L[-1] + h)

        # except ZeroDivisionError as e:
        #     print(e)
        #     return L, L_mas, H, U, integral, DL, kappas

        u_new = erk(u, h, F_last, p)
        U.append(u_new)

        F_new = F(u_new)  # rhs on next step
        
        # real kappa on current step
        kappa = np.linalg.norm(F_new - F_last) / h
        kappas.append(kappa)
        F_last = F_new

        u = u_new
        
        if u_new[0] >= case.T:
            break

    H = np.array(H)
    kappas = np.array(kappas)
    U = np.array(U)
    L = np.array(L)

    integral = (kappas**0.4 * H).sum()
    # DL = np.sqrt((((U[1:, 1] - UL(L[1:], lambd, u0))**2 +
    #                (U[1:, 0] - TL(L[1:], lambd, u0))**2) /
    #               (UL(L[1:], lambd, u0)**2 + 
    #                TL(L[1:], lambd, u0)**2) * H).sum() / H.sum())

    case.solutions.append(U)
    case.lengths.append(L)
    case.steps.append(H)
    case.integrals.append(integral)
    # case.errors.append(DL)
    
    case.Nmin.append(2 * Nmin)
    case.Nmax.append(2 * Nmax)
    
    case.grid_sizes.append(H.size)
    
    
def thicken(steps):
    '''
    Quasi-uniformly thicken a mesh
    that is given as a numpy vector of steps.
    '''
    N = steps.size
    steps_new = np.zeros(2 * N)
    
    if N == 1:
        steps_new[0] = steps[0] / 2
        steps_new[1] = steps[0] / 2
    else:
        steps_new[0] = steps[0] * steps[0]**.5 / (
            steps[0]**.5 + steps[1]**.5)
        steps_new[1] = steps[0] * steps[1]**.5 / (
            steps[0]**.5 + steps[1]**.5)
        steps_new[-2] = steps[-1] * steps[-2]**.5 / (
            steps[-2]**.5 + steps[-1]**.5)
        steps_new[-1] = steps[-1] * steps[-1]**.5 / (
            steps[-2]**.5 + steps[-1]**.5)

        if N > 2:
            for i in range(1, N - 1):
                steps_new[2 * i] = steps[i] * steps[i - 1]**.25 / (
                    steps[i - 1]**.25 + steps[i + 1]**.25)
                steps_new[2 * i + 1] = steps[i] * steps[i + 1]**.25 / (
                    steps[i - 1]**.25 + steps[i + 1]**.25)
    return steps_new


def stage1(case):
    print('Stage 1...')
    counter = 0
    dist = case.crit1 + 1
    while dist > case.crit1:
        '''Apply Euler scheme and double Nmin and Nmax
        until the mesh becomes adaptive.'''
        counter += 1

        run_iterations(case)

        # print('N =', case.steps[-1].size, '\n')

        if counter >= 2:
            dist = 0
            H_old = case.steps[-2]
            H = case.steps[-1]
            N = min(H_old.size, H.size // 2)
            for n in range(N):
                ksi = (H[2 * n] + H[2 * n + 1]) / H_old[n]
                dist = dist + (ksi**0.5 - ksi**(-0.5))**2 * H[2 * n]
            dist /= H.sum()
            dist **= .5

            case.rich_errors.append(richardson(
                case.solutions[-2],
                case.solutions[-1],
                H_old,
                case.approx_order_1))


def stage2(case):
    print('Stage 2...')
        
    k = 1
    while case.rich_errors[-1] > case.crit2 or k <= 3:        
        H = case.steps[-1]
        H_new = thicken(H)
        H_old = H
        H = H_new

        u = case.u0
        U = [u]
        L = [0]
        for h in H:
            L.append(L[-1] + h)
            u = erk(u, h, F(u), case.approx_order_2)
            U.append(u)
            
        L = np.array(L)
        U = np.array(U)
            
        case.grid_sizes.append(H.size)
        case.solutions.append(U)
        case.lengths.append(L)
        case.steps.append(H)
        if len(case.grid_sizes) > 1:
            case.rich_errors.append(richardson(
                case.solutions[-2],
                case.solutions[-1],
                H_old,
                case.approx_order_2))
        
        k += 1


def run(t0, T, u0):    
    case = SolverData(t0, T, u0)
    
    stage1(case)
    
    # Index of the last grid of stage 1.
    case.switch = len(case.grid_sizes) - 1
    print('''Switch between the strategies is done after '''
          f'''{case.switch + 1}-th grid''')

    stage2(case)
    
    case.plot_error()
    
    t = case.solutions[-1][:, 0]
    u = case.solutions[-1][:, 1:]
    return t, u
                        
