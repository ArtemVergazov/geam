# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 16:32:52 2019

@author: Artem
"""

import numpy as np
import warnings


def erk1(u, h, F_n, lambda_):
    return u + h * F_n


def erk2(u, h, F_n, lambda_):
    a2 = 2 / 3
    b1, b2 = 1 / 4, 3 / 4
    w1 = F_n
    w2 = F(u + h * a2 * w1, lambda_)
    return u + h * (b1 * w1 + b2 * w2)


def erk3(u, h, F_n, lambda_):
    b1, b2, b3 = 2 / 9, 3 / 9, 4 / 9
    a21, a32 = 1 / 2, 3 / 4
    w1 = F_n
    w2 = F(u + h * a21 * w1, lambda_)
    w3 = F(u + h * a32 * w2, lambda_)
    return u + h * (b1 * w1 + b2 * w2 + b3 * w3)


def erk4(u, h, F_n, lambda_):
    b1, b2, b3, b4 = 1 / 6, 1 / 3, 1 / 3, 1 / 6
    a21, a32, a43 = 1 / 2, 1 / 2, 1
    w1 = F_n
    w2 = F(u + h * a21 * w1, lambda_)
    w3 = F(u + h * a32 * w2, lambda_)
    w4 = F(u + h * a43 * w3, lambda_)
    return u + h * (b1 * w1 + b2 * w2 + b3 * w3 + b4 * w4)


def erk(u, h, F_n, lambda_, p):
    if p == 1:
        return erk1(u, h, F_n, lambda_)
    elif p == 2:
        return erk2(u, h, F_n, lambda_)
    elif p == 3:
        return erk3(u, h, F_n, lambda_)
    elif p == 4:
        return erk4(u, h, F_n, lambda_)
    else:
        raise ValueError('Invalid order of approximation p = {}!'.format(p))


# def write_to_file(case):
#     print('write to file')
#     filename = str(lambda_) + '.txt'
#     with open(filename, "w") as file:
#         print('{0} {1} {2} {3}'.format('N'.rjust(10), 'err'.rjust(10), 'crit'.rjust(10), 'rich'.rjust(10)), file=file)
#         print('{0:10.4f} {1:10.4f}'.format(data['Grids'][0], data['L error'][0]), file=file)
#         for i in range(len(data['Richardson'])):
#             print('{0:10.4f} {1:10.4f} {2:10.4f} {3:10.4f}'.format(
#                 data['Grids'][i + 1], data['L error'][i + 1],
#                 data['crit2'][i], data['Richardson'][i]), file=file)


def scalar(x, y):
    # x, y of type np.ndarray
    assert(x.size == y.size)
    return (x * y).sum()


def norm(vec: np.ndarray):
    # vec - np.ndarray
    return (scalar(vec, vec))**0.5


# правая часть ОДУ
def f(u, lambda_):
    # return np.array([1., np.tan(lambda_ * u[1])])
    return np.array([1., np.sinh(lambda_ * u[1])])


# правая часть после перехода к длине дуги
def F(u, lambda_):
    # f_n = f(u)
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            res = np.array([1 / np.cosh(lambda_ * u[1]),
                            np.tanh(lambda_ * u[1])])
        except Warning:
            res = np.array([0., np.inf])
    return res


def UL(L, lambda_, u0):
    '''
    Exact solution u(l). Must support vectorizing.

    Parameters
    ----------
    L : float or np.array
        Value or vector of values of arc length.

    Returns
    -------
    float or np.array
        Value or vector of values of solution at L.
    '''
    A = np.exp(lambda_ * L) * np.sinh(lambda_ * u0[1])
    return 1 / lambda_ * (np.log(A) + np.log(1 + np.sqrt(1 + A**(-2))))


def UT(T, lambda_, u0):
    '''
    Exact solution u(t). Must support vectorizing.

    Parameters
    ----------
    T : float or np.array
        Value or vector of values of t.

    Returns
    -------
    float or np.array
        Value or vector of values of solution at t.
    '''
    # return (1 / lambda_) * np.arcsin(np.exp(lambda_ * t) * \
    #    np.sin(lambda_ * u0[1]))
    B = np.exp(lambda_ * T) * np.tanh(lambda_ * u0[1] / 2)
    return 1 / lambda_ * np.log((1 + B) / (1 - B))


def TL(L, lambda_, u0):
    '''
    Exact solution t(l). Must support vectorizing.

    Parameters
    ----------
    L : float or np.array
        Value or vector of values of arc length.

    Returns
    -------
    float or np.array
        Value or vector of values of t at L.
    '''
    A = np.exp(lambda_ * L) * np.sinh(lambda_ * u0[1])
    return 1 / lambda_ * np.log(
        np.tanh(np.log(A) + np.log(1 + np.sqrt(A**(-2) + 1)) / 2) / np.tanh(
        lambda_ * u0[1] / 2))


def richardson(u1, u2, H, p):
    # H must be array of steps from previous grid.
    u2 = u2[::2]
    min_size = min(u1.size, u2.size)
    u1 = u1[:min_size]
    u2 = u2[:min_size]
    H = H[:min_size]
    
    R = (u1 - u2) / (2**p - 1)
    R = np.sqrt((R**2 * H).sum() / H.sum())
    return R


def run_iterations(case):
    # Retrieve settings.
    lambda_ = case.lambda_
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
    F_last = F(u, lambda_)
    kappa = 1.  # on the first step.
    kappas = []
    H = []  # steps
    L = [0]
    
    while True:

        # try:
        # Extrapolate kappa from previous step.
        h = (Nmin / L_old + Nmax * kappa**0.4 / int_old)**(-1)
        H.append(h)
        L.append(L[-1] + h)
            # print('h =', h)
        # except ZeroDivisionError as e:
        #     print(e)
        #     return L, L_mas, H, U, integral, DL, kappas

        u_new = erk(u, h, F_last, lambda_, p)
        U.append(u_new)

        F_new = F(u_new, lambda_)  # rhs on next step
        kappa = norm((F_new - F_last) / h)  # real kappa on current step
        kappas.append(kappa)
        F_last = F_new

        u = u_new
        
        if u_new[1] >= 1 / lambda_ * np.arcsinh(
                (lambda_ + (lambda_**2 - 4)**0.5) / 2):
            break

    H = np.array(H)
    kappas = np.array(kappas)
    U = np.array(U)
    # Do not write the first node L = 0 because it does not matter for DL.
    L = np.array(L)

    integral = (kappa**0.4 * H).sum()
    DL = np.sqrt((((U[1:, 1] - UL(L[1:], lambda_, u0))**2 +
                   (U[1:, 0] - TL(L[1:], lambda_, u0))**2) /
                  (UL(L[1:], lambda_, u0)**2 + 
                   TL(L[1:], lambda_, u0)**2) * H).sum() / H.sum())

    case.solutions.append(U)
    case.lengths.append(L)
    case.steps.append(H)
    case.integrals.append(integral)
    case.errors.append(DL)
    
    case.Nmin.append(2 * Nmin)
    case.Nmax.append(2 * Nmax)
    
    case.grid_sizes.append(H.size)


def stage1(case):

    counter = 0
    dist = case.crit1 + 1
    while dist > case.crit1:
        '''Apply Euler scheme and double Nmin and Nmax
        until the mesh becomes adaptive.'''
        counter += 1
        print('{0}-й прогон Эйлера'.format(counter))

        run_iterations(case)

        print('N =', case.steps[-1].size)

        if counter >= 2:
            dist = 0
            H_old = case.steps[-2]
            H = case.steps[-1]
            N = min(len(H_old), len(H[::2]))
            for n in range(N):
                ksi = (H[2 * n] + H[2 * n + 1]) / H_old[n]
                dist = dist + (ksi**0.5 - ksi**(-0.5))**2 * H[2 * n]
            dist /= H.sum()
            dist **= 0.5

            case.rich_errors.append(richardson(
                case.solutions[-2][:, 1],
                case.solutions[-1][:, 1],
                H_old,
                case.approx_order_1))


def stage2(case):
    
    if case.grid_sizes[-1] >= case.num2:
        raise ValueError('Stage 1 ended with too many steps!')
    
    while case.grid_sizes[-1] < case.num2:
        lambda_ = case.lambda_
        u0 = case.u0
        H = case.steps[-1]
        N = H.size
        H_new = np.zeros(2 * N)
        
        if N == 1:
            H_new[0] = H[0] / 2
            H_new[1] = H[0] / 2
        else:
            H_new[0] = H[0] * H[0]**0.5 / (H[0]**0.5 + H[1]**0.5)
            H_new[1] = H[0] * H[1]**0.5 / (H[0]**0.5 + H[1]**0.5)
            H_new[-2] = H[-1] * H[-2]**0.5 / (H[-2]**0.5 + H[-1]**0.5)
            H_new[-1] = H[-1] * H[-1]**0.5 / (H[-2]**0.5 + H[-1]**0.5)

            if N > 2:        
                for i in range(1, N - 1):
                    H_new[2 * i] = H[i] * H[i - 1]**0.25 / \
                        (H[i - 1]**0.25 + H[i + 1]**0.25)
                    H_new[2 * i + 1] = H[i] * H[i + 1]**0.25 / \
                        (H[i - 1]**0.25 + H[i + 1]**0.25)

        H_old = H
        H = H_new

        u = case.u0
        U = [u]
        L = [0]
        for h in H:
            L.append(L[-1] + h)
            u = erk(u, h, F(u), lambda_, case.approx_order_2)
            U.append(u)
            
        L = np.array(L)
        U = np.array(U)
        DL = np.sqrt((((U[1:, 1] - UL(L[1:], lambda_, u0))**2 +
                       (U[1:, 0] - TL(L[1:], lambda_, u0))**2) /
                      (UL(L[1:], lambda_, u0)**2 + 
                       TL(L[1:], lambda_, u0)**2) * H).sum() / H.sum())
            
        case.grid_sizes.append(H.size)
        case.solutions.append(U)
        case.lengths.append(L)
        case.steps.append(H)
        case.errors.append(DL)
        if len(case.grid_sizes) > 1:
            case.rich_errors.append(richardson(
                case.solutions[-2][:, 1],
                case.solutions[-1][:, 1],
                H_old,
                case.approx_order_2))


def run(cases):
    
    for case in cases:
        
        print('Расчет задачи с lambda =', case.lambda_)
        print('Построение адаптивной сетки\n')
        
        stage1(case)
        
        case.switch = case.grid_sizes.size
        msg = 'Переход с первого этапа на второй происходит на'
        msg += case.switch
        msg += 'сетке'
        print(msg)

        stage2(case)
    
    return cases
                    
