# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 16:32:52 2019

@author: Artem
"""

import matplotlib.pyplot as plt
import numpy as np


def erk1(u, h, F_n):
    return u + h*F_n


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


def write_to_file(data):
    print('write to file')
    filename = str(lambda_) + '.txt'
    with open(filename, "w") as file:
        print('{0} {1} {2} {3}'.format('N'.rjust(10), 'err'.rjust(10), 'crit'.rjust(10), 'rich'.rjust(10)), file=file)
        print('{0:10.4f} {1:10.4f}'.format(data['Grids'][0], data['L error'][0]), file=file)
        for i in range(len(data['Richardson'])):
            print('{0:10.4f} {1:10.4f} {2:10.4f} {3:10.4f}'.format(
                data['Grids'][i + 1], data['L error'][i + 1],
                data['crit2'][i], data['Richardson'][i]), file=file)


# u - вектор-функция, u[0] - функция-время после автономизации,
# u[1:] - искомая вектор-функция
a = 0.0  # левый конец отрезка по времени
# Nmax - число шагов по всему отрезку,
# Nmin - число шагов на регулярных участках,
# kappa - кривизна на предыдущем шаге,
# считаем шаг h с учетом кривизны решения
# берем длину дуги, кривизну на первом шаге и интеграл с потолка
Lmax = 1
Nmin = 5
Nmax = 10


# правая часть ОДУ
def f(u):
    # return np.array([1., u[1] * (u[0] - u[1]) / 0.5])
    # return np.array([1., np.tan(lambda_*u[1])])
    # return np.array([
    #        1., -lambda_ * np.cos(u[0]) * (u[1]**2 - a**2)**2 / (u[1]**2 + a**2)
    #        ])
    # return np.array([1., 2*u[0]])
    return np.array([1., np.sinh(lambda_ * u[1])])


def scalar(x, y):
    # x, y of type np.ndarray
    assert(x.size == y.size)
    return (x * y).sum()


def norm(vector):
    # vector - np.ndarray
    return (scalar(vector, vector)) ** 0.5


# правая часть после перехода к длине дуги
def F(u):
    f_n = f(u)
    return f_n / norm(f_n)


# точное решение
def UL(L):
    # return (1/lambda_) * 2 * np.atan(np.exp(lambda_*L) * np.tan(lambda_*u0[1]/2))
    A = np.exp(lambda_ * L) * np.sinh(lambda_ * u0[1])
    return 1 / lambda_ * (np.log(A) + np.log(1 + np.sqrt(1 + A ** (-2))))


def UT(t):
    return (1 / lambda_) * np.arcsin(np.exp(lambda_ * t) * np.sin(lambda_ * u0[1]))
#    B = np.exp(lambda_*t) * np.tanh(lambda_ * u0[1] / 2)
#    return 1/lambda_ * np.log((1 + B) / (1 - B))


def TL(L):
    A = np.exp(lambda_ * L) * np.sinh(lambda_ * u0[1])
    return 1 / lambda_ * np.log(np.tanh(np.log(A + np.sqrt(A**2 + 1)) / 2) / np.tanh(lambda_ * u0[1] / 2))


print('Построение адаптивной сетки\n')


def run_iterations(L, L_mas, L_mas_old, H, U, integral, old_integral, Nmin, Nmax, solver):
    u = u0
    F_n1 = F(u)

    # с потолка
    kappa = 1
    kappas = [kappa]

    U.append([])
    U[-1].append(u)

    DL = 0

    while True:

        try:
            h = (Nmin / L_mas_old[-1] + Nmax * kappa ** (2 / 5) * old_integral ** (-1)) ** (-1)
        except ZeroDivisionError as e:
            print(e)
            return L, L_mas, H, U, integral, DL, kappas

        F_n = F_n1
        u_new = solver['scheme'](u, h, F_n)

        #        if u_new[0] > T:
        #
        #            u_new[1] = u[1] + F_n[1] * (T - u[0])
        #            u_new[0] = T
        #            h = ((T - u[0])**2 + (u_new[1] - u[1])**2) ** 0.5
        #
        #            L += h
        #            L_mas.append(L)
        #            H.append(h)
        #
        #            #L2
        #            DL += (u_new[1] - UL(L)) ** 2 * h
        #            D += (u_new[1] - UL(L)) ** 2
        #            #DT += (u_new[1] - UT(u_new[0])) ** 2 * h
        #
        #            U[-1].append(u_new)
        #            integral += kappa ** (2/5) * h
        #            break

        L += h

        H.append(h)
        L_mas.append(L)
        U[-1].append(u_new)

        # L2
        try:
            DL += ((u_new[1] - UL(L))**2 + (u_new[0] - TL(L))**2) / (UL(L)**2 + TL(L)**2) * h
        except OverflowError as e:
            print(e)
            return L, L_mas, H, U, integral, DL, kappas
        # DT += (u_new[1] - UT(u_new[0])) ** 2 * h

        F_n1 = F(u_new)
        kappa = norm((F_n1 - F_n) / h)
        kappas.append(kappa)

        #        w3 = F(u)
        #        kappa = norm((2*w2 - 2*w3) / h)

        integral += kappa ** (2 / 5) * h
        u = u_new

        if u_new[1] >= 1 / lambda_ * np.arcsinh((lambda_ + (lambda_ ** 2 - 4) ** 0.5) / 2):
            break

    DL /= sum(H)
    DL **= 0.5

    return L, L_mas, H, U, integral, DL, kappas


def stage1(lambda_, solver, break_crit):
    Nmin = 6
    Nmax = 20
    # настоящий интеграл (посчитаем его после первого прогона)
    integral = 1

    U = []

    # текущая длина дуги
    L_mas = [1]
    # массив шагов
    H = []
    # массив чисел шагов
    N_mas = []
    # массив интегралов
    integral_mas = []

    crit_mas = [break_crit * 10]
    DL_mas = []
    DT_mas = []
    R_mas = []

    print('\nlambda =', lambda_)
    counter = 0

    # Эйлерим и удваиваем Nmax и Nmin, пока сетка не станет адаптивной
    while crit_mas[-1] > break_crit:
        counter += 1
        print('{0}-й прогон Эйлера'.format(counter))

        L = 0

        L_mas_old = L_mas
        L_mas = [0]
        H_old = H
        H = []

        old_integral = integral
        integral = 0

        L, L_mas, H, U, integral, DL, kappas = run_iterations(L, L_mas, L_mas_old, H, U, integral, old_integral, Nmin,
                                                              Nmax, solver)
        Nmin *= 2
        Nmax *= 2

        if counter >= 2:
            crit_mas.append(0)
            N = min(len(H_old), len(H) // 2)
            for n in range(N):
                ksi = (H[2 * n] + H[2 * n + 1]) / H_old[n]
                crit_mas[-1] += (ksi ** 0.5 - ksi ** (-0.5)) ** 2 * H[2 * n]
            crit_mas[-1] /= L
            crit_mas[-1] **= 0.5

        N_mas.append(len(H))

        integral_mas.append(integral)

        if counter >= 2:
            R_mas.append(richardson(U, H, solver['p']))

        DL_mas.append(DL)

    return U, H, DL_mas, DT_mas, N_mas, integral_mas, L_mas, crit_mas[1:], R_mas, kappas


def stage2(lambda_, u0, T, U, DL_mas, N_mas, R_mas, solver, H):

    while np.log10(N_mas[-1]) < 4:

        H_new = [0 for i in range(2 * len(H))]
        for i in range(1, len(H) - 1):
            H_new[2 * i] = H[i] * H[i - 1] ** 0.25 / (H[i - 1] ** 0.25 + H[i + 1] ** 0.25)
            H_new[2 * i + 1] = H[i] * H[i + 1] ** 0.25 / (H[i - 1] ** 0.25 + H[i + 1] ** 0.25)
        H_new[0] = H[0] * H[0] ** 0.5 / (H[0] ** 0.5 + H[1] ** 0.5)
        H_new[1] = H[0] * H[1] ** 0.5 / (H[0] ** 0.5 + H[1] ** 0.5)
        H_new[-2] = H[-1] * H[-2] ** 0.5 / (H[-2] ** 0.5 + H[-1] ** 0.5)
        H_new[-1] = H[-1] * H[-1] ** 0.5 / (H[-2] ** 0.5 + H[-1] ** 0.5)

        H = H_new

        U.append([u0])
        DL_mas.append(0)
        u = u0
        L = 0
        for h in H:
            L = L + h
            u = solver['scheme'](u, h, F(u))
            U[-1].append(u)
            DL_mas[-1] = DL_mas[-1] + ((u[1] - UL(L))**2 + (u[0] - TL(L))**2) / (UL(L)**2 + TL(L)**2) * h
        DL_mas[-1] /= sum(H)
        DL_mas[-1] **= 0.5
        N_mas.append(len(H))

        R_mas.append(richardson(U, H, solver['p']))

    return U, DL_mas, N_mas, R_mas


def richardson(U, H, p):
    R = [(U[-2][j][1] - U[-1][2 * j][1]) / (2 ** p - 1) for j in range(min(len(U[-2]), len(U[-1]) // 2))]
    R = [R[i] ** 2 * H[2 * i] for i in range(len(R))]
    R = (sum(R)) ** 0.5
    return R


Lambdas = [10, 100, 1000, 10000, 1e5]
break_crit = 1e-1

# lambda_ = Lambdas[-1]
lambda_ = 100
u0 = np.array([0, 0.001])
T_pole = (-1 / lambda_) * np.log(np.sin(lambda_ * u0[1]))  # полюс производной
T = T_pole - 0.01 * T_pole
L_extr = 1 / lambda_ * np.log(1 / np.sinh(lambda_ * u0[1]))
L_MAX = 2 * L_extr

L_array = []

data = {}

solver_stage1 = {'scheme': erk2, 'p': 2}
solver_stage2 = {'scheme': erk4, 'p': 4}

errors = plt.figure()
for lambda_ in Lambdas:
    u0 = np.array([0, 1 / lambda_ * np.arcsinh((lambda_ - (lambda_ ** 2 - 4) ** 0.5) / 2)])
    T_pole = (-1 / lambda_) * np.log(np.sin(lambda_ * u0[1]))  # полюс производной
    print(T_pole)
    U, H, DL_mas, DT_mas, N_mas, integral_mas, L_mas, crit_mas, R_mas, kappas = stage1(lambda_, solver_stage1,
                                                                                       break_crit)
    switch = len(N_mas) - 1
    print('Переход с первого этапа на второй происходит на', switch + 1, 'сетке')

    U, DL_mas, N_mas, R_mas = stage2(lambda_, u0, T, U, DL_mas, N_mas, R_mas, solver_stage2, H)
    L_array.append(L_mas[-1])
    L_extr = 1 / lambda_ * np.log(1 / np.sinh(lambda_ * u0[1]))
    data[lambda_] = {
        'u0': u0,
        'L_extr': L_extr,
        'Solutions': U,
        'L error': np.log10(DL_mas),
        'T error': np.log10(DT_mas),
        'Grids': np.log10(N_mas),
        'L array': L_mas,
        'switch': switch,
        'crit2': np.log10(crit_mas),
        'Richardson': np.log10(R_mas),
        'Kappa': kappas}

    line1 = np.array([-j for j in np.log10(N_mas)])  # for erk1 on stage1
    line_p = np.array([-solver_stage2['p']*j for j in np.log10(N_mas)])  # for stage2
    fig = plt.figure()

    # # crit
    # crit, = plt.plot(np.log10(N_mas[1:]), np.log10(crit_mas), 'ko-', label='crit')

    line1, = plt.plot(np.log10(N_mas), line1, 'co-', label='45 degree line')
    line_p, = plt.plot(np.log10(N_mas), line_p, 'co-', label='tg(alpha) = -p')

    err, = plt.plot(np.log10(N_mas), np.log10(DL_mas), 'yo-', label='err')
    # plt.plot(np.log10(N_mas), np.log10(DT_mas), 'mo-')
    # plt.plot(np.log10(N_mas), np.log10(D_mas), 'go-')

    rich, = plt.plot(np.log10(N_mas[1:]), np.log10(R_mas), 'r^-', label='richardson')

    plt.title('Lambda = ' + str(lambda_))
    plt.xlabel('lg(N)')
    plt.ylabel('lg(err)')
    plt.legend(handles=[line1, line_p, err, rich])
    plt.show()
    fig.savefig('err' + str(lambda_) + '.png')

    X = [U[-1][i][0] for i in range(len(U[-1]))]
    Y = [U[-1][i][1] for i in range(len(U[-1]))]
    fig1 = plt.figure()
    plt.plot(X, Y)
    plt.title('Lambda = ' + str(lambda_))
    plt.xlabel('t')
    plt.ylabel('u(t)')
    plt.show()
    fig1.savefig('graph' + str(lambda_) + '.png')

plt.figure()

line10, = plt.plot(data[10]['Grids'], data[10]['L error'], '-ob', linewidth=5, markersize=10, label='lambda = 10')
switch = data[10]['switch']
yswitch = data[10]['L error'][switch]
plt.vlines(data[10]['Grids'][switch], yswitch - .7, yswitch + .7, colors='b')

line100, = plt.plot(data[100]['Grids'], data[100]['L error'], '-og', linewidth=5, markersize=10, label='lambda = 100')
switch = data[100]['switch']
yswitch = data[100]['L error'][switch]
plt.vlines(data[100]['Grids'][switch], yswitch - .7, yswitch + .7, colors='g')

line1000, = plt.plot(data[1000]['Grids'], data[1000]['L error'], '-or', linewidth=5, markersize=10, label='lambda = 1000')
switch = data[1000]['switch']
yswitch = data[1000]['L error'][switch]
plt.vlines(data[1000]['Grids'][switch], yswitch - .7, yswitch + .7, colors='r')

line10000, = plt.plot(data[10000]['Grids'], data[10000]['L error'], '-om', linewidth=5, markersize=10, label='lambda = 10000')
switch = data[10000]['switch']
yswitch = data[10000]['L error'][switch]
plt.vlines(data[10000]['Grids'][switch], yswitch - .7, yswitch + .7, colors='m')

line100000, = plt.plot(data[100000]['Grids'], data[100000]['L error'], '-oc', linewidth=5, markersize=10, label='lambda = 100000')
switch = data[100000]['switch']
yswitch = data[100000]['L error'][switch]
plt.vlines(data[100000]['Grids'][switch], yswitch - .7, yswitch + .7, colors='c')

plt.legend(handles=[line10, line100, line1000, line10000, line100000])
plt.xlabel('lg(N)')
plt.ylabel('lg(error)')
plt.title('Calculations with ERK2 -> ERK4')
plt.savefig('D:/Docs/MSU/NumericalMethods/Programs/GEAM/20092020/ERK1.png')
