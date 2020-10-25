# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 01:19:55 2020

@author: kuzne
"""

from solver_data import SolverData
from GEAM import run
import numpy as np
# default solver params are set in SolverData.__init__

case = SolverData(1000)
run([case])

# line_1 = np.array([-j for j in np.log10(N_mas)])  # for erk1 on stage1
# line_p = np.array([-solver_stage2['p']*j for j in np.log10(N_mas)])  # for stage2

# line1, = plt.plot(np.log10(N_mas), line1, 'co-', label='45 degree line')
# line_p, = plt.plot(np.log10(N_mas), line_p, 'co-', label='tg(alpha) = -p')
# err, = plt.plot(np.log10(N_mas), np.log10(DL_mas), 'yo-', label='err')
# rich, = plt.plot(np.log10(N_mas[1:]), np.log10(R_mas), 'r^-', label='richardson')

# plt.title('Lambda = ' + str(lambda_))
# plt.xlabel('lg(N)')
# plt.ylabel('lg(err)')
# plt.legend(handles=[line1, line_p, err, rich])
# plt.show()
# fig.savefig('err' + str(lambda_) + '.png')

# X = [U[-1][i][0] for i in range(len(U[-1]))]
# Y = [U[-1][i][1] for i in range(len(U[-1]))]
# fig1 = plt.figure()
# plt.plot(X, Y)
# plt.title('Lambda = ' + str(lambda_))
# plt.xlabel('t')
# plt.ylabel('u(t)')
# plt.show()
# fig1.savefig('graph' + str(lambda_) + '.png')        

# line10, = plt.plot(data[10]['Grids'], data[10]['L error'], '-ob', linewidth=5, markersize=10, label='lambda = 10')
# switch = data[10]['switch']
# yswitch = data[10]['L error'][switch]
# plt.vlines(data[10]['Grids'][switch], yswitch - .7, yswitch + .7, colors='b')

# line100, = plt.plot(data[100]['Grids'], data[100]['L error'], '-og', linewidth=5, markersize=10, label='lambda = 100')
# switch = data[100]['switch']
# yswitch = data[100]['L error'][switch]
# plt.vlines(data[100]['Grids'][switch], yswitch - .7, yswitch + .7, colors='g')

# line1000, = plt.plot(data[1000]['Grids'], data[1000]['L error'], '-or', linewidth=5, markersize=10, label='lambda = 1000')
# switch = data[1000]['switch']
# yswitch = data[1000]['L error'][switch]
# plt.vlines(data[1000]['Grids'][switch], yswitch - .7, yswitch + .7, colors='r')

# line10000, = plt.plot(data[10000]['Grids'], data[10000]['L error'], '-om', linewidth=5, markersize=10, label='lambda = 10000')
# switch = data[10000]['switch']
# yswitch = data[10000]['L error'][switch]
# plt.vlines(data[10000]['Grids'][switch], yswitch - .7, yswitch + .7, colors='m')

# line100000, = plt.plot(data[100000]['Grids'], data[100000]['L error'], '-oc', linewidth=5, markersize=10, label='lambda = 100000')
# switch = data[100000]['switch']
# yswitch = data[100000]['L error'][switch]
# plt.vlines(data[100000]['Grids'][switch], yswitch - .7, yswitch + .7, colors='c')

# line1000000, = plt.plot(data[1000000]['Grids'], data[1000000]['L error'], '-ok', linewidth=5, markersize=10, label='lambda = 1000000')
# switch = data[1000000]['switch']
# yswitch = data[1000000]['L error'][switch]
# plt.vlines(data[1000000]['Grids'][switch], yswitch - .7, yswitch + .7, colors='k')

# plt.legend(handles=[line10, line100, line1000, line10000, line100000, line1000000])
# plt.xlabel('lg(N)')
# plt.ylabel('lg(error)')
# plt.title('Calculations with ERK2 -> ERK4')
# plt.savefig('D:/Docs/MSU/NumericalMethods/Programs/GEAM/20092020/ERK1.png')
