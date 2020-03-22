# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 09:05:12 2020

@author: Matthew
"""

import scipy.io as spio
import numpy as np
from Cit_par import *
import matplotlib.pyplot as plt
import control.matlab as control
from CG_per_time import instantanious_weight

matlab = spio.loadmat('matlab.mat')
time = matlab['flightdata'][0][0][47][0][0][0].transpose()
data = time   
for i in range(47):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis = 1)
    
#======================================================================================

# EOM in state space form
"""
The A matrix is set up as was explained in the flight dynamics project manual.
This means that the non dimensional components in the state vector have been made dimensional.
Thus: u_hat -> u    multiply the first column of matrix A by  1/V      (symmetric)
      qc/V  -> q    multiply the last column of matrix A by   c/V     (symmetric)
      pb/2V -> p    multiply the third column of matrix A by  b/2V   (asymmetric)
      rb/2V -> r    multiply the last column of matrix A by   b/2V    (asymmetric)
"""

# Assumption
#CXq = 0

# Symmetric
C1 = np.matrix([[-2 * muc * c / V0, 0, 0, 0],
                [0, (CZadot - 2 * muc) * c / V0, 0, 0],
                [0, 0, -c / V0, 1],
                [0, Cmadot * c / V0, 0, -2 * muc * KY2 * c / V0]])
C1[:, 0] = (1/V0) * C1[:,0]
C1[:, 3] = (c/V0) * C1[:, 3]

C1_inv = np.linalg.inv(C1)

C2 = np.matrix([[CXu, CXa, CZ0, CXq],
                [CZu, CZa, -CX0, CZq + 2 * muc],
                [0, 0, 0, 1],
                [Cmu, Cma, 0, Cmq]])
C2[:, 0] = (1/V0) * C2[:,0]
C2[:, 3] = (c/V0) * C2[:, 3]

C3 = np.matrix([[-CXde],
                [-CZde],
                [0],
                [-Cmde]])

A_sym = - C1_inv * C2

B_sym = C1_inv * C3

C_sym = np.identity(4)

D_sym  = np.matrix([[0],
                    [0],
                    [0],
                    [0]])

symmetric = control.ss(A_sym,B_sym,C_sym,D_sym)     # Set up a system
control.damp(symmetric)
control.pzmap(symmetric)

# Asymmetric
B1 = np.matrix([[(CYbdot - 2 * mub) * (c/V0), 0, 0, 0],
                [0, -c/(2*V0), 0, 0],
                [0, 0, -4 * mub * KX2 * (c/V0), 4 * mub * KXZ * (c/V0)],
                [Cnbdot * (c/V0), 0, 4 * mub * KXZ * (c/V0), -4 * mub * KZ2 * (c/V0)]])
B1[:, 2] = (b/(2*V0)) * B1[:,2]
B1[:, 3] = (b/(2*V0)) * B1[:, 3]

B1_inv = np.linalg.inv(B1)

B2 = np.matrix([[CYb, CL, CYp, CYr - 4 * mub],
                [0, 0, 1, 0],
                [Clb, 0, Clp, Clr],
                [Cnb, 0, Cnp, Cnr]])
B2[:, 2] = (b/(2*V0)) * B2[:,2]
B2[:, 3] = (b/(2*V0)) * B2[:, 3]

B3 = np.matrix([[-CYda, -CYdr],
                [0, 0],
                [-Clda, -Cldr],
                [-Cnda, -Cndr]])

A_asym = - B1_inv * B2

B_asym = B1_inv * B3

C_asym = np.identity(4)

D_asym = np.matrix([[0, 0],
                    [0, 0],
                    [0, 0],
                    [0, 0]])

asymmetric = control.ss(A_asym, B_asym, C_asym, D_asym)

# control.damp(asymmetric)
# control.pzmap(asymmetric)

X0 = np.array([0.01, 5.43746, 9.97352 , -0.10809])


# t = np.arange(0, 15.1, 0.1)
# y_out, t_out = control.impulse(asymmetric, t, input = 1)
# plt.plot(t_out, y_out[:,3])
# plt.show()
