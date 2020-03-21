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
Thus: u_hat -> u    multiply the first column of matrix A by V      (symmetric)
      qc/V  -> q    multiply the last column of matrix A by V/c     (symmetric)
      pb/2V -> p    multiply the third column of matrix A by 2V/b   (asymmetric)
      rb/2V -> r    multiply the last column of matrix A by 2V/b    (asymmetric)
"""

"""
The following needs to be checked:
    CL in the reader is a weight component but in Cit_car the lift coefficient
    How are the changes in V0 taken into account
"""
# Assumption
CXq = 0

# Symmetric
A_sym = np.matrix([[(V0 * CXu) / (c * 2 * muc ), (V0 * CXa) / (c * 2 * muc ), (V0 * CZ0) / (c * 2 * muc ), (V0 * CXq) / (c * 2 * muc )],
                    [(V0 * CZu) / (c * (2 * muc - CZadot)), (V0 * CZa) / (c * (2 * muc - CZadot)), (-V0 * CX0) / (c * (2 * muc - CZadot)), (V0 * (2 * muc + CZq)) / (c * (2 * muc - CZadot))],
                    [0, 0, 0, V0/c],
                    [V0 * ((Cmu + CZu * (Cmadot / (2 * muc - CZadot))) / (c * 2 * muc * KY2)), V0 * ((Cma + CZa * (Cmadot / (2 * muc - CZadot))) / (c * 2 * muc * KY2)), -V0 * ((CX0 * (Cmadot / (2 * muc - CZadot))) / (c * 2 * muc * KY2)), V0 * ((Cmq + Cmadot * ((2 * muc + CZq) / (2 * muc - CZadot))) / (c * 2 * muc * KY2))]])
#A_sym[:,0] =  (1/V0) * A_sym[:,0]            #multiply first column with V
#A_sym[:,3] = (c/V0) * A_sym[:,3]

B_sym = np.matrix([[(V0 * CXde) / (c * 2 * muc)], 
                    [(V0 * CZde) / (c *( 2 * muc - CZadot))],
                    [0],
                    [V0 * ((Cmde + CZde * (Cmadot / (2 * muc - CZadot))) / (c * 2 * muc * KY2))]])

C_sym = np.identity(4)

D_sym  = np.matrix([[0],
                    [0],
                    [0],
                    [0]])
symmetric = control.ss(A_sym,B_sym,C_sym,D_sym)
eigenvalues = np.linalg.eig(A_sym)
control.damp(symmetric)
# control.pzmap(symmetric)

# Asymmetric
A_asym = np.matrix([[(V0 * CYb)/(b * 2 * mub), (V0 * CL)/(b * 2 * mub), (V0 * CYp)/(b * 2 * mub), (V0 * (CYr - 4 * mub))/(b * 2 * mub)],
                    [0, 0, 2 * V0 / b, 0],
                    [(V0 * (Clb * KZ2 + Cnb * KXZ)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), 0, (V0 * (Clp * KZ2 + Cnp * KXZ)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), (V0 * (Clr * KZ2 + Cnr * KXZ)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2))],
                    [(V0 * (Clb * KXZ + Cnb * KX2)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), 0, (V0 * (Clp * KXZ + Cnp * KX2)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), (V0 * (Clr * KXZ + Cnr * KX2)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2))]])
# A_asym[:, 2] = (b / 2 * V0) * A_asym[:,2]
# A_asym[:, 3] = (b / 2 * V0) * A_asym[:,3]

B_asym = np.matrix([[0, (V0 * (CYr - 4 * mub))/(b * 2 * mub)],
                    [0, 0],
                    [(V0 * (Clda * KZ2 + Cnda * KXZ)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), (V0 * (Cldr * KZ2 + Cndr * KXZ)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2))],
                    [(V0 * (Clda * KXZ + Cnda * KX2)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), (V0 * (Cldr * KXZ + Cndr * KX2)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2))]])

C_asym = np.identity(4)

D_asym = np.matrix([[0, 0],
                    [0, 0],
                    [0, 0],
                    [0, 0]])

asymmetric = control.ss(A_asym, B_asym, C_asym, D_asym)

# control.damp(asymmetric)
# control.pzmap(asymmetric)

X0 = np.array([V0, 5.43746, 9.97352 , -0.10809])


# t = np.arange(0, 60.1, 0.1)
# y_out, t_out = control.initial(symmetric, t, X0)
# plt.plot(t_out, y_out[:,0])
# plt.show()