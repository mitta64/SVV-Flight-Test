# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 09:05:12 2020

@author: Matthew
"""

import scipy.io as spio
import numpy as np
from Cit_par import *
import control.matlab as control

#=====================================================================================


matlab = spio.loadmat('FTISxprt-20200306_flight3.mat')
time = matlab['flightdata'][0][0][48][0][0][0].transpose()
data = time   
for i in range(48):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis = 1)
    

#======================================================================================


def instantanious_weight(t,bem,block,data, weight_seats): # contains hardcoded passenger weights

    fuel_use = data[:, 14] + data[:, 15]  # first  is left engine 2nd is right engine fuel use

    fuel_present = block - fuel_use[np.where(t == data[:,0])[0]]

    current_total = bem + block + np.sum(weight_seats[:,0])  - fuel_use[np.where(t == data[:,0])[0]]

    zero_fuel_mass = np.sum(weight_seats[:,0])  + bem
    ramp_mass = bem + block + np.sum(weight_seats[:,0])

    return current_total, fuel_present, zero_fuel_mass, ramp_mass

def fuel_moment(t,bem,block,weight_seats):
    # load the moment data

    path = 'moment_data.csv'
    file = open(path, "r")
    moment_weight = np.genfromtxt(path, delimiter=",", skip_header=1)
    file.close()

    moment_weight[:,1] = 100 * moment_weight[:,1]

    current_weight, fuel_present , zero_fuel,ramp_mass = instantanious_weight(t, bem, block, data, weight_seats)

    fuel_moment = 0

    for i in range(np.shape(moment_weight)[0]):
        if i == (np.shape(moment_weight)[0] - 1):
            a = 0

        else:
            if fuel_present >= moment_weight[i, 0] and fuel_present <= moment_weight[i + 1, 0]:
                x1 = moment_weight[i, 0]
                x2 = moment_weight[i + 1, 0]
                y1 = moment_weight[i, 1]
                y2 = moment_weight[i + 1, 1]

                gradient = (y2 - y1) / (x2 - x1)
                y_comp = y1 - gradient * x1

                fuel_moment = gradient * fuel_present + y_comp


    return fuel_moment

def moment_equilibrium(weight_seats,current_weight,zero_fuel,ramp_mass,fuel_moment_datum,bem,arm_bem,mac):
    m_seats = np.sum(weight_seats[:, 0] * weight_seats[:,1])      # mass in lbs times distance to datum in inch
    m_bem = bem*arm_bem
    m_block = 11705.50 *100 # From the table read of at 4100 lbs

    total_moment_datum = m_seats + m_bem + fuel_moment_datum
    xcg_datum = total_moment_datum / current_weight

    # standard locations
    xcg_datum_bem = m_bem/bem
    xcg_datum_0fuel = (m_seats + m_bem) / zero_fuel
    xcg_datum_ramp = (m_seats + m_bem + m_block)/ramp_mass

    # cg normalized wrt MAC
    #xcg = (xcg_datum - 261.45) * 0.0254 #xcg in distance from start MAC
    xcg = (xcg_datum - 261.45)*0.0254*100/(mac*0.0254) #xcg is in %c MAC


    return xcg_datum_bem,xcg_datum_0fuel,float(xcg_datum_ramp),float(xcg_datum),float(xcg)

# ========================================= MAIN =======================================================================


block = 4100          # in lbs
bem = 9165            # in lbs
arm_bem = 291.65      # in inch
mac = 80.98           # in inch
t = 9                 # in sec



# seating defined
weight_seats = np.array([[80, 131], [102, 131], [66, 214], [81, 214], [99, 251], [64, 251], [81, 288], [99, 288], [88, 170]])
weight_seats[:, 0] = 2.20462262 * weight_seats[:, 0]  # convert weight in kg to lbs


# weight determination
current_weight,fuel_present, zero_fuel_mass,ramp_mass = instantanious_weight(t, bem, block, data, weight_seats)

# moment determined wrt datum caused by fuel consumption
fuel_moment_datum = fuel_moment(t,bem,block,weight_seats)


# cg contains: BEM cg wrt datum, Zero fuel cg wrt datum, Ramp cg wrt datum, current cg wrt datum, cg wrt MAC
cg = moment_equilibrium(weight_seats,current_weight,zero_fuel_mass,ramp_mass,fuel_moment_datum,bem,arm_bem,mac)



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
A_sym[:,0] =  V0 * A_sym[:,0]            #multiply first column with V
A_sym[:,3] = (V0/c) * A_sym[:,3]

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


# Asymmetric
A_asym = np.matrix([[(V0 * CYb)/(b * 2 * mub), (V0 * CL)/(b * 2 * mub), (V0 * CYp)/(b * 2 * mub), (V0 * (CYr - 4 * mub))/(b * 2 * mub)],
                    [0, 0, 2 * V0 / b, 0],
                    [(V0 * (Clb * KZ2 + Cnb * KXZ)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), 0, (V0 * (Clp * KZ2 + Cnp * KXZ)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), (V0 * (Clr * KZ2 + Cnr * KXZ)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2))],
                    [(V0 * (Clb * KXZ + Cnb * KX2)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), 0, (V0 * (Clp * KXZ + Cnp * KX2)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2)), (V0 * (Clr * KXZ + Cnr * KX2)) / (b * 4 * mub * (KX2 * KZ2 - KXZ**2))]])
A_asym[:, 2] = (2 * V0 / b) * A_asym[:,2]
A_asym[:, 3] = (2 * V0 / b) * A_asym[:,3]

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

