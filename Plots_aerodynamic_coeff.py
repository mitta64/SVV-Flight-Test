# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 15:28:43 2020

@author: Matthew
"""
import scipy.io as spio
import numpy as np
from Parameters_at_altitude import V_TAS, Rho
from CG_per_time import ramp_mass
from Cit_par import S

# ======================================================================================
    
def weight_array(ramp_mass):
    fuel_used = 0.453592 * np.matrix([475, 524, 553, 577, 598, 640])
    weight = np.full((1,6), ramp_mass) - fuel_used
        
    
    return weight

def lift_coeff(weight, V_C, altitude, temp_measured_total, wingspan):
    CL = np.empty(6)
    

    for i in range(np.size(weight)):
        CL[i] = weight[0,i] / (0.5 * Rho(altitude[0,i]) * (V_TAS(V_C[0,i], altitude[0,i], temp_measured_total[0,i]))**2 *S)
        
    
    return CL


# ======================================================================================
    
# weight changes over time
weight = weight_array(ramp_mass)

# Lift coefficient 
V_C = 0.51444 * np.matrix([247, 214, 192, 154, 129, 118])
altitude = 0.3048 * np.matrix([11990, 11980, 11990, 11990, 11990, 11990])
temp_measured_total = np.add(273, np.matrix([-2.8, -5.5, -7.2, -9.8, -12, -12]))
CL = lift_coeff(weight, V_C, altitude, temp_measured_total, S)
# ======================================================================================



