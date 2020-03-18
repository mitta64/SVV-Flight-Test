# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 15:28:43 2020

@author: Matthew
"""
import matplotlib.pyplot as plt
import numpy as np
from Parameters_at_altitude import V_TAS, Rho
from CG_per_time import ramp_mass
from Cit_par import S, g, A

# ======================================================================================
    
def weight_array(ramp_mass):
    fuel_used =  np.matrix([475, 524, 553, 577, 598, 640])
    weight = 0.453592 * (np.full((1,6), ramp_mass) - fuel_used)
        
    
    return weight

def lift_coeff(weight, V_C, altitude, temp_measured_total, wingspan):
    CL = np.empty(6)
    

    for i in range(np.size(weight)):
        CL[i] = (weight[0,i] * g) / (0.5 * Rho(altitude[0,i]) * (V_TAS(V_C[0,i], altitude[0,i], temp_measured_total[0,i]))**2 *S)
           
    return CL

def drag_coeff(thrust, V_C, altitude, temp_measured_total, wingspan):
    CD = np.empty(6)
    

    for i in range(np.size(weight)):
        CD[i] = (thrust[0,i] ) / (0.5 * Rho(altitude[0,i]) * (V_TAS(V_C[0,i], altitude[0,i], temp_measured_total[0,i]))**2 *S)
            
    return CD


    


# ======================================================================================
    
# weight changes over time
weight = weight_array(ramp_mass)

# Angle of attack [deg]
alpha = np.array([1.7, 2.7, 3.5, 6.1, 8.9, 11]) 

# Thrust
thrust = np.matrix([3100.32, 2191.92, 1961.72, 1733.292, 2301.2, 2077.536])

# Calibrated airspeed [m/s]
V_C = 0.51444 * np.subtract(np.matrix([247, 214, 192, 154, 129, 118]), 2)

# Altitude [m]
altitude = 0.3048 * np.matrix([11990, 11980, 11990, 11990, 11990, 11990])

# Total temp [K]
temp_measured_total = np.add(273.15, np.matrix([-2.8, -5.5, -7.2, -9.8, -12, -12]))

# Lift coefficient 
CL = lift_coeff(weight, V_C, altitude, temp_measured_total, S)
CL = np.array(CL)
    #np.reshape(CL, (1,6))


# Drag coefficient
CD = drag_coeff(thrust, V_C, altitude, temp_measured_total, S)
CD = np.array(CD)
#np.reshape(CD, (1,6))



# Plot CL over alpha

#plt.plot(alpha, CL)
#plt.plot(CL**2, CD)


# CD0 & Oswald efficiency
plt.scatter(CL**2, CD)
z = np.polyfit((CL**2).flatten(), CD.flatten(), 1)
p = np.poly1d(z)

plt.plot(CL**2,p(CL**2),"r--")
plt.title("y=%.6fx+%.6f"%(z[0],z[1])) 

CD0 = z[1]
e = 1 / (np.pi * A * z[0])
print(CD0, e)
plt.show()

# ======================================================================================



