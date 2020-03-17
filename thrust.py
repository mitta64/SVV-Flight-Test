# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 13:45:58 2020

@author: HJ Hoogendoorn

File for Thrust calculations
"""

from Matlab_converter import data
import numpy as np

pressure_alt    = data[:,37]
Temp            = data[:,36]
fuelflowleft    = data[:,4]
fuelflowright   = data[:,5]
mach            = data[:,40]

#deltaT calculations
delta_T = []
for i in range(0,len(Temp)):
    T_ISA = 15.15 -0.0065*(pressure_alt[i]*0.3048)
    delta_T.append(Temp[i]-T_ISA)
    


delta_T = np.array(delta_T)




#creating .dat file for Thrust.exe
np.set_printoptions(suppress=True)
DataOut=np.column_stack((pressure_alt*0.3048,
                         mach,
                         delta_T,
                         fuelflowleft*0.45359237/3600,
                         fuelflowright*0.45359237/3600))
np.savetxt('matlab.dat', DataOut,  fmt='%10.6f',delimiter="")
