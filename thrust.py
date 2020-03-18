# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 13:45:58 2020

@author: HJ Hoogendoorn

File for Thrust calculations

Thrust.exe only works for a selected range, e.g. the altitude should at least be 1000 meters and the mach number > 0
"""

from Parameters_at_altitude import *
import numpy as np

#open data file
data = np.genfromtxt('inputthrust.dat')

#reading the columns of interest out of data file
pressure_alt    = data[:,0]
Temp            = data[:,6]
fuelflowleft    = data[:,3]
fuelflowright   = data[:,4]
IAS             = data[:,1]

#deltaT calculations
delta_T = []
for i in range(0,len(Temp)):
    T_ISA = 15.15 -0.0065*(pressure_alt[i]*0.3048)
    delta_T.append(Temp[i]-T_ISA)

IAS = IAS - 2
machn = []
pressure_alt = pressure_alt*0.3048
for i in range(0,len(IAS)):
    machn.append(mach(IAS[i], pressure_alt[i]))


delta_T = np.array(delta_T)
machn = np.array(machn)




#creating .dat file for Thrust.exe
#np.set_printoptions(suppress=True)
DataOut=np.column_stack((pressure_alt*0.3048,
                         machn,
                         delta_T,
                         fuelflowleft*0.45359237/3600,
                         fuelflowright*0.45359237/3600))
np.savetxt('matlab.dat', DataOut,  fmt='%1.6f',delimiter=' ')
