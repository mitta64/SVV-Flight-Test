# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 16:21:23 2020

@author: Matthew
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.io as spio
from scipy.signal import find_peaks

matlab = spio.loadmat('FTISxprt-20200306_flight3.mat')
time = matlab['flightdata'][0][0][48][0][0][0].transpose()
data = time
for i in range(48):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis=1)

###Retrieving the assymetric eiegnmotions###
t_asym = 3515 #[sec], the time at which the asym eigenmotions are conducted
t_endasym = 3555
i = list(data[:,0]).index(t_asym) #finding the index at which the asym motions start
iend = list(data[:,0]).index(t_endasym)

#yaw has 29 index in matlabfile
yawasym = data[:,29]

V_tas = data[:,43] * 0.514444   # true airspeed in m/s

# Computes only complex eigenvalues  
def Flightdata_eigenvalues(state, t_begin, t_end):
    # Compute the period
    period = (find_peaks(state)[0][2] - find_peaks(state)[0][1]) * 0.1 # The timestep is 0.1, therefore multiply by this
    indices = find_peaks(state)[0][1:7] # Leave out the first peak since it is during the input phase
    peaks = []
    for i in indices:
        peaks.append(state[i])
    
    
    time = 0.1 * indices
    # Define the function that needs to fit the data points
    def find_function_peaks(x, a, b):
        return b*np.exp(a*x)
    
    # Curve_fit finds the coefficients a, b of the defined function to fit the data of the peaks
    a, b = curve_fit(find_function_peaks, time, peaks)[0]
        
    # Obtain values of the exponential function such that a plot can be made
    t = np.arange(time[0], time[5] + 0.01, 0.001)
    y = []
    for i in t:
        f = b * np.exp(a * i)
        y.append(f)
    
    # Find the time to damp to half the amplitud: T_half
    T_half = t[np.where(y <= 0.5 * np.max(y))[0][0]] - time[0]
    
    #Indices for V_tas
    index_Vtas_Thalf = t_begin + int(np.rint(t[np.where(y <= 0.5 * np.max(y))[0][0]] * 10))
    index_Vtas_period = t_begin + indices[0]
    
    # Compute parts of eigenvalues
    eig_real = (np.log(0.5) * b) / (T_half * V_tas[index_Vtas_Thalf])
    eig_complex = (2 * np.pi * b) / (period * V_tas[index_Vtas_period])
    eig = complex(eig_real, eig_complex)
    
    return y, t, time, peaks, eig

y, t, time, peaks, eig = Flightdata_eigenvalues(yawasym[i:iend], i, iend)
print(eig)

#Check the exponential function that fits the peaks
plt.scatter(time, peaks)
plt.plot(t, y)
plt.show()
