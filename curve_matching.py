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
from Cit_par import b, c

matlab = spio.loadmat('FTISxprt-20200306_flight3.mat')
time = matlab['flightdata'][0][0][48][0][0][0].transpose()
data = time
for i in range(48):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis=1)
#===============================================================
def createlist(starttime, endtime, x):
	return np.arange(0, endtime-starttime, x)

#pitch_rate_asym = data[:, 28]
# t_asym = 3205 #[sec], the time at which the asym eigenmotions are conducted
# t_endasym = 3377
# find_peaks(state, height = 0.5, width = 5)[0]
#===============================================================
# yawasym = data[:,29]   #yaw has 29 index in matlabfile
# t_asym = 3515 #[sec], the time at which the asym eigenmotions are conducted
# t_endasym = 3555
#===============================================================
roll_rate_asym = data[:, 27]
# t_asym = 3705 #[sec], the time at which the asym eigenmotions are conducted
# t_endasym = 3730
#===============================================================

###Retrieving the assymetric eiegnmotions###
t_asym = 3705 #[sec], the time at which the asym eigenmotions are conducted
t_endasym = 3730
i = list(data[:,0]).index(t_asym) #finding the index at which the asym motions start
iend = list(data[:,0]).index(t_endasym)

plt.plot(createlist(0, t_endasym-t_asym, 0.1), roll_rate_asym[i:iend], label = 'Pitch rate: flight data')
plt.show()



#print(find_peaks(pitch_rate_asym[i:iend], height = 0.5, width = 5)[0])


V_tas = data[:,43] * 0.514444   # true airspeed in m/s

# Computes only complex eigenvalues  
def Flightdata_eigenvalues(state, t_begin, t_end):
    # Compute the period
    period = (find_peaks(state)[0][2] - find_peaks(state)[0][1]) * 0.1 # The timestep is 0.1, therefore multiply by this
    indices = find_peaks(state)[0]#[1:5] # Leave out the first peak since it is during the input phase

    peaks = []
    for i in indices:
        peaks.append(state[i])
    # print(peaks)
    # print(indices)
    
    time = 0.1 * indices
    # Define the function that needs to fit the data points
    def find_function_peaks(x, a, d):
        return d*np.exp(-a*x)
    
    # Curve_fit finds the coefficients a, b of the defined function to fit the data of the peaks
    a, d = curve_fit(find_function_peaks, time, peaks)[0]
    print(a,d)
    # Obtain values of the exponential function such that a plot can be made
    t = np.arange(time[0], time[-1] + 0.001, 0.001)
    y = []
    for i in t:
        f = d * np.exp(-a * i)
        y.append(f)
    
    # Find the time to damp to half the amplitud: T_half
    T_half = (t[np.where(y <= 0.5 * np.max(y))[0][0]] - time[0]) *0.1
    
    #Indices for V_tas
    index_Vtas_Thalf = t_begin + int(np.rint(t[np.where(y <= 0.5 * np.max(y))[0][0]] * 10))
    index_Vtas_period = t_begin + indices[0]
    #print(index_Vtas_period)
    # Compute parts of eigenvalues
    eig_real = (np.log(0.5) * b) / (T_half * V_tas[index_Vtas_Thalf])
    eig_complex = (2 * np.pi * b) / (period * V_tas[index_Vtas_period])
    eig = complex(eig_real, eig_complex)
    damp = (-eig_real) / (np.sqrt(eig_real**2 + eig_complex**2))
    return T_half, period, damp, y, t, time, peaks, eig

T_half, period, damp, y, t, time, peaks, eig = Flightdata_eigenvalues(roll_rate_asym[i:iend], i, iend)
print(eig, T_half, period, damp)

#Check the exponential function that fits the peaks
plt.scatter(time, peaks)
plt.plot(t, y)
plt.show()
