# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 09:05:12 2020

@author: Matthew
"""

import scipy.io as spio
import numpy as np

#=====================================================================================


matlab = spio.loadmat('matlab.mat')
time = matlab['flightdata'][0][0][47][0][0][0].transpose()
data = time   
for i in range(47):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis = 1)
    

#======================================================================================



def instantanious_weight(t,bem,block,data): # contains hardcoded passenger weights


    weight_seats = np.array([[10, 131], [20,131], [30,214], [40,214], [50,251], [60,251], [70,288], [80,288], [100,170]])


    # 1st column  = mass 2nd is the distance in inches to front of aircraft
    weight_bagage = np.array([[10,74],[20,321],[30,338]])


    fuel_use = data[:, 14] + data[:, 15]  # first  is left engine 2nd is right engine fuel use

    total  =    bem + block + np.sum(weight_seats[:,0]) + np.sum(weight_bagage[:,0]) - fuel_use[np.where(t == data[:,0])[0]]

    return total



block = 10
bem = 10
t = 10

current_weight = instantanious_weight(t,bem,block,data) # output in pounds
print(current_weight)




