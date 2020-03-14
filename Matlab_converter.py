# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 09:05:12 2020

@author: Matthew
"""

import scipy.io as spio
import numpy as np

#=====================================================================================


matlab = spio.loadmat('FTISxprt-20200306_flight3.mat')
time = matlab['flightdata'][0][0][48][0][0][0].transpose()
data = time   
for i in range(48):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis = 1)
    

#======================================================================================





def instantanious_weight(t,bem,block,data): # contains hardcoded passenger weights


    weight_seats = np.array([[80, 131], [102,131], [66,214], [81,214], [99,251], [64,251], [81,288], [99,288], [88,170]])
    weight_seats[:,0] = 2.20462262 * weight_seats[:,0] # convert weight in kg to lbs

    # 1st column  = mass 2nd is the distance in inches to front of aircraft
    weight_bagage = np.array([[0,74],[0,321],[0,338]]) # bagage set to 0


    fuel_use = data[:, 14] + data[:, 15]  # first  is left engine 2nd is right engine fuel use

    fuel_present = block - fuel_use[np.where(t == data[:,0])[0]]

    total  =   bem + block + np.sum(weight_seats[:,0]) + np.sum(weight_bagage[:,0]) - fuel_use[np.where(t == data[:,0])[0]]

    return total, fuel_present

block = 4100          # in lbs
bem = 9165            # in lbs
t = 10              # in sec



# load the moment data

path = 'moment_data.csv'
file = open(path, "r")
moment_weight = np.genfromtxt(path, delimiter=",", skip_header=1)
file.close()

current_weight , fuel_present = instantanious_weight(t,bem,block,data)


moment = 0
for i in range(np.shape(moment_weight)[0]):
    if i == (np.shape(moment_weight)[0]-1):
        a=0

    else:
        if fuel_present >= moment_weight[i,0] and fuel_present <= moment_weight[i+1,0]:
            x1 = moment_weight[i,0]
            x2 = moment_weight[i+1,0]
            y1 = moment_weight[i,1]
            y2 = moment_weight[i+1,1]

            gradient = (y2-y1)/(x2-x1)
            y_comp = y1-gradient*x1
            moment = gradient*current_weight+y_comp

print(moment)



