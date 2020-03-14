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





def instantanious_weight(t,bem,block,data, weight_seats,weight_bagage): # contains hardcoded passenger weights

    fuel_use = data[:, 14] + data[:, 15]  # first  is left engine 2nd is right engine fuel use

    fuel_present = block - fuel_use[np.where(t == data[:,0])[0]]

    current_total  =   bem + block + np.sum(weight_seats[:,0]) + np.sum(weight_bagage[:,0]) - fuel_use[np.where(t == data[:,0])[0]]

    zero_fuel_mass = np.sum(weight_seats[:,0]) + np.sum(weight_bagage[:,0]) + bem
    ramp_mass = bem + block + np.sum(weight_seats[:,0]) + np.sum(weight_bagage[:,0])

    return current_total, fuel_present, zero_fuel_mass, ramp_mass

def fuel_moment(t,bem,block,weight_seats,weight_bagage):
    # load the moment data

    path = 'moment_data.csv'
    file = open(path, "r")
    moment_weight = np.genfromtxt(path, delimiter=",", skip_header=1)
    file.close()

    current_weight, fuel_present , zero_fuel,ramp_mass = instantanious_weight(t, bem, block, data, weight_seats, weight_bagage)

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
                fuel_moment = (gradient * current_weight + y_comp)*100

    return fuel_moment

def moment_equilibrium(weight_seats,weight_baggage,current_weight,zero_fuel,ramp_mass,fuel_moment_datum,bem,arm_bem):
    m_seats = np.sum(weight_seats[:,0] * weight_seats[:,1])      # mass in lbs times distance to datum in inch
    m_baggage = np.sum(weight_baggage[:,0] * weight_baggage[:,1])   # mass in lbs times distance to datum in inch
    m_bem = bem*arm_bem

    total_moment_datum = m_seats + m_baggage + m_bem + fuel_moment_datum
    xcg_datum = total_moment_datum / current_weight

    # standard locations
    xcg_datum_bem = m_bem/bem
    xcg_datum_0fuel = (m_seats + m_baggage + m_bem) / zero_fuel
    xcg_datum_ramp = (m_seats+m_baggage+m_bem+fuel_moment_datum)/ramp_mass

    return xcg_datum,xcg_datum_bem,xcg_datum_0fuel,xcg_datum_ramp

# ========================================= MAIN =======================================================================


block = 4100          # in lbs
bem = 9165            # in lbs
arm_bem = 291.65      # in inch
t = 9                 # in sec

weight_seats = np.array([[80, 131], [102, 131], [66, 214], [81, 214], [99, 251], [64, 251], [81, 288], [99, 288], [88, 170]])
weight_seats[:, 0] = 2.20462262 * weight_seats[:, 0]  # convert weight in kg to lbs

# 1st column  = mass 2nd is the distance in inches to front of aircraft
weight_baggage = np.array([[0, 74], [0, 321], [0, 338]])  # baggage set to 0

current_weight,fuel_present, zero_fuel_mass,ramp_mass = instantanious_weight(t, bem, block, data, weight_seats, weight_baggage)
fuel_moment_datum = fuel_moment(t,bem,block,weight_seats,weight_baggage)

xcg_datum = moment_equilibrium(weight_seats,weight_baggage,current_weight,zero_fuel_mass,ramp_mass,fuel_moment_datum,bem,arm_bem)
print(xcg_datum)