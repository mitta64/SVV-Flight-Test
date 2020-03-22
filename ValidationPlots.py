# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 17:32:11 2020

@author: Thomas
"""

#import
from CG_per_time import data
import numpy as np
import matplotlib.pyplot as plt
from Initial_state_simulation import initial_repsonse, A_asym, x0, mass, A_sym

def createlist(starttime, endtime, x):
	return np.arange(0, endtime-starttime, x)


###Retrieving the assymetric eiegnmotions###
t_asym = 3515 #[sec], the time at which the asym eigenmotions are conducted
t_endasym = 3555
i = list(data[:,0]).index(t_asym) #finding the index at which the asym motions start
iend = list(data[:,0]).index(t_endasym)
#17-19 control surface inputs
da = data[:,17]
de = data[:,18]
dr = data[:,19]

#yaw has 29 index in matlabfile
yawasym = data[:,29]
time = data[:,0]



#Determining the output of the numerical model
#timeused = time[t_asym] - time[t_asym][0] #making sure the function uses the same time as the flight test to calc, starting from 0
timeused = createlist(t_asym, t_endasym, 0.1)
yawratenum = initial_repsonse(A_asym, timeused, x0, mass)

plt.plot(createlist(0, t_endasym-t_asym, 0.1), yawasym[i:iend], label = 'Yaw rate: flight data')
plt.plot(createlist(0, t_endasym-t_asym, 0.1), dr[i:iend], label = 'Rudder input')
plt.plot(timeused, yawratenum[3], label = 'Yaw rate: numerical model') #plotting the yaw rate against the time of the validation data.
plt.legend()
plt.show()

