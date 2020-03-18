# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:02:14 2020

@author: HJ Hoogendoorn

Check weight during phugoid and short period

Phugoid 0:54
Short Period 0:56
"""

from Matlab_converter import data

time    = data[:,0]
f_usedl = data[:,14]*0.45359237 #convert to kg
f_usedr = data[:,15]*0.45359237 #convert to kg
f_flowl = data[:,4]*0.45359237/3600 #convert to kg/s
f_flowr = data[:,5]*0.45359237/3600 #convert to kg/s



W0 = 14937*0.45359237 #kg
Weight_flow = []
for i in range(0,len(time)):
    Weight_flow.append(W0 - (f_flowl[i] + f_flowr[i])/10)
    W0 = W0-(f_flowl[i]+f_flowr[i])/10
    
W0 = 14937*0.45359237 #kg
Weight = []
for i in range(0,len(time)):
    '''
    Phugoid 54-56
    Short Period 56-57
    '''
    if time[i] > 56*60 and time[i] < 57*60: #for short period
        Weight.append(W0-f_usedl[i] - f_usedr[i])

        
W_avg = Weight[0] + (Weight[-1] - Weight[0])/2
        
    
