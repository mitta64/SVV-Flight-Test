# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 17:32:11 2020

@author: Thomas
"""

#import
from Matlab_converterold import data
import numpy as np
import matplotlib.pyplot as plt


###Retrieving the assymetric eiegnmotions###
t_asym = 57*60 #[sec], the time at which the asym eigenmotions are conducted
t_endasym = 3739
i = list(data[:,0]).index(t_asym) #finding the index at which the asym motions start
iend = list(data[:,0]).index(t_endasym)
#16-18 control surface inputs
da = data[:,17]
de = data[:,18]
dr = data[:,19]

#yew has 29 index in matlabfile
yawasym = data[:,29]
time = data[:,0]

#plt.plot(time[i:iend], yawasym[i:iend])
#plt.plot(time[i:iend], dr[i:iend])
plt.plot(time[i::], yawasym[i::])
plt.plot(time[i::], dr[i::])
plt.show()

