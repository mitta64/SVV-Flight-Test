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
t_asym = 3600 + 1*60 #[sec], the time at which the asym eigenmotions are conducted
t_endasym = 3600 + 7*60
i = list(data[:,0]).index(t_asym) #finding the index at which the asym motions start
iend = list(data[:,0]).index(t_endasym)
#27-29 control surface inputs
da = data[:,17]
de = data[:,18]
dr = data[:,19]

yawasym = data[:,29]
time = data[:,0]

plt.plot(time, yawasym)
plt.plot(time[i::], dr[i::])
plt.show()

