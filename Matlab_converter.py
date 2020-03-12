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

weight_seats = np.array([[10, 131], [20,131], [30,214], [40,214], [50,251], [60,251], [70,288], [80,288], [100,170]])
