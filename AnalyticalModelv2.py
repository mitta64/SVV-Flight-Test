# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 08:30:03 2020

@author: Thomas
"""

#Imports
from Cit_parold import *
from math import *
import numpy as np


#Short period
A = 4*muc**2*KY2
B = -2*muc*(KY2*CZa + Cmadot + Cmq)
C = CZa*Cmq - 2*muc*Cma
print(A,B,C)

#General complex eigenvalue form
Lambda1 = -(B/(2*A)) - 1j*(sqrt(4*A*C - B**2) / (2*A))
Lambda2 = -(B/(2*A)) + 1j*(sqrt(4*A*C - B**2) / (2*A))

print(Lambda1, Lambda2)


