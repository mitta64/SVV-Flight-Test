# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 08:30:03 2020

@author: Thomas
"""

#Imports
from Cit_parBijbelValues import *
from math import *
import numpy as np
import matplotlib.pyplot as plt


#Short period
A = 4*muc**2*KY2
B = -2*muc*(KY2*CZa + Cmadot + Cmq)
C = CZa*Cmq - 2*muc*Cma
print(A,B,C)

#General complex eigenvalue form
Lambda1 = -(B/(2*A)) - 1j*(sqrt(4*A*C - B**2) / (2*A))
Lambda2 = -(B/(2*A)) + 1j*(sqrt(4*A*C - B**2) / (2*A))

correctL1 = -3.9161*10**-2 - 1j*3.7971*10**-2
correctL2 = -3.9161*10**-2 + 1j*3.7971*10**-2

print(Lambda1, Lambda2)
print(correctL1, correctL2)

plt.plot(Lambda1.real, Lambda1.imag, marker="x")
plt.plot(Lambda2.real, Lambda2.imag, marker="x")
plt.plot(correctL1.real, correctL1.imag, marker="o")
plt.plot(correctL2.real, correctL2.imag, marker="o")
plt.show()


