# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#Imports
from Cit_parBijbelValues import *
#from Matlab_converter import data
from math import *
import numpy as np
#import matplotlib.pyplot as plt

###Symmetric EoM###

#Short Period (sp)
"""using simplification a as mentioned in the reader"""
#The time array from the flight data
# = data[:,0]

def ShortPeriodSymplified():
    """
    Parameters
    ----------
    -
    
    Returns: All Short Period characteristic parameters
    -------
    lambda1: The eigenvalue calced with - discriminant 
    lambda2: The eigenvalue calced with + discriminant 
    Thalf: the time to half damp of the period
    P: period
    Chalf: Amount of periods untill daped to half amplitude for lambda 1
    delta: logarithmic decrement lambda1
    damp: damping ratio lambda1
    
    """
    #Calculation of eigenvalues
    A = 4*muc**2*KY2
    B = -2*muc*(KY2*CZa + Cmadot + Cmq)
    C = CZa*Cmq - 2*muc*Cma
    
    lambda1 = (-B - (np.sqrt(4*A*C - B**2))*1j ) / (2*A)
    lambda2 = (-B + (np.sqrt(4*A*C - B**2))*1j ) / (2*A)
    
    #half time to damp
    Thalf = (np.log(0.5)/lambda1.real)*(c/V0)
    
    #Period
    P = (2*np.pi/lambda2.imag)*(c/V0)
    
    #Number of periods to half damped amplitude
    Chalf = (np.log(0.5)/2*np.pi)*(lambda2.imag/lambda2.real)
    
    #Logarithmic decrement
    delta = 2*np.pi*(lambda1.real/lambda1.imag)
    
    #damping ratio
    damp = -lambda1.real / (np.sqrt(lambda1.real**2 + lambda1.imag**2))
    
    
    return lambda1, lambda2, Thalf, P, Chalf, delta, damp

#Phugoid
def PhugoidSymplified():
    """
    Parameters
    ----------
    -
    
    Returns: All Short Period characteristic parameters
    -------
    lambda1: The eigenvalue calced with - discriminant 
    lambda2: The eigenvalue calced with + discriminant 
    Thalf: the time to half damp of the period
    Chalf1: Amount of periods untill daped to half amplitude for lambda 1
    Chalf2: "" "" for lambda2
    delta1: logarithmic decrement lambda1
    delta2: logarithmic decrement lambda2
    damp1: damping ratio lambda1
    damp2: damping ratio lambda2
    
    """
    #Calculation of eigenvalues
    A = 2*muc*(CZa*Cmq - 2*muc*Cma)
    B = 2*muc*(CXu*Cma - Cmu*CXa) + Cmq*(CZu*CXa - CXu*CZa)
    C = CZ0*(Cmu*CZa - CZu*Cma)
    
    lambda1 = (-B - (np.sqrt(4*A*C - B**2))*1j ) / (2*A)
    lambda2 = (-B + (np.sqrt(4*A*C - B**2))*1j ) / (2*A)
    
    #half time to damp
    Thalf = (np.log(0.5)/lambda1.real)*(c/V0)
    
    #Period
    P = (2*np.pi/lambda2.imag)*(c/V0)
    
    #Number of periods to half damped amplitude
    Chalf = (np.log(0.5)/2*np.pi)*(lambda2.imag/lambda2.real)
    
    #Logarithmic decrement
    delta = 2*np.pi*(lambda1.real/lambda1.imag)
    
    #damping ratio
    damp = -lambda1.real / (np.sqrt(lambda1.real**2 + lambda1.imag**2))
    
    
    return lambda1, lambda2, Thalf, P, Chalf, delta, damp

    
###Asymmetric Eigenmotions###

#Heavily damped aperiodic rolling motion
def HeavilyDampedApriodicRoll():
    lambdaha = Clp / (4*mub*KX2)
    
    #Time to half damped
    Thalf = (np.log(0.5) / lambdaha)*(c/60)
    
    #Time constant
    tao = -(1/lambdaha)*(c/60)
    
    return lambdaha, Thalf, tao

#Dutch Roll
def DutchRoll():
    #A = 8*mub**2*KZ2
    #B = -2*mub*(Cnr + 2*KZ2*CYb)
    #C = 4*mub*Cnb + CYb*Cnr
    
    A = -2*mub*KZ2
    B = 0.5*Cnr
    C = -Cnb
    
    lambda1 = (-B - (np.sqrt(4*A*C - B**2))*1j ) / 2*A
    lambda2 = (-B + (np.sqrt(4*A*C - B**2))*1j ) / 2*A
    
    #half time to damp
    Thalf = -(np.log(0.5)/lambda1.real)
    
    #Number of periods to half damped amplitude
    Chalf1 = -(np.log(0.5)/2*np.pi)*(lambda1.imag/lambda1.real)
    Chalf2 = -(np.log(0.5)/2*np.pi)*(lambda2.imag/lambda2.real)
    
    #Logarithmic decrement
    delta1 = 2*np.pi*(lambda1.real/lambda1.imag)
    delta2 = 2*np.pi*(lambda2.real/lambda2.imag)
    
    #damping ratio
    damp1 = -lambda1.real / np.sqrt(lambda1.real**2 + lambda1.imag**2)
    damp2 = -lambda1.real / np.sqrt(lambda2.real**2 + lambda2.imag**2)
    
    return lambda1, lambda2, Thalf, Chalf1, Chalf2, delta1, delta2, damp1, damp2
    
    
        




