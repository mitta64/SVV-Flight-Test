# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#Imports
from Cit_parBijbelValues import *
from math import *
import numpy as np

###Symmetric EoM#######################################################

#Short Period (sp)
"""using simplification a as mentioned in the reader"""

def ShortPeriodSymplified():
    """
    Returns: All Short Period characteristic parameters
    -------
    lambda1: The eigenvalue calced with - discriminant 
    lambda2: The eigenvalue calced with + discriminant 
    Thalf: the time to half damp of the period
    P: period
    Chalf: Amount of periods untill daped to half amplitude 
    delta: logarithmic decrement 
    damp: damping ratio 
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
    """YOU CAN ONLY USE THE PERIOD TO VERIFY THIS MODE, THE AMPLITUDE AND THUS DAMPING GIVES UNACCURATE VALUES"""
    """
    Returns: All Short Period characteristic parameters
    -------
    lambda1: The eigenvalue calced with - discriminant 
    lambda2: The eigenvalue calced with + discriminant 
    Thalf: the time to half damp of the period
    P: period
    Chalf: Amount of periods untill daped to half amplitude 
    delta: logarithmic decrement 
    damp: damping ratio 
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

    
###Asymmetric Eigenmotions#############################################

#Heavily damped aperiodic rolling motion
def HeavilyDampedApriodicRoll():    
    """ 
    Returns: All Short Period characteristic parameters
    -------
    lambdaha: The eigenvalue of the heavily damped aperiod roll,  calced with - discriminant 
    Thalf: the time to half damp of the period
    tao: time constant of this eigenmode
    """
    lambdaha = Clp / (4*mub*KX2)
    
    #Time to half damped
    Thalf = (np.log(0.5) / lambdaha)*(b/V0)
    
    #Time constant
    tao = -(1/lambdaha)*(b/V0)
    
    return lambdaha, Thalf, tao

#Dutch Roll
def DutchRoll():
    """
    Returns: All Short Period characteristic parameters
    -------
    lambda1: The eigenvalue calced with - discriminant 
    lambda2: The eigenvalue calced with + discriminant 
    Thalf: the time to half damp of the period
    P: period
    Chalf: Amount of periods untill daped to half amplitude
    delta: logarithmic decrement 
    damp: damping ratio
    """
    A = 8*mub**2*KZ2
    B = -2*mub*(Cnr + 2*KZ2*CYb)
    C = 4*mub*Cnb + CYb*Cnr
    
    lambda1 = (-B - (np.sqrt(4*A*C - B**2))*1j ) / (2*A)
    lambda2 = (-B + (np.sqrt(4*A*C - B**2))*1j ) / (2*A)
    
    #half time to damp
    Thalf = (np.log(0.5)/lambda1.real)*(b/V0)
    
    #Period
    P = (2*np.pi/lambda2.imag)*(b/V0)
    
    #Number of periods to half damped amplitude
    Chalf = (np.log(0.5)/2*np.pi)*(lambda2.imag/lambda2.real)
    
    #Logarithmic decrement
    delta = 2*np.pi*(lambda1.real/lambda1.imag)
    
    #damping ratio
    damp = -lambda1.real / (np.sqrt(lambda1.real**2 + lambda1.imag**2))
    
    
    return lambda1, lambda2, Thalf, P, Chalf, delta, damp

def ApriodicSpiral():  
    """DONT USE TO VERIFY EIGENVALUES, AS THE LAMBDA IS POS HERE, BUT SHOULD BE NEGATIVE"""
    """
    Returns: All Short Period characteristic parameters
    -------
    lambdaha: The eigenvalue of the heavily damped aperiod roll,  calced with - discriminant 
    Thalf: the time to half damp of the period
    tao: time constant of this eigenmode
    """
    
    lambdaspiral = (2*CL*(Clb*Cnr - Cnb*Clr)) / (Clp*(CYb*Cnr + 4*mub*Cnb) - Cnp*(CYb*Clr + 4*mub*Clb))
    
    #Time to half damped
    Thalf = (np.log(0.5) / lambdaspiral)*(b/V0)
    
    #Time constant
    tao = -(1/lambdaspiral)*(b/V0)
    
    return lambdaspiral, Thalf, tao

def DutchandAperiodRoll():
    """    
    Returns: All Short Period characteristic parameters
    -------
    lambda1: The eigenvalue of the roll of this mode 
    lambda2: The eigenvalue of the dutch roll, calced with - discriminant 
    lambda3: The eigenvalue of the dutch roll, calced with + discriminant 
    Thalf: the time to half damp of the period of the rolling motion
    Thalfdutch: the time to half damp of the period of the Dutch roll
    P: period of the period Dutch roll
    Chalf: Amount of periods untill daped to half amplitude for Dutch roll
    delta: logarithmic decrement for Dutch roll
    damp: damping ratio, Dutch Roll
    
    """
    A = 4*mub**2*(KX2*KZ2 - KXZ**2)
    B = -mub*((Clr + Cnp)*KXZ + Cnr*KX2 + Clp*KZ2)
    C = 2*mub*(Clb*KXZ + Cnb*KX2) + 0.25*(Clp*Cnr - Cnp*Clr)
    D = 0.5*(Clb*Cnp - Cnb*Clp)
    
    eigenvals = np.roots((A,B,C,D))
    print(eigenvals, len(eigenvals))
    lambda1 = eigenvals[0]
    lambda2 = eigenvals[2]
    lambda3 = eigenvals[1]
    
    #Time to half damped
    Thalfroll = (np.log(0.5) / lambda1)*(b/V0)
    
    #half time to damp of imaginary eigenvalues
    Thalfdutch = (np.log(0.5)/lambda2.real)*(b/V0)
    
    #Period
    P = (2*np.pi/lambda2.imag)*(b/V0)
    
    #Number of periods to half damped amplitude
    Chalf = (np.log(0.5)/2*np.pi)*(lambda2.imag/lambda2.real)
    
    #Logarithmic decrement
    delta = 2*np.pi*(lambda1.real/lambda2.imag)
    
    #damping ratio
    damp = -lambda2.real / (np.sqrt(lambda2.real**2 + lambda2.imag**2))
    
    
    return lambda1, lambda2, lambda3, Thalfroll, Thalfdutch, P, Chalf, delta, damp
    
    
        




