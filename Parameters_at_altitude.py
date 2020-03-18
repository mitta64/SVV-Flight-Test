#Functions to determine the static parpameters at altitude
import numpy as np


def hp0(alt):
	#static coefficients ISA
	lmbda  = -0.0065         # temperature gradient in ISA [K/m]
	Temp0  = 288.15          # temperature at sea level in ISA [K]
	R      = 287.05          # specific gas constant [m^2/sec^2K]
	g      = 9.81            # [m/sec^2] (gravity constant)
	p0     = 101325.0        # pressure at sea level [Pa]
	
	
	hp    = p0 * ((1 + (lmbda * alt)/(Temp0))**((-g)/(lmbda*R)))# pressure altitude in the stationary flight condition [m]
	
	return hp




def mach(V_C, altitude):
 	rho0   = 1.2250          # air density at sea level [kg/m^3] 
 	p0     = 101325.0              # pressure at sea level [Pa]
 	gamma  = 1.4
 	a = 2/(gamma - 1)
 	p = hp0(altitude)
 	b = p0/p	
 	c = (gamma -1)/(2*gamma)
 	d = rho0/p0
 	e = V_C**2
 	f = gamma/(gamma -1)
 	g = (gamma -1)/gamma
 	Mach   = np.sqrt((a)*(((1+b*(((1+c*d*e)**f) - 1))**g)-1))

 	return Mach




def hTemp0(V_Cal, altitude, temp_measured_total):
	gamma  = 1.4
	Mach   = mach(V_Cal, altitude)
	hTemp0 = (273.15 + temp_measured_total) / (1 + ((gamma - 1)/2) * Mach**2)
	
	return hTemp0 



def SoS(V_cal, altitude, temp_measured_total):
	R      = 287.05          # specific gas constant [m^2/sec^2K]
	gamma  = 1.4
	HT     = hTemp0(V_cal, altitude, temp_measured_total)
	a      = np.sqrt(gamma * R * HT)
	
	return a



def V_TAS(V_C, altitude, temp_measured_total):
	M     = mach(V_C, altitude)
	a     = SoS(V_C, altitude, temp_measured_total)
	V_TAS = M*a         # true airspeed in the stationary flight condition [m/sec]
	
	return V_TAS

def Rho(altitude):
	rho0   = 1.2250          # air density at sea level [kg/m^3] 
	lmbda  = -0.0065         # temperature gradient in ISA [K/m]
	Temp0  = 288.15          # temperature at sea level in ISA [K]
	R      = 287.05          # specific gas constant [m^2/sec^2K]
	g      = 9.81            # [m/sec^2] (gravity constant)
	
	rho    = rho0 *((1+(lmbda * altitude / Temp0))**(-((g/(lmbda*R)))-1)) 

	return rho

def V_EAS(V_C, altitude, temp_measured_total):
	rho0 = 1.2250
	Vt   = V_TAS(V_C, altitude, temp_measured_total)	
	r    = Rho(altitude)
	V_E = Vt * np.sqrt((r/rho0))
	
	return V_E



#==============================================================================
#TEST OF FUNCTIONS
#==============================================================================

# hp = hp0(3048)
# print(hp)

# Mach = mach(274.4, 0.0)
# print(Mach)

# htemp0 = hTemp0(274.4, 0, 15)
# print(htemp0)

# a = SoS(274.4, 0, 15)
# print(a)

# V = V_TAS(274.4, 0, 15)
# print(V)
 	
# rho = Rho(2000)
# print(rho)	


