# Citation 550 - Linear simulation
import numpy as np


# Stationary flight condition

#hp0    =    13000*0.3048	      # pressure altitude in the stationary flight condition [m]
V0     =    59.9         # true airspeed in the stationary flight condition [m/sec]
alpha0 =    0         # angle of attack in the stationary flight condition [rad]
th0    =    0         # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      =    4547.8      # mass [kg]

# aerodynamic properties
#e      =     1        # Oswald factor [ ]
#CD0    =     1        # Zero lift drag coefficient [ ]
#CLa    =     1        # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    =    -0.4300       # longitudinal stabilty [ ]
Cmde   =    -1.5530       # elevator effectiveness [ ]

# Aircraft geometry

S      = 24.2	          # wing area [m^2]
#Sh     = 0.2 * S         # stabiliser area [m^2]
#Sh_S   = Sh / S	          # [ ]
lh     = 5.5   # tail length [m]
c      = 2.022	          # mean aerodynamic cord [m]
xcg    = 0.3 * c
lh_c   = lh / c	          # [ ]
b      = 13.36	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
#A      = b ** 2 / S      # wing aspect ratio [ ]
#Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
#Vh_V   = 1	          # [ ]
#ih     = -2 * np.pi / 180   # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity

rho0   = 1.2250          # air density at sea level [kg/m^3] 
lmbda = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

# air density [kg/m^3]  
#rho    = rho0 *(((1+(lmbda * hp0 / Temp0)))**(-((g / (lmbda*R)))))   
W      = m * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia

muc    = 102.7
mub    = 15.5
KX2    = 0.012
KZ2    = 0.037
KXZ    = 0.002
KY2    = 0.980

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
#CNwa   = CLa                    # Wing normal force slope [ ]
#CNha   = 2 * np.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
#depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

CL = 1.1360              # Lift coefficient [ ]
#CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e) # Drag coefficient [ ]

# Stabiblity derivatives

CX0    = 0
CXu    = -0.2199
CXa    = +0.4653		# Positive! (has been erroneously negative since 1993) 
CXadot = 0
CXq    = 0
CXde   = 0

CZ0    = -1.1360
CZu    = -2.2720
CZa    = -5.16
CZadot = -1.43
CZq    = -3.86
CZde   = -0.6238

Cmu    = 0
Cmadot = -3.7
Cmq    = -7.04

CYb    = -0.9896
CYbdot =  0     
CYp    = -0.087
CYr    = +0.4300
CYda   = 0
CYdr   = +0.3037

Clb    = -0.0772
Clp    = -0.3444
Clr    = +0.28
Clda   = -0.2349
Cldr   = +0.0286

Cnb    =  +0.1638
Cnbdot =   0     
Cnp    =  -0.0108
Cnr    =  -0.1930
Cnda   =  +0.0286
Cndr   =  -0.1261