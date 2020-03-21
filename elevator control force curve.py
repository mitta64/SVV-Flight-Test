# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 13:46:10 2020

@author: rikkr
"""

import matplotlib.pyplot as plt
import numpy as np
from Parameters_at_altitude import V_EAS

#ELEVATOR CONTROL FORCE CURVE

#The following arrays define the velocity and the force on the control column

Vias = np.array([152,143,133,122,164,172,183]) #values are given for the indicated airspeed in kts
Fe = np.array([0,-14,-32,-47,28,42,73]) #values are given for the control force in N
Altitude = np.array([13880,14000,14000,14180,13300,12900,12420]) #altitude of the aircraft given in ft
Ttrue = np.array([-14.2,-15,-15.5,-16.5,-12,-10.5,-9.2]) #true air temperature given in degrees C

Vcas = 0.51444 * np.subtract(Vias, 2)

Veas = V_EAS(Vcas, Altitude, Ttrue)

Vseas = sorted(Veas)
Feas = sorted(Fe)

plt.title("Elevator control force curve")
plt.plot(Vseas, Feas)
plt.xlabel("$V_{e} [m/s]$")
plt.ylabel("$F_{e} [N]$")
plt.gca().invert_yaxis()
plt.grid()
#z = np.polyfit(Veas.flatten(), Fe.flatten(), 1)
#p = np.poly1d(z)
#plt.plot(Veas,p(Veas),"r--")
#plt.title("y=%.6fx+%.6f"%(z[0],z[1])) 
plt.show()