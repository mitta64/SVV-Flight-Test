"""
program determines the cg position with with respect to time

"""

import scipy.io as spio
import numpy as np
from Parameters_at_altitude import *
import matplotlib.pyplot as plt

# =====================================================================================


matlab = spio.loadmat('FTISxprt-20200306_flight3.mat')
time = matlab['flightdata'][0][0][48][0][0][0].transpose()
data = time
for i in range(48):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis=1)


# ======================================================================================


def instantanious_weight(t, bem, block, data, weight_seats):  # contains hardcoded passenger weights

    fuel_use = data[:, 14] + data[:, 15]  # first  is left engine 2nd is right engine fuel use

    fuel_present = block - fuel_use[np.where(t == data[:, 0])[0]]

    current_total = bem + block + np.sum(weight_seats[:, 0]) - fuel_use[np.where(t == data[:, 0])[0]]
    
    zero_fuel_mass = np.sum(weight_seats[:, 0]) + bem
    ramp_mass = bem + block + np.sum(weight_seats[:, 0])

    return current_total, fuel_present, zero_fuel_mass, ramp_mass


# def weight_array(time, bem, block, data, weight_seats):
    
#     fuel_use = data[:, 14] + data[:, 15]  # first  is left engine 2nd is right engine fuel use   
#     j = 0
#     weight = np.empty(np.size(time))
#     for i in time:   
        
#         current_total = bem + block + np.sum(weight_seats[:, 0]) - fuel_use[np.where(i == data[:, 0])[0]]
#         weight[j ] = current_total
#         j += 1
        
    
#     return weight

def fuel_moment(t, bem, block, weight_seats):
    # load the moment data

    path = 'moment_data.csv'
    file = open(path, "r")
    moment_weight = np.genfromtxt(path, delimiter=",", skip_header=1)
    file.close()

    moment_weight[:, 1] = 100 * moment_weight[:, 1]

    current_weight, fuel_present, zero_fuel, ramp_mass = instantanious_weight(t, bem, block, data, weight_seats)

    fuel_moment = 0

    for i in range(np.shape(moment_weight)[0]):
        if i == (np.shape(moment_weight)[0] - 1):
            a = 0

        else:
            if fuel_present >= moment_weight[i, 0] and fuel_present <= moment_weight[i + 1, 0]:
                x1 = moment_weight[i, 0]
                x2 = moment_weight[i + 1, 0]
                y1 = moment_weight[i, 1]
                y2 = moment_weight[i + 1, 1]

                gradient = (y2 - y1) / (x2 - x1)
                y_comp = y1 - gradient * x1

                fuel_moment = gradient * fuel_present + y_comp

    return fuel_moment


def moment_equilibrium(weight_seats, current_weight, zero_fuel, ramp_mass, fuel_moment_datum, bem, arm_bem, mac):
    m_seats = np.sum(weight_seats[:, 0] * weight_seats[:, 1])  # mass in lbs times distance to datum in inch
    m_bem = bem * arm_bem
    m_block = 11705.50 * 100  # From the table read of at 4100 lbs

    total_moment_datum = m_seats + m_bem + fuel_moment_datum
    xcg_datum = total_moment_datum / current_weight

    # standard locations
    xcg_datum_bem = m_bem / bem
    xcg_datum_0fuel = (m_seats + m_bem) / zero_fuel
    xcg_datum_ramp = (m_seats + m_bem + m_block) / ramp_mass

    # cg normalized wrt MAC
    xcg = (xcg_datum - 261.45) * 0.0254 #xcg in distance from start MAC
    #xcg = (xcg_datum - 261.45) * 0.0254 * 100 / (mac * 0.0254)  # xcg is in %c MAC

    return xcg_datum_bem, xcg_datum_0fuel, float(xcg_datum_ramp), float(xcg_datum), float(xcg)



# ========================================= MAIN =======================================================================


block = 4100  # in lbs
bem = 9165  # in lbs
arm_bem = 291.65  # in inch
mac = 80.98  # in inch
t = 3257.8  # in sec

# seating defined
weight_seats = np.array(
    [[80, 131], [102, 131], [66, 214], [81, 214], [99, 251], [64, 251], [81, 288], [99, 288], [88, 170]])
weight_seats[:, 0] = 2.20462262 * weight_seats[:, 0]  # convert weight in kg to lbs

# weight determination
current_weight, fuel_present, zero_fuel_mass, ramp_mass = instantanious_weight(t, bem, block, data, weight_seats)

mass = instantanious_weight(time,bem, block, data, weight_seats)[0]
print(mass[32302])
# moment determined wrt datum caused by fuel consumption
fuel_moment_datum = fuel_moment(t, bem, block, weight_seats)

# cg contains: BEM cg wrt datum, Zero fuel cg wrt datum, Ramp cg wrt datum, current cg wrt datum, cg wrt MAC
cg = moment_equilibrium(weight_seats, current_weight, zero_fuel_mass, ramp_mass, fuel_moment_datum, bem, arm_bem, mac)
cg1 = cg[-1]


# cg shift R3 to 131 inch

t = 3164

weight_seats = np.array([[80, 131], [102, 131], [66, 214], [81, 214], [99, 251], [64, 251], [81, 288], [99, 131], [88, 170]])
weight_seats[:, 0] = 2.20462262 * weight_seats[:, 0]  # convert weight in kg to lbs

current_weight, fuel_present, zero_fuel_mass, ramp_mass = instantanious_weight(t, bem, block, data, weight_seats)
fuel_moment_datum = fuel_moment(t, bem, block, weight_seats)
cg = moment_equilibrium(weight_seats, current_weight, zero_fuel_mass, ramp_mass, fuel_moment_datum, bem, arm_bem, mac)
cg2 = cg[-1]

d_cg = (cg2-cg1)

d_elevator = -0.4 #deg

mac = 2.0569
s = 30.0


v_tas = V_TAS((155-2)*0.514444444,3919.728,-12.2)
rho = Rho(3919.728)


d_elevator = d_elevator*(np.pi/180)
C_N = current_weight/(0.5*rho*v_tas**2*s)

cm_de = -(1/d_elevator)*C_N*(d_cg/mac)


# Cm_alpha determination

slope_trimcurve = -0.4562 # from excel plot
cm_alpha = -cm_de*slope_trimcurve
print(cm_de)
print(cm_alpha)
#
# ======================================================================================
# Plot vtas over t, eigenmotions
# VTAS = data[:,43][29911:40112]
# t = np.arange(3000, 4020.1, 0.1)
# plt.plot(t, VTAS)
# plt.show()