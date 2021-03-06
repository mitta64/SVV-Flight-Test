import scipy.io as spio
import numpy as np
from Cit_par_phogoid import *
import matplotlib.pyplot as plt
from CG_per_time import instantanious_weight

matlab = spio.loadmat('FTISxprt-20200306_flight3.mat')
time = matlab['flightdata'][0][0][48][0][0][0].transpose()
data = time
for i in range(48):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis=1)

print('loaded')

time = data[:,0]

block = 4100  # in lbs
bem = 9165  # in lbs
arm_bem = 291.65  # in inch
mac = 80.98  # in inch

# seating defined
weight_seats = np.array([[80, 131], [102, 131], [66, 214], [81, 214], [99, 251], [64, 251], [81, 288], [99, 288], [88, 170]])
weight_seats[:, 0] = 2.20462262 * weight_seats[:, 0]  # convert weight in kg to lbs

mass = instantanious_weight(time,bem, block, data, weight_seats)[0]


def initial_repsonse(mode,t0,duration,x0,input,mass,velocity):
    # t0 in sec, accurate to 1 decimal place
    # duration is in sec
    # if mode == 1 --> symetric EOM used
    # if mode == 0 --> asymetric EOM used

    # A = system matrix
    # B = control matrix
    # C = output matrix
    # D = gain matrix

    x = x0
    dt = 0.1

    y1 = []
    y2 = []
    y3 = []
    y4 = []

    inital_index = int(t0*10)

    for i in range(duration*10):
        # Redefine matices
        # symetric
        #V0 = velocity[inital_index+i]
        #mub = mass[inital_index+i] / (rho * S * b)
        #muc = mass[inital_index+i] / (rho * S * c)

        C1 = np.matrix([[-2 * muc * c / V0, 0, 0, 0],
                        [0, (CZadot - 2 * muc) * c / V0, 0, 0],
                        [0, 0, -c / V0, 1],
                        [0, Cmadot * c / V0, 0, -2 * muc * KY2 * c / V0]])
        C1[:, 0] = (1 / V0) * C1[:, 0]
        C1[:, 3] = (c / V0) * C1[:, 3]

        C1_inv = np.linalg.inv(C1)

        C2 = np.matrix([[CXu, CXa, CZ0, CXq],
                        [CZu, CZa, -CX0, CZq + 2 * muc],
                        [0, 0, 0, 1],
                        [Cmu, Cma, 0, Cmq]])
        C2[:, 0] = (1 / V0) * C2[:, 0]
        C2[:, 3] = (c / V0) * C2[:, 3]

        C3 = np.matrix([[-CXde],
                        [-CZde],
                        [0],
                        [-Cmde]])

        A_sym = - C1_inv * C2
        B_sym = - C1_inv * C3

        # asymetric
        B1 = np.matrix([[(CYbdot - 2 * mub) * (b / V0), 0, 0, 0],
                        [0, -b / (2 * V0), 0, 0],
                        [0, 0, -4 * mub * KX2 * (b / V0), 4 * mub * KXZ * (b / V0)],
                        [Cnbdot * (b / V0), 0, 4 * mub * KXZ * (b / V0), -4 * mub * KZ2 * (b / V0)]])
        B1[:, 2] = (b / (2 * V0)) * B1[:, 2]
        B1[:, 3] = (b / (2 * V0)) * B1[:, 3]

        B1_inv = np.linalg.inv(B1)

        B2 = np.matrix([[CYb, CL, CYp, CYr - 4 * mub],
                        [0, 0, 1, 0],
                        [Clb, 0, Clp, Clr],
                        [Cnb, 0, Cnp, Cnr]])
        B2[:, 2] = (b / (2 * V0)) * B2[:, 2]
        B2[:, 3] = (b / (2 * V0)) * B2[:, 3]

        B3 = np.matrix([[-CYda, -CYdr],
                        [0, 0],
                        [-Clda, -Cldr],
                        [-Cnda, -Cndr]])

        A_asym = - B1_inv * B2

        B_asym = - B1_inv * B3

        if mode == 1:
            A = A_sym
            B = B_sym
        else:
            A = A_asym
            B = B_asym

        y1.append(float(x[0]))
        y2.append(float(x[1]))
        y3.append(float(x[2]))
        y4.append(float(x[3]))

        #print(np.linalg.eig(A))
        input_vector = np.array([[input[0,inital_index+i]],[input[1,inital_index+i]]])
        x_dot = np.dot(A,x) + np.dot(B,input_vector)
        x = x + dt*x_dot

    return y1,y2,y3,y4

# ================================ MAIN ============================================================

velocity = data[:,43]

#Sym state vector parameters
AOA = data[:,1]
pitch = data[:,23]
pitchrate = data[:,28]
control_input_sym = data[:,18]

#Getting input in required formcontrol_input_asym


aileron = np.reshape(data[:,17],(52561,1))
rudder = np.reshape(data[:,19],(52561,1))
control_input_asym = np.concatenate((aileron.T,rudder.T),axis=0)


#Asym state vector parameters
#beta = data[:,]
rollangle = data[:,22]
rollrate = data[:,27]
yawrate = data[:,29]

t_initial = 3680 #sec # =3207 for phogoid
v_init = velocity[int(t_initial*10)] # =185 for phogoid
u_flight = data[:,43]-v_init
duration = 40 #sec = 200 for phogoid


x0_sym = np.array([[velocity[int(t_initial*10)]-v_init],[AOA[int(t_initial*10)]],[pitch[int(t_initial*10)]],[pitchrate[int(t_initial*10)]]])
x0_asym = np.array([[0],[rollangle[int(t_initial*10)]],[rollrate[int(t_initial*10)]],[yawrate[int(t_initial*10)]]])
'''
plt.subplot(2, 1, 1)
plt.plot(time[31000:42000],rollangle[31000:42000])
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(time[31000:42000],aileron[31000:42000])
plt.plot(time[31000:42000],rollrate[31000:42000])
#plt.plot(time[31000:42000],data[31000:42000,42]) # is the calibrated airspeed
plt.grid(True)

plt.show()
'''

y1,y2,y3,y4 = initial_repsonse(0,t_initial,duration,x0_asym,control_input_asym,mass,velocity)
shift = 170
errorlist = abs(y3-rollrate[int(t_initial*10):int((t_initial+duration)*10)])
avgerror = sum(abs(errorlist))/len(errorlist)
maxerror = abs(max(errorlist))
print(avgerror, maxerror)

y3 = np.array(y3)[shift:]
y3 += -1.5

time = time[int(t_initial*10+shift):int((t_initial+duration)*10)]
time =np.linspace(0,len(time)/100,len(time))




plt.plot(time,y3,label='Rollrate_numerical')
#plt.plot(time[0:],y2[0:],label='AOA')
#plt.plot(time[0:],y3[0:],label='pitch')
#plt.plot(time[0:],y4[0:],label='pitchrate')

#compare against flight data


plt.plot(time,rollrate[int(t_initial*10+shift):int((t_initial+duration)*10)],label='Rollrate_flight')
#plt.plot(time,control_input_asym[1,int(t_initial*10+shift):int((t_initial+duration)*10)],label = 'rudder_input')
plt.plot(time,control_input_asym[0,int(t_initial*10+shift):int((t_initial+duration)*10)],label = 'Aileron_input')
#plt.plot(time,y1[0:]-u_flight[32502:34502],label='Absolute Error')
plt.xlabel('Time [s]')
plt.ylabel('Rollrate [deg/s] & Aileron input [deg]')
plt.legend()
plt.grid(True)
plt.show()



'''
C1 = np.matrix([[-2 * muc * c / V0, 0, 0, 0],
                [0, (CZadot - 2 * muc) * c / V0, 0, 0],
                [0, 0, -c / V0, 1],
                [0, Cmadot * c / V0, 0, -2 * muc * KY2 * c / V0]])
C1[:, 0] = (1 / V0) * C1[:, 0]
C1[:, 3] = (c / V0) * C1[:, 3]

C1_inv = np.linalg.inv(C1)

C2 = np.matrix([[CXu, CXa, CZ0, CXq],
                [CZu, CZa, -CX0, CZq + 2 * muc],
                [0, 0, 0, 1],
                [Cmu, Cma, 0, Cmq]])
C2[:, 0] = (1 / V0) * C2[:, 0]
C2[:, 3] = (c / V0) * C2[:, 3]

C3 = np.matrix([[-CXde],
                [-CZde],
                [0],
                [-Cmde]])

A_sym = - C1_inv * C2

B1 = np.matrix([[(CYbdot - 2 * mub) * (c / V0), 0, 0, 0],
                [0, -c / (2 * V0), 0, 0],
                [0, 0, -4 * mub * KX2 * (c / V0), 4 * mub * KXZ * (c / V0)],
                [Cnbdot * (c / V0), 0, 4 * mub * KXZ * (c / V0), -4 * mub * KZ2 * (c / V0)]])
B1[:, 2] = (b / (2 * V0)) * B1[:, 2]
B1[:, 3] = (b / (2 * V0)) * B1[:, 3]

B1_inv = np.linalg.inv(B1)

B2 = np.matrix([[CYb, CL, CYp, CYr - 4 * mub],
                [0, 0, 1, 0],
                [Clb, 0, Clp, Clr],
                [Cnb, 0, Cnp, Cnr]])
B2[:, 2] = (b / (2 * V0)) * B2[:, 2]
B2[:, 3] = (b / (2 * V0)) * B2[:, 3]

B3 = np.matrix([[-CYda, -CYdr],
                [0, 0],
                [-Clda, -Cldr],
                [-Cnda, -Cndr]])

A_asym = - B1_inv * B2

print(np.linalg.eig(A_sym))
print(np.linalg.eig(A_asym))
'''