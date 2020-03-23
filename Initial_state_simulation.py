import scipy.io as spio
import numpy as np
from Cit_par_phogoid import *
import matplotlib.pyplot as plt
from CG_per_time import mass

matlab = spio.loadmat('FTISxprt-20200306_flight3.mat')
time = matlab['flightdata'][0][0][48][0][0][0].transpose()
data = time
for i in range(48):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis=1)

print('loaded')

#A = np.array([[-1,1,0,0],[-8,0,0,0],[1,3,6,5],[1,7,6,7]])
#B = np.array([[2],[2],[0],[1]])
#C = np.array([[1,0],[0,1]])
#D = np.array([[0],[0]])



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
        V0 = velocity[32308+i]
        mub = mass[32308+i]/ (rho * S * b)
        muc = mass[32308+i] / (rho * S * c)

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
        B_sym = C1_inv * C3

        # asymetric
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
        B_asym = B1_inv * B3

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

        x_dot = np.dot(A, x) + input[inital_index + i]*B
        #x_dot = np.dot(A,x) + np.dot(B,np.array([[input[inital_index+i]]]))

        x = x + dt*x_dot

    return y1,y2,y3,y4

# ================================ MAIN ============================================================

time = data[:,0]
velocity = data[:,43]

AOA = data[:,1]
pitch = data[:,23]
pitchrate = data[:,28]
control_input_sym = data[:,18]
print("elevator",control_input_sym[0])
control_input_sym = np.radians(control_input_sym)


t_initial = 3415.3 #sec # =3207 for phogoid
v_init = velocity[int(t_initial*10)] # =185 for phogoid
print(v_init)
u_flight = data[:,43]-v_init
duration = 200 #sec = 200 for phogoid


x0= np.array([[0],[np.radians(AOA[int(t_initial*10)])],[np.radians(pitch[int(t_initial*10)])],[np.radians(pitchrate[int(t_initial*10)])]])
'''
plt.subplot(2, 1, 1)
plt.plot(time[31000:42000],pitch[31000:42000])
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(time[31000:42000],data[31000:42000,43])
#plt.plot(time[31000:42000],data[31000:42000,42]) # is the calibrated airspeed
plt.grid(True)

plt.show()
'''
y1,y2,y3,y4 = initial_repsonse(1,t_initial,duration,x0,control_input_sym,mass,velocity)

time = time[0:duration*10]

plt.plot(time,y4,label='u_numerical')
#plt.plot(time[0:],y2[0:],label='AOA')
#plt.plot(time[0:],y3[0:],label='pitch')
#plt.plot(time[0:],y4[0:],label='pitchrate')

#compare against flight data

plt.plot(time,pitchrate[int(t_initial*10):int((t_initial+duration)*10)],label='pitch_flight')
#plt.plot(time,control_input_sym[int(t_initial*10):int((t_initial+duration)*10)]*57.2957795,label = 'control_input')
#plt.plot(time,y1[0:]-u_flight[32502:34502],label='Absolute Error')
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