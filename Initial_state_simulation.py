import scipy.io as spio
import numpy as np
from Cit_par_phogoid import *
import matplotlib.pyplot as plt
from CG_per_time import mass
'''
matlab = spio.loadmat('matlab.mat')
time = matlab['flightdata'][0][0][47][0][0][0].transpose()
data = time
for i in range(47):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis = 1)

np.savetxt("data.csv", data, delimiter=",")
'''

path = 'data.csv'
file = open(path, "r")
data = np.genfromtxt(path, delimiter=",", skip_header=0)
file.close()

print('loaded')

#A = np.array([[-1,1,0,0],[-8,0,0,0],[1,3,6,5],[1,7,6,7]])
#B = np.array([[2],[2],[0],[1]])
#C = np.array([[1,0],[0,1]])
#D = np.array([[0],[0]])

time = data[:,0]
velocity = data[:,41]
u_flight = data[:,41]-160
AOA = data[:,1]
pitch = data[:,22]
pitchrate = data[:,27]

control_input_sym = np.radians(data[:,17])


#plt.plot(time[31000:35000],velocity[31000:35000],time[31000:35000],pitch[31000:35000])
#plt.grid(True)
#plt.show()


x0_phogoid = np.array([[velocity[32502]-161],[np.radians(AOA[32502])],[np.radians(pitch[32502])],[np.radians(pitchrate[32502])]])

# correct system matrix, asymmetric EOMs


# Asymetric EOM


#print(np.linalg.eig(A_sym))


def initial_repsonse(mode,t,x0,input,mass):
    # if mode == 1 --> symetric EOM used
    # if mode == 0 --> asymetric EOM used

    # A = system matrix
    # B = control matrix
    # C = output matrix
    # D = gain matrix

    # Redefine matices
    #symetric
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

    #asymetric
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

    x = x0
    dt = 0.1

    y1 = []
    y2 = []
    y3 = []
    y4 = []
    length = np.shape(t)[0]

    for i in range(2000):

        y1.append(float(x[0]))
        y2.append(float(x[1]))
        y3.append(float(x[2]))
        y4.append(float(x[3]))

        x_dot = np.dot(A,x) + input[32502+i]*B

        x = x + dt*x_dot

    return y1,y2,y3,y4

y1,y2,y3,y4 = initial_repsonse(1,time,x0_phogoid,control_input_sym,mass)

time = time[0:2000]

plt.plot(time[0:],y1[0:],label='u_numerical')
#plt.plot(time[0:],y2[0:],label='AOA')
#plt.plot(time[0:],y3[0:],label='pitch')
#plt.plot(time[0:],y4[0:],label='pitchrate')

#compare against flight data

plt.plot(time,u_flight[32502:34502],label='u_flight')
plt.plot(time,control_input_sym[32502:34502]*57.2957795,label = 'control_input')
plt.plot(time,y1[0:]-u_flight[32502:34502],label='Absolute Error')
plt.legend()
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