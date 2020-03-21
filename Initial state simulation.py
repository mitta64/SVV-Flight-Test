import scipy.io as spio
import numpy as np
from Cit_par import *
import matplotlib.pyplot as plt
import control.matlab as control
from CG_per_time import instantanious_weight

matlab = spio.loadmat('matlab.mat')
time = matlab['flightdata'][0][0][47][0][0][0].transpose()
data = time
for i in range(47):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis = 1)


A = np.array([[-1,1],[-8,0]])
B = np.array([[2],[2]])
C = np.array([[1,0],[0,1]])
D = np.array([[0],[0]])
t = np.linspace(10,60,500)
x0 = np.array([[1],[2]])

print(np.linalg.eig(A))

weight_seats = np.array([[80, 131], [102, 131], [66, 214], [81, 214], [99, 251], [64, 251], [81, 288], [99, 288], [88, 170]])
weight_seats[:, 0] = 2.20462262 * weight_seats[:, 0]  # convert weight in kg to lbs
block = 4100  # in lbs
bem = 9165  # in lbs

# compute weight for each time t
time = data[:,0]
mass = instantanious_weight(time,bem, block, data, weight_seats)[0]

def initial_repsonse(A,B,C,D,t,x0,t0,mass):
    # A = system matrix
    # B = control matrix
    # C = output matrix
    # D = gain matrix
    x = x0
    dt = 0.1
    y1 = []
    y2 = []


    t = t0
    for i in range(100):

        y1.append(float(x[0]))
        y2.append(float(x[1]))


        x_dot = np.dot(A,x)
        x = x + dt*x_dot

        # update A
        #A[1][1] = mass[int(t*10)]*0.0001
        t = t+dt
    return y1,y2

y1,y2 = initial_repsonse(A,B,C,D,t,x0,9.0,mass)


x = np.linspace(0,100,100)

plt.plot(x[0:],y1[0:],x[0:],y2[0:])
plt.show()