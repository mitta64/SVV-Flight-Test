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

x0 = np.array([[4],[0.1]])

weight_seats = np.array([[80, 131], [102, 131], [66, 214], [81, 214], [99, 251], [64, 251], [81, 288], [99, 288], [88, 170]])
weight_seats[:, 0] = 2.20462262 * weight_seats[:, 0]  # convert weight in kg to lbs
block = 4100  # in lbs
bem = 9165  # in lbs

# compute weight for each time t
time = data[:,0]
mass = instantanious_weight(time,bem, block, data, weight_seats)[0]
#velocity = data[:,13]

def initial_repsonse(A,B,C,D,t,x0,mass):
    # A = system matrix
    # B = control matrix
    # C = output matrix
    # D = gain matrix
    x = x0
    dt = 0.1
    y1 = []
    y2 = []
    length = np.shape(t)[0]

    for i in t:

        y1.append(float(x[0]))
        y2.append(float(x[1]))


        x_dot = np.dot(A,x)
        x = x + dt*x_dot

        # update A
        A[1][1] = mass[int(np.where(data[:,0] == i)[0])]*0.00001

    return y1,y2

y1,y2 = initial_repsonse(A,B,C,D,time,x0,mass)




plt.plot(time[0:20],y1[0:20],time[0:20],y2[0:20])
plt.show()