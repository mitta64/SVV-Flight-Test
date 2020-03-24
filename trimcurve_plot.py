import numpy as np
from matplotlib import pyplot as plt
import scipy.io as spio
from Cit_par_phogoid import *

data = np.array([[3.9,	0.6],[4.6,0.2],[5.2,0],[6.2,-0.5],[7.1,-0.9],[8.4,-1.5],[10,-2.2],])
fig,ax = plt.subplots()
ax.plot(data[:,0],data[:,1], marker = 'v')
ax.grid(True)
plt.xlabel( 'Angle of Attack [deg]')
plt.ylabel("Elevator deflection [deg]")
ax.axhline(y=0, color='k')
plt.show()


matlab = spio.loadmat('matlab.mat')
time = matlab['flightdata'][0][0][47][0][0][0].transpose()
data = time
for i in range(47):
    data = np.concatenate((data, matlab['flightdata'][0][0][i][0][0][0]), axis = 1)

time = data[32400:33600,0]
time = (time)/60
elevator = data[32400:33600,17]
plt.plot(time,elevator)
plt.xlabel('time[min]')
plt.ylabel('Elevator deflection [deg]')
plt.show()


