import numpy as np
from matplotlib import pyplot as plt

data = np.array([[3.9,	0.6],[4.6,0.2],[5.2,0],[6.2,-0.5],[7.1,-0.9],[8.4,-1.5],[10,-2.2],])
fig,ax = plt.subplots()
ax.plot(data[:,0],data[:,1], marker = 'v')
ax.grid(True)
plt.xlabel( 'Angle of Attack [deg]')
plt.ylabel("Elevator deflection [deg]")
ax.axhline(y=0, color='k')
plt.show()