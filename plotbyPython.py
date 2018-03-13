#!/bin/python
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import pylab as p
import matplotlib
import matplotlib.pyplot as plt

indata    = np.loadtxt("ETVirtual.dat")
print(indata.shape) 
boxCol    = indata[:,[0]]
timeCol   = indata[:,[1]]
energyCol = indata[:,[2]]

fig = p.figure()
ax=p3.Axes3D(fig)
#ax.plot_wireframe(boxCol, timeCol, energyCol)
#ax.plot_surface(boxCol, timeCol, energyCol)
#boxCol, timeCol = p.meshgrid(boxCol, timeCol)
#ax.contourf3D(boxCol, timeCol, energyCol)

ax.set_xlabel('vbox')
ax.set_ylabel('Time')
ax.set_zlabel('Energy')
#fig.add_axes(ax)
p.show()

