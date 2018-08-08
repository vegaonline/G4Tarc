import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.animation as anim

indata = np.loadtxt("neutSpec.dat")
print(indata.shape)

timecol = indata[:,[0]]
print(timecol.shape)
energycol = indata[:,[1]]
print(energycol.shape)

print(np.isnan(energycol).any())
print(not np.isfinite(energycol).any())

#plt.hist(energycol, normed=True, bins=100)
