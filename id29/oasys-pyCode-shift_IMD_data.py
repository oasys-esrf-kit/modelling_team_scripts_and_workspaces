import numpy as np
from srxraylib.plot.gol import plot
import matplotlib.pyplot as plt

# plot(*positional_parameters, title="", xtitle="", ytitle="", xrange=None, yrange=None, show=1, legend=None, legend_position=None, color=None, marker=None, linestyle=None, xlog=False, ylog=False, figsize=None)

file = './Pd1.36_B4C0.8_C0.8-11keV-th1.12deg.txt'

data = np.loadtxt(file)


e_shift = 60
data2 = data.copy()
data2[:,0] += e_shift

plot(data[:,0],data[:,1],data2[:,0],data2[:,1], yrange = (0.001,1), ylog = True, legend=['before shift', 'after'])
plt.xlim([10000,12000])

np.savetxt(file[:-4] + f'-shifted_{e_shift:.0f}eV' + '.txt', data2)