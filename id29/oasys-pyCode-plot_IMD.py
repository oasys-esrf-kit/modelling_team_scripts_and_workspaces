import numpy as np
from srxraylib.plot.gol import plot
import matplotlib.pyplot as plt

# plot(*positional_parameters, title="", xtitle="", ytitle="", xrange=None, yrange=None, show=1, legend=None, legend_position=None, color=None, marker=None, linestyle=None, xlog=False, ylog=False, figsize=None)

file = './Pd1.36_B4C0.8_C0.8-11keV-th1.12deg.txt'

data = np.loadtxt(file)

plot(data[:,0],data[:,1], yrange = (0.001,1), ylog = True)
plt.xlim([10000,12000])
