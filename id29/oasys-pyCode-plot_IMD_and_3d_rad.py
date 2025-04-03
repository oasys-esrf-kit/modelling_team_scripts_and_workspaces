import numpy as np
from srxraylib.plot.gol import plot
import matplotlib.pyplot as plt

# plot(*positional_parameters, title="", xtitle="", ytitle="", xrange=None, yrange=None, show=1, legend=None, legend_position=None, color=None, marker=None, linestyle=None, xlog=False, ylog=False, figsize=None)


#file_IMD = './Pd1.36_B4C0.8_C0.8-11keV-th1.12deg.txt'
file_IMD = './Pd1.36_B4C0.8_C0.8-11keV-th1.12deg-shifted_200eV.txt'
data_IMD = np.loadtxt(file_IMD)

#plot(data[:,0],data[:,1], yrange = (0.001,1), ylog = True)

rad = in_object_1.get_contents("xoppy_data")[0]
e = in_object_1.get_contents("xoppy_data")[1]
h = in_object_1.get_contents("xoppy_data")[2]
v = in_object_1.get_contents("xoppy_data")[3]

spec = np.trapz(rad, v)
spec = np.trapz(spec, h)


# Plots
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('energy (eV)')
ax1.set_ylabel('flux (ph/s/0.1%bw', color=color)
ax1.plot(e, spec, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('reflectivity (-)', color=color)  # we already handled the x-label with ax1
ax2.plot(data_IMD[:,0], data_IMD[:,1], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
