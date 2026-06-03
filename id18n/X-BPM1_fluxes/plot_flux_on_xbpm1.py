from srxraylib.plot.gol import plot
import numpy
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import host_subplot


e = in_object_1.get_content('xoppy_data')[:,0]
flux = in_object_1.get_content('xoppy_data')[:,1]

e1 = in_object_2.get_content('xoppy_data')[:,0]
r1 = in_object_2.get_content('xoppy_data')[:,2]
#r2 = in_object_3.get_content('xoppy_data')[:,2]
#r3 = in_object_4.get_content('xoppy_data')[:,2]


# Plots
host = host_subplot(111)
par = host.twinx()

host.set_xlabel("Energy (eV)")
host.set_ylabel("Flux (ph/s/0.1%bw)")
par.set_ylabel("Reflectivity")

p1, = host.plot(e, flux, label="Flux")
p2, = par.plot(e1, r1, label="R1")
#p3, = par.plot(e, r2, label="R2")
#p4, = par.plot(e, r3, label="R3")

host.legend(labelcolor="linecolor")

host.yaxis.label.set_color(p1.get_color())
par.yaxis.label.set_color(p2.get_color())

plt.show()