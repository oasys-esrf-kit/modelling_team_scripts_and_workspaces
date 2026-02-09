import os
import copy
import sys
import math
import codecs
from string import *
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm


inESi = open("HistoSilicon.dat","r")
ESi = inESi.readlines()[0:]
inESi.close()
del inESi
if len(ESi)==0:
	print ("There is no absorbtion in the silicon")
else:
	BIN = int(input("Number of bins in the histogram: "))
	vmin = float(input("min. value of colorscale: "))
	vmax = float(input("max. value of colorscale: "))
	x_range = float(input("xrange (cm): "))
	#z_range = float(input("zrange (cm): "))
	lESi=[[],[],[],[]]
	for h in range(len(ESi)):
		SESi = ESi[h].split()
		lESi[0].append(float(SESi[0]))  # energy
		lESi[1].append(float(SESi[2]))  # y
		lESi[2].append(float(SESi[3]))  # z
	
	h, xedges, yedges, image = plt.hist2d(lESi[1], lESi[2], BIN, density=True, weights=lESi[0], vmin = vmin, vmax = vmax, range=[[-x_range/2, x_range/2], [-6, 0]])
	plt.xlabel('y [cm]')
	plt.ylabel('z [cm]')
	plt.title(r'$\mathrm{Energy\ distribution\ in\ Silicon}$')
	plt.colorbar()
	plt.show()
	del SESi
	del lESi
