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

inED = open("HistoDetectorUp.dat","r")
ED = inED.readlines()[0:]
inED.close()
del inED
if len(ED)==0:
	print ("There is no transmitted radiation!!!")
else:
	BIN = int(input("Number of bins in the histogram: "))
	lED = [[],[],[],[]]
	for h in range(len(ED)):
		SED = ED [h].split()
		lED[0].append(float(SED[0]))
		lED[1].append(float(SED[2]))
	n, bins, patches = plt.hist(lED[1], BIN, density=True, weights=lED[0],facecolor='red')
	plt.xlabel('y [cm]')
	plt.ylabel('E [eV]')
	plt.title(r'$\mathrm{Energy\ distribution\ in\ the\ upper\ Detector}$')
	plt.show()
	del lED
	del SED
