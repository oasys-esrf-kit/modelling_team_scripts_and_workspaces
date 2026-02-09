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

inEDD = open("HistoDetectorDown.dat","r")
EDD = inEDD.readlines()[0:]
inEDD.close()
del inEDD
if len(EDD)==0:
	print ("There is no transmitted radiation!!!")
else:
	BIN = int(input("Number of bins in the histogram: "))
	lEDD = [[],[],[],[]]
	for h in range(len(EDD)):
		SEDD = EDD [h].split()
		lEDD[0].append(float(SEDD[0]))
		lEDD[1].append(float(SEDD[2]))
	n, bins, patches = plt.hist(lEDD[1], BIN, density=True, weights=lEDD[0],facecolor='red')
	plt.xlabel('y [cm]')
	plt.ylabel('E [eV]')
	plt.title(r'$\mathrm{Energy\ distribution\ in\ the\ lower\ Detector}$')
	plt.show()
	del lEDD
	del SEDD
