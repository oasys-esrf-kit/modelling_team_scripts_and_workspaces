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


inECu = open("HistoCopper.dat","r")
ECu = inECu.readlines()[0:]
inECu.close()
del inECu
if len(ECu)==0:
	print ("There is no absorbtion in the copper")
else:
	BIN = int(input("Number of bins in the histogram: "))
	lECu=[[],[],[],[]]
	for h in range(len(ECu)):
		SECu = ECu[h].split()
		lECu[0].append(float(SECu[0]))
		lECu[1].append(float(SECu[1]))
	n, bins, patches = plt.hist(lECu[1], BIN, density=True, weights=lECu[0],facecolor='red')
	plt.xlabel('x [cm]')
	plt.ylabel('E [eV]')
	plt.title(r'$\mathrm{Energy\ distribution\ in\ Copper}$')
	plt.show()
	del SECu
	del lECu
