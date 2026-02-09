import os
import copy
import sys
import math
import codecs
from string import *

#************************ Importing files **************************************

fromDir = sys.argv[0]
fromDir = fromDir[:-9]
#print ("fromDir: **"+fromDir+"**")
#print ("Working directory is: "+os.getcwd() )
if (fromDir != ""): 
    if (fromDir != os.getcwd() ) : 
        print ("Copying files from directory: "+fromDir)
        os.system("cp "+fromDir+"fortranCode/penEasy.x .")
        os.system("cp "+fromDir+"*.mat .")
        os.system("cp "+fromDir+"*.gpl .")
        os.system("cp "+fromDir+"*.tmp .")
        os.system("cp "+fromDir+"pescaoANALYSIS .")
        os.system("cp "+fromDir+"*.py .")
        os.system("cp "+fromDir+"Spectrum.dat .")
        os.system("cp "+fromDir+"README.txt .")
        #os.system("cp -r "+fromDir+"fortranCode .")

#******************** Removing files of the last run ****************************** 

try:
	os.remove("penEasy.out")
except OSError:
	pass	
try:
	os.remove("tallyEnergyDeposition.dat")
except OSError:
        pass
try:
	os.remove("tallyPulseHeightSpectrum.dat")
except OSError:
        pass
try:
	os.remove("SourceSpectrum.dat")
except OSError:
        pass
try:
	os.remove("HistoSilicon.dat")
except OSError:
        pass
try:
	os.remove("HistoDetectorUp.dat")
except OSError:
        pass
try:
	os.remove("HistoCopper.dat")
except OSError:
        pass
try:
	os.remove("HistoDetectorDown.dat")
except OSError:
        pass
try:
	os.remove("OutputValues.dat")
except OSError:
        pass
try:
	os.remove("SourceSpectrum.dat")
except OSError:
	pass


######################################################
#                                                    #
#                Input parameters                    #
#						     #
######################################################

num_particle = input("Number of particles [e.g. 10000]: ")
num_particle = int(num_particle)
if abs(math.fmod(num_particle,1))>1e-06:
	print ("Integer positive numbers!")
	sys.exit()
if num_particle <= 0:
   print ("Please put positive numbers!")
   num_particle = float(input("Number of particles [e.g. 10000]: "))
   if num_particle <= 0:
      print ("Sorry, I said that you have to put positive values!")
      sys.exit()	
DistanceDetUp = float(input("Distance between the upper detector plane and the silicon one (cm) [e.g. 10cm]: "))
if DistanceDetUp <= 0:
   print ("Please put positive numbers!")
   DistanceDetUp = float(input("Distance between the upper detector plane and the silicon one (cm) [e.g. 10cm]: "))
   if DistanceDetUp <= 0:
      print ("Sorry, I said that you have to put positive values!")
      sys.exit()
DistanceDetDown = float(input("Distance between the lower detector plane and the silicon one (cm) [e.g. 10cm]: "))
if DistanceDetDown <= 0:
   print ("Please put positive numbers!")
   DistanceDetDown = float(input("Distance between the lower detector plane and the silicon one (cm) [e.g. 10cm]: "))
   if DistanceDetDown <= 0:
      print ("Sorry, I said that you have to put positive values!")
      sys.exit()
DistanceCuSi = float(input("Distance between the copper plane and the silicon one (cm) [e.g. 5cm]: "))
if DistanceCuSi <= 0:
   print ("Please put positive numbers!")
   DistanceCuSi = float(input("Distance between the copper plane and the silicon one (cm) [e.g. 5cm]: "))
   if DistanceCuSi <= 0:
      print ("Sorry, I said that you have to put positive values!")
      sys.exit()
ThiknessSi = float(input("Thikness of the silicon plane (cm) [e.g. 10cm]: "))
if ThiknessSi <= 0:
   print ("Please put positive numbers!")
   ThiknessSi = float(input("Thikness of the silicon plane (cm) [e.g. 10cm]: "))
   if ThiknessSi <= 0:
      print ("Sorry, I said that you have to put positive values!")
      sys.exit()
ThiknessCu = float(input("Thikness of the copper plane (cm) [e.g. 0.5cm]: "))
if ThiknessCu <= 0:
   print ("Please put positive numbers!")
   ThiknessCu = float(input("Thikness of the silicon plane (cm) [e.g. 0.5cm]: "))
   if ThiknessCu <= 0:
      print ("Sorry, I said that you have to put positive values!")
      sys.exit()
theta = float(input("Theta [deg] (grazing) [e.g. 0.17 deg]: "))
phi = float(input("Angle of semiaperture of the source (0 -180)  [deg] [e.g. 0 ]: "))
source_position = eval(input("Location of the source in Cartesian coordinates (cm) [e.g. 0,0,0 ]: "))
Bin = input("Number of bin of the plots [e.g. 100]: ")
Bin=int(Bin)
if abs(math.fmod(Bin,1))>1e-06:
	print ("Integer positive numbers!")
	sys.exit()
if Bin <= 0:
   print ("Please put integer positive numbers!")
   Bin = int(input("Number of bin of the plots [e.g. 100]: "))
   if Bin <= 0:
      print ("Sorry, I said that you have to put integer positive values!")
      sys.exit()


#********************** Input parameters in the geometry ***********************************
	
# indata1 = open("phantom.geo.tmp","r")
# data1 = indata1.read()
# indata1.close()
# del indata1
# try:
# 	os.remove("phantom.geo")
# except OSError:
#         pass
# infile1 = open("phantom.geo","w")
# 
# DistanceSiDetDown = ThiknessSi + DistanceDetDown
# PlaneCuUp = ThiknessCu + DistanceCuSi
# 
# distanceSiDetUp = "%6.15e" % DistanceDetUp
# distanceSiDetDown = "%6.15e" % DistanceSiDetDown
# distanceCuSi = "%6.15e" % DistanceCuSi
# thiknessSi = "%6.15e" % ThiknessSi
# planeCuUp = "%6.15e" % PlaneCuUp
# 
# data1 = data1.replace("PPPPPPPPPPPPPPPPPPPPP",distanceSiDetUp.replace("e+","E+"))
# data1 = data1.replace("HHHHHHHHHHHHHHHHHHHHH",distanceSiDetDown.replace("e+","E+"))
# data1 = data1.replace("DDDDDDDDDDDDDDDDDDDDD",distanceCuSi.replace("e+","E+"))
# data1 = data1.replace("YYYYYYYYYYYYYYYYYYYYY",thiknessSi.replace("e+","E+"))
# data1 = data1.replace("CCCCCCCCCCCCCCCCCCCCC",planeCuUp.replace("e+","E+"))
# infile1.write(data1)
# infile1.close()
# del data1
# del infile1

num = "%6.2e" % num_particle
costheta = math.cos((theta*math.pi)/180)
sintheta = (math.sin((theta*math.pi)/180))*(-1)
Phi = "%.8f" % phi
x_SOURCE = "%.8f" % source_position[0]
y_SOURCE = "%.8f" % source_position[1]
z_SOURCE = "%.8f" % source_position[2]
DSmaxSi = ThiknessSi/10
DSMaxSi = "%6.3e" % DSmaxSi
DSmaxCu = ThiknessCu/10
DSMaxCu = "%6.3e" % DSmaxCu


passo = [0,3,6,9]
outfile = open("OutputValues.dat","w")
out=['Absorbed energy in the silicon = ','no','no','Scattered energy to the upper detector = ','no','no','Scattered energy to the lower detector = ','no','no','Scattered energy to the Copper = ']

choice = int(input("\nType of source? \n[1] Monochromatic \n[2] From file \n"))

######################################################
#                                                    #
#              Monochronatic Source                  #
#						     #
######################################################

if choice == 1:
	beamEnergy = float(input("Beam Energy (eV) [e.g. 10000eV]: "))
	if beamEnergy <= 0:
   		print ("Please put positive values of energy!")
   		beamEnergy = float(input("Beam Energy (eV) [e.g. 10000eV]: "))
   		if beamEnergy <= 0:
      			print ("Sorry, I said that you have to put positive values!")
      			sys.exit()
# ******************* Init writing of the outout file fill in with numerical results *************************
   		outfile.write("\n")
   		outfile.write('**************************************************'+"\n")
   		outfile.write('***********************************************'+"\n")
   		outfile.write("\n")
   		outfile.write(str('Total Incident Energy = ')+str(beamEnergy*num_particle)+" eV"+"\n")
   		outfile.write("\n")


# ******************* Substituting input parameters in input penEasy *********************************************
        	
	indata = open("penEasy.in.tmp","r")
	data = indata.read()
	indata.close()
	del indata
	try:
		os.remove("penEasy.in")
	except OSError:
		pass
	infile = open("penEasy.in","w")  #input penEasy
	         
	arg = "%6.2e" % beamEnergy
        	
	data = data.replace("XYXYX",arg.replace("e+0","e"))
	data = data.replace("ZZZZZ",num.replace("e+0","e"))
	data = data.replace("PP",str(Bin))
	data = data.replace("CCCCCCCC",str(costheta))
	data = data.replace("SSSSSSS",str(sintheta))
	data = data.replace("CONE",str(Phi))
	data = data.replace("POSX", str(x_SOURCE))
	data = data.replace("POSY", str(y_SOURCE))
	data = data.replace("POSZ", str(z_SOURCE))
	data = data.replace("DSMaxSi",DSMaxSi.replace("e+0","e"))
	data = data.replace("DSMaxCu",DSMaxCu.replace("e+0","e"))
	infile.write(data)
	infile.close()
	del data
	del infile


#************** Running penEasy ******************************************************************************

	print ("./penEasy.x < penEasy.in > penEasy.out")
	print ('...penEasy is running...wait please...')
	#os.system("script -c './penEasy.x < penEasy.in' penEasy.out")
	os.system("./penEasy.x < penEasy.in > penEasy.out")


#************** Calculating re-scattered and absorbed radiation ***********************************************
	
	incol = open("tallyEnergyDeposition.dat","r")
	col = incol.readlines()[3:]
	incol.close()
	del incol
	for c in passo:
		riga = col[c].split()
		col[c] = (float(riga[1])*num_particle, float(riga[2])*num_particle)
		outfile.write(str(out[c])+str(col[c][0])+"+-"+str(col[c][1])+" eV"+"\n")

	outfile.write("\n")
	PrSi = (col[0][0]/(beamEnergy*num_particle))*100
	ErrPrSi = (col[0][1]/(beamEnergy*num_particle))*100
	prSi ="%.3f" % PrSi 
	ErrprSi ="%.3f" % ErrPrSi
	outfile.write("Absorbed power fraction in the silicon: "+prSi+"+-"+ErrprSi+" %"+"\n")
	PrDetUp = (col[3][0]/(beamEnergy*num_particle))*100
	ErrPrDetUp = (col[3][1]/(beamEnergy*num_particle))*100
	prDetUp ="%.3f" % PrDetUp 
	ErrprDetUp ="%.3f" % ErrPrDetUp
	outfile.write("Absorbed power fraction in the upper detector: "+prDetUp+"+-"+ErrprDetUp+" %"+"\n")
	PrDetDown = (col[6][0]/(beamEnergy*num_particle))*100
	ErrPrDetDown = (col[6][1]/(beamEnergy*num_particle))*100
	prDetDown ="%.3f" % PrDetDown 
	ErrprDetDown ="%.3f" % ErrPrDetDown
	outfile.write("Absorbed power fraction in the lower detector: "+prDetDown+"+-"+ErrprDetDown+" %"+"\n")
	PrCu = (col[9][0]/(beamEnergy*num_particle))*100
	ErrPrCu = (col[9][1]/(beamEnergy*num_particle))*100
	prCu ="%.3f" % PrCu
	ErrprCu ="%.3f" % ErrPrCu
	outfile.write("Absorbed power fraction in the copper: "+prCu+"+-"+ErrprCu+" %"+"\n")
	outfile.write("\n")
	outfile.write('***********************************************'+"\n")
	outfile.write('**************************************************'+"\n") 
	del col
	del riga
	outfile.close()
	del outfile
	
	
#************** End of Monochromatic Section ********************************************************************


######################################################
#                                                    #
#                Spectrum Source                     #
#						     #
######################################################

elif choice == 2:

	fileSpectrum = input("Please enter the name of the spectrum file (.dat format): ")
	x = os.path.isfile(fileSpectrum)
	if x == False:
		print ("Sorry! This file doesn't exist!")
		sys.exit()
	try:
		os.remove("InputPenEasy.dat")
	except OSError:
		pass
	try:
		os.remove("output.in")
	except OSError:
		pass


#************** Creating input spectrum file for penEasy modifying user spectrum file ********************

	outfileS = open("InputPenEasy.dat","w")
        

	Incolo = open(fileSpectrum,"r")
	incolo = Incolo.readlines()
	Incolo.close()
	del Incolo
	colo = [[],[]]
	for c in range(len(incolo)):
		col = incolo[c].split()
		colo[0].append(float(col[0])) 
		colo[1].append(float(col[1])/(colo[0][c]/1000))

	step=(colo[0][1]-colo[0][0])/2

	for c in range(len(colo[0])):
    		colo[0][c] = float(colo[0][c]-step)
    		outfileS.write(str(colo[0][c])+"     "+str(colo[1][c])+"\n")
	outfileS.close()
	del outfileS

	step2 = colo[0][len(colo[0])-1]+2*step
	maxE = colo[0][len(colo[0])-1]+step

	del colo
	del incolo
	del col

	inprobtmp = open("penEasy.in.Spectrum.tmp","r")
	probtmp = inprobtmp.read()
	inprobtmp.close()
	del inprobtmp
	try:
		os.remove("penEasy.Spectrum.in")
	except OSError:
		pass


#********************** Substituting input parameters in penEasy *********************************************

	prob = open("penEasy.Spectrum.in","w")
	probtmp = probtmp.replace("EEEE",str(step2))
	probtmp = probtmp.replace("ZZZZZ",num.replace("e+0","e"))
	probtmp = probtmp.replace("YYYYY",str(maxE))	
	probtmp = probtmp.replace("CCCCCCCC",str(costheta))
	probtmp = probtmp.replace("SSSSSSS",str(sintheta))
	probtmp = probtmp.replace("CONE",str(Phi))
	probtmp = probtmp.replace("POSX", str(x_SOURCE))
	probtmp = probtmp.replace("POSY", str(y_SOURCE))
	probtmp = probtmp.replace("POSZ", str(z_SOURCE))
	probtmp = probtmp.replace("DSMaxSi",DSMaxSi.replace("e+0","e"))
	probtmp = probtmp.replace("DSMaxCu",DSMaxCu.replace("e+0","e"))
	probtmp = probtmp.replace("PP",str(Bin))
	prob.write(probtmp)
	prob.close()
	del probtmp
	del prob


#*********************** Creating a penEasy input file output.in ***********************************************************

	penel = open("penEasy.Spectrum.in","r")
	outS = open("output.in","w")  # penEasy input file

	line = penel.readline()
	while   line != 'Energy(eV)  Probabilty\n':
        	outS.write(line)
        	line = penel.readline()

	outS.write('Energy(eV)  Probabilty'+"\n")

	inlineu = open("InputPenEasy.dat","r")
	Lineu = inlineu.readlines()
	inlineu.close()
	del inlineu
	lineu = [[],[]]
	for c in range(len(Lineu)):
		lin = Lineu[c].split()
		lineu[0].append(lin[0])
		lineu[1].append(lin[1])
		outS.write(str(lineu[0][c])+" 	  "+str(lineu[1][c])+"\n")

	while line != '# >>>> END OF FILE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n':
		line = penel.readline()
		outS.write(line)
	
	outS.close()
	penel.close()
	del penel
	del outS
	del lineu
	del lin
	del Lineu


#************** Running penEasy ******************************************************************************

	print ('...penEasy is running...wait please...')
	#os.system("script -c './penEasy.x <output.in' penEasy.out")
	os.system("./penEasy.x <output.in > penEasy.out")


#************** Calculating re-scattered and absorbed radiation ***********************************************

	inenergySource = open("SourceSpectrum.dat","r")
	energySource = inenergySource.readlines()[0:]
	inenergySource.close()
	del inenergySource
	for e in range(len(energySource)):
		energySource[e] = float(energySource[e].split()[0])
	Etot=sum(energySource) 
	del energySource
	ETot = "%6.2e" % Etot 
	outfile.write("\n")
	outfile.write('**************************************************'+"\n") 
	outfile.write('***********************************************'+"\n") 
	outfile.write("\n")
	outfile.write("Total Incident Energy: "+ETot+" eV"+"\n")
	Emean=Etot/num_particle
	EMean = "%6.2e" % Emean 
	outfile.write("Mean energy of each particle: "+EMean+" eV"+"\n"+"\n")

	convfact1 = 1000.0*1.6021765e-19

	inPowerSource = open(fileSpectrum,"r")
	PowerSource = inPowerSource.readlines()[0:]
	inPowerSource.close()
	del inPowerSource
	for i in range(len(PowerSource)):
		PowerSource[i] = float(PowerSource[i].split()[1])*convfact1*2*step 
	PowerIn = sum(PowerSource)
	powerIn = "%.3f" % PowerIn
	del PowerSource 
	
	incol = open("tallyEnergyDeposition.dat","r")
	col = incol.readlines()[3:]
	incol.close()
	del incol
	for c in passo:
		riga = col[c].split()
		col[c] = (float(riga[1])*num_particle, float(riga[2])*num_particle)
		outfile.write(str(out[c])+str(col[c][0])+"+-"+str(col[c][1])+" eV"+"\n")
	outfile.write("\n")

	PowerSi = PowerIn*(col[0][0]/Etot)
	ErrPowerSi = PowerIn*(col[0][1]/Etot)
	Powersi = "%6.3e" % PowerSi
	ErrPowersi = "%6.3e" % ErrPowerSi
	PowerDetUp = PowerIn*(col[3][0]/Etot)
	ErrPowerDetUp = PowerIn*(col[3][1]/Etot)
	PowerdetUp = "%6.3e" % PowerDetUp
	ErrPowerdetUp = "%6.3e" % ErrPowerDetUp
	PowerDetDown = PowerIn*(col[6][0]/Etot)
	ErrPowerDetDown = PowerIn*(col[6][1]/Etot)
	PowerdetDown = "%6.3e" % PowerDetDown
	ErrPowerdetDown = "%6.3e" % ErrPowerDetDown
	PowerCu = PowerIn*(col[9][0]/Etot)
	ErrPowerCu = PowerIn*(col[9][1]/Etot)
	Powercu = "%6.3e" % PowerCu
	ErrPowercu = "%6.3e" % ErrPowerCu
	outfile.write("\n")
	outfile.write("Total Power in Source: "+powerIn+" W"+"\n")
	outfile.write("\n")
	outfile.write("Absorbed power in the silicon: "+Powersi+"+-"+ErrPowersi+" W"+"\n")
	outfile.write("Scattered power to the upper detector: "+PowerdetUp+"+-"+ErrPowerdetUp+" W"+"\n")
	outfile.write("Scattered power to the lower detector: "+PowerdetDown+"+-"+ErrPowerdetDown+" W"+"\n")
	outfile.write("Scattered power to the copper: "+Powercu+"+-"+ErrPowercu+" W"+"\n")
	outfile.write("\n")

	fractPowerSi = 100*col[0][0]/Etot
	ErrfractPowerSi = 100*col[0][1]/Etot
	fractPowerSil = "%.3f" % fractPowerSi
	ErrfractPowerSil = "%.3f" % ErrfractPowerSi
	fractPowerDetUp = 100*col[3][0]/Etot
	ErrfractPowerDetUp = 100*col[3][1]/Etot
	fractPowerDeteUp = "%.3f" % fractPowerDetUp
	ErrfractPowerDeteUp = "%.3f" % ErrfractPowerDetUp
	fractPowerDetDown =100*col[6][0]/Etot
	ErrfractPowerDetDown =100*col[6][1]/Etot
	fractPowerDeteDown = "%.3f" % fractPowerDetDown
	ErrfractPowerDeteDown = "%.3f" % ErrfractPowerDetDown
	fractPowerCu = 100*col[9][0]/Etot
	ErrfractPowerCu = 100*col[9][1]/Etot
	fractPowerCop = "%.3f" % fractPowerCu
	ErrfractPowerCop = "%.3f" % ErrfractPowerCu

	outfile.write("Absorbed energy fraction in the silicon: "+fractPowerSil+"+-"+ErrfractPowerSil+" %"+"\n")
	outfile.write("Absorbed energy fraction in the upper detector: "+fractPowerDeteUp+"+-"+ErrfractPowerDeteUp+" %"+"\n")
	outfile.write("Absorbed energy fraction in the lower detector: "+fractPowerDeteDown+"+-"+ErrfractPowerDeteDown+" %"+"\n")
	outfile.write("Absorbed energy fraction in the copper: "+fractPowerCop+"+-"+ErrfractPowerCop+" %"+"\n")
	outfile.write("\n")
	outfile.write('***********************************************'+"\n")
	outfile.write('**************************************************'+"\n") 
	outfile.close()
	del col
	del riga
	del outfile


#************** End of Spectrum Section ********************************************************************
	
else: 
	print ('Wrong choice!!!')
	sys.exit()


#****************** Check: if it is True make histograms*****************************************************
 	
intemps = open("penEasy.out","r")
temps = intemps.read()
intemps.close()
del intemps

t=0

if "Have a nice day" in temps:
	t=10


######################################################
#                                                    #
#                Creating histograms                 #
#						     #
######################################################

#******** Upper Detector *************************   
	outfileDetUp = open("HistoDetectorUp.dat","w")
	incolD = open("DepositedEnergyDetectorUp.dat","r")
	colD = incolD.readlines()[1:]
	incolD.close()
	del incolD
	for c in range(len(colD)):
		colD[c] = colD[c].split()[1:5]
		outfileDetUp.write(str(colD[c][0])+"  "+str(colD[c][1])+"  "+str(colD[c][2])+"  "+str(colD[c][3])+"\n")
	del colD
	outfileDetUp.close()
	del outfileDetUp
	# try:
		# os.remove("DepositedEnergyDetectorUp.dat")
        # pass
	# except OSError:
		# pass

#******** Silicon *************************
	outfileSi = open("HistoSilicon.dat","w")
	incolS = open("DepositedEnergySilicon.dat","r")
	colS = incolS.readlines()[1:]
	incolS.close()
	del incolS
	for c in range(len(colS)):
		colS[c] = colS[c].split()[1:5]
		outfileSi.write(str(colS[c][0])+"  "+str(colS[c][1])+"  "+str(colS[c][2])+"  "+str(colS[c][3])+"\n")
	outfileSi.close()
	del colS
	del outfileSi
	try:
		os.remove("DepositedEnergySilicon.dat") 
	except:
		pass

#******** Lower Detector *************************	
	outfileDetDown = open("HistoDetectorDown.dat","w")
	incolDD = open("DepositedEnergyDetectorDown.dat","r")
	colDD = incolDD.readlines()[1:]
	incolDD.close()
	del incolDD 
	for c in range(len(colDD)):
		colDD[c] = colDD[c].split()[1:5]
		outfileDetDown.write(str(colDD[c][0])+"  "+str(colDD[c][1])+"  "+str(colDD[c][2])+"  "+str(colDD[c][3])+"\n")
	outfileDetDown.close()
	del colDD  
	del outfileDetDown 
	try:
		os.remove("DepositedEnergyDetectorDown.dat")
	except OSError:
		pass

#******** Copper *************************	
	outfileCu = open("HistoCopper.dat","w")
	incolC = open("DepositedEnergyCopper.dat","r")
	colC = incolC.readlines()[1:]
	incolC.close()
	del incolC
	for c in range(len(colC)):
		colC[c] = colC[c].split()[1:5]
		outfileCu.write(str(colC[c][0])+"  "+str(colC[c][1])+"  "+str(colC[c][2])+"  "+str(colC[c][3])+"\n")
	outfileCu.close()
	del colC
	del outfileCu
	try:
		os.remove("DepositedEnergyCopper.dat") 
	except OSError:
		pass
	
	outfilePLOT = open("EnergyHisto.dat","w")
	inplot = open("tallyPulseHeightSpectrum.dat","r")
	plot = inplot.readlines()
	inplot.close()
	del inplot
	pl = ""
	for tmp in range(len(plot)):
		if plot[tmp].find("#C Have a nice day.")==0 : 
			pl = tmp-7
	incol = open("tallyPulseHeightSpectrum.dat","r")
	col = incol.readlines()
	incol.close()
	del incol
	for c in range(9,pl):
		col[c] = col[c].split()[0:13]
		outfilePLOT.write(str(col[c][0])+"  "+str(col[c][1])+"  "+str(col[c][2])+"  "+str(col[c][3])+"  "+str(col[c][4])+"  "+str(col[c][5])+"  "
	+str(col[c][6])+"  "+str(col[c][7])+"  "+str(col[c][8])+"  "+str(col[c][9])+"  "+str(col[c][10])+"  "+str(col[c][11])+"\n")
	outfilePLOT.write(str(col[pl-1][1])+"  "+str(col[pl-1][1])+"  "+str(col[pl-1][2])+"  "+str(col[pl-1][3])+"  "+str(col[pl-1][4])+"  "
	+str(col[pl-1][5])+"  "+str(col[pl-1][6])+"  "+str(col[pl-1][7])+"  "+str(col[pl-1][8])+"  "+str(col[pl-1][9])+"  "+str(col[pl-1][10])
	+"  "+str(col[pl-1][11])+"\n")		
	del col
	outfilePLOT.close()

	
	print ("")
	print ("All is right!")


	os.system("./pescaoANALYSIS")


#************************** if it is False exit ***********************************************************    

elif t == 0:
	print ("")
	print ("Sorry! penEasy didn't run in a good way...let change the input parameters!")
	print ("")

del temps
#try:
#	os.remove("penEasy.out")
#except OSError:
#	pass	

#************************ End *******************************************************************************
