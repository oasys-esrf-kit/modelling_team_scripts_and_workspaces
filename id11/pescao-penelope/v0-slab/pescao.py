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
Thikness = float(input("Thickness of the silicon plane (cm) [e.g. 10cm]: "))
if Thikness <= 0:
   print ("Please put positive numbers!")
   Thikness = float(input("Thickness of the silicon plane (cm) [e.g. 10cm]: "))
   if Thikness <= 0:
      print ("Sorry, I said that you have to put positive values!")
      sys.exit()
Distance = float(input("Distance between the detector plane and the silicon one (cm) [e.g. 10cm]: "))
if Distance <= 0:
   print ("Please put positive numbers!")
   Distance = float(input("Distance between the detector plane and the silicon one (cm) [e.g. 10cm]: "))
   if Distance <= 0:
      print ("Sorry, I said that you have to put positive values!")
      sys.exit()
theta = float(input("Theta [deg] (grazing) [0.17 deg]: "))
Bin = int(input("Number of bins of the plots [e.g. 100]: "))
if abs(math.fmod(Bin,1))>1e-06:
	print ("Integer positive numbers!")
	sys.exit()
if Bin <= 0:
   print ("Please put integer positive numbers!")
   Bin = int(input("Number of bins of the plots [e.g. 100]: "))
   if Bin <= 0:
      print ("Sorry, I said that you have to put integer positive values!")
      sys.exit()


#********************** Input parameters in the geometry ***********************************            
	
indata1 = open("phantom.geo.tmp","r")
data1 = indata1.read()
indata1.close()
del indata1

try:
	os.remove("phantom.geo")
except OSError:
	pass
infile1 = open("phantom.geo","w")

thi = "%6.15e" % Thikness
dis = "%6.15e" % Distance
data1 = data1.replace("YYYYYYYYYYYYYYYYYYYYY",thi.replace("e+","E+"))
data1 = data1.replace("PPPPPPPPPPPPPPPPPPPPP",dis.replace("e+","E+"))
infile1.write(data1)
infile1.close()
del data1
del infile1



num = "%6.2e" % num_particle
costheta = math.cos((theta*math.pi)/180)
sintheta = (math.sin((theta*math.pi)/180))*(-1)
DSmax = Thikness/10
DSmaxSi = "%6.3e" % DSmax

passo = [0,3]
outfile = open("OutputValues.dat","w")
out=['Absorbed energy in the silicon = ','no','no','Scattered energy to the detector = ']


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


# ******************* Init writing of the outout file fill in with numerical results *********************************************  

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
	infile = open("penEasy.in","w")  
        
	arg = "%6.2e" % beamEnergy
        	
	data = data.replace("XYXYX",arg.replace("e+0","e"))
	data = data.replace("ZZZZZ",num.replace("e+0","e"))
	data = data.replace("dSmax",DSmaxSi.replace("e+0","e"))
	data = data.replace("PP",str(Bin))
	data = data.replace("CCCCCCCC",str(costheta))
	data = data.replace("SSSSSSS",str(sintheta))
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
	PrSiErr = (col[0][1]/(beamEnergy*num_particle))*100
	prSi ="%.2f" % PrSi
	prSiErr = "%.2f" % PrSiErr
	outfile.write("Absorbed power fraction in the silicon: "+str(prSi)+"+-"+str(prSiErr)+" %"+"\n")
	PrDet = (col[3][0]/(beamEnergy*num_particle))*100
	PrDetErr = (col[3][1]/(beamEnergy*num_particle))*100
	prDet ="%.2f" % PrDet
	prDetErr = "%.2f" % PrDetErr  
	outfile.write("Absorbed power fraction in the detector: "+str(prDet)+"+-"+str(prDetErr)+" %"+"\n")
	outfile.write("\n")
	outfile.write('***********************************************'+"\n")
	outfile.write('**************************************************'+"\n")
	outfile.close()
	del outfile
	del col
	del riga
	
	
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
	except:
		pass
	try:
		os.remove("output.in")
	except:
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
	maxE = colo[0][len(colo[0])-1]+2*step

	del colo
	del incolo
	del col


#********************** Substituting input parameters in penEasy *********************************************
	
	inprobtmp = open("penEasy.in.Spectrum.tmp","r")
	probtmp = inprobtmp.read()
	inprobtmp.close()
	del inprobtmp

	try:
		os.remove("penEasy.Spectrum.in")
	except OSError:
		pass
	prob = open("penEasy.Spectrum.in","w")
	probtmp = probtmp.replace("EEEE",str(step2))
	probtmp = probtmp.replace("ZZZZZ",num.replace("e+0","e"))
	probtmp = probtmp.replace("dSmax",DSmaxSi.replace("e+0","e"))
	probtmp = probtmp.replace("YYYYY",str(maxE))	
	probtmp = probtmp.replace("CCCCCCCC",str(costheta))
	probtmp = probtmp.replace("SSSSSSS",str(sintheta))
	probtmp = probtmp.replace("PP",str(Bin))
	prob.write(probtmp)
	prob.close()
	del probtmp
	del prob


#*********************** Creating a penEasy input file output.in ***********************************************************

	penel = open("penEasy.Spectrum.in","r")
	outS = open("output.in","w")  # penEasy input file

	line = penel.readline()
	while line != 'Energy(eV)  Probabilty\n':
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


	while   line != '# >>>> END OF FILE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n':
		line = penel.readline()
		outS.write(line)
	
	outS.close()
	penel.close()
	del outS
	del penel
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
	outfile.write("Total Incident Energy: "+str(ETot)+" eV"+"\n")
	Emean=Etot/num_particle
	EMean = "%6.2e" % Emean
	outfile.write("Mean energy of each particle: "+str(EMean)+" eV"+"\n"+"\n")
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
	outfile.write("\n")
	outfile.write("Total Power in Source: "+str(powerIn)+" W"+"\n")
	outfile.write("\n")

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
	outfile.write("Absorbed power in the silicon: "+str(Powersi)+"+-"+str(ErrPowersi)+" W"+"\n")
	PowerDet = PowerIn*(col[3][0]/Etot)
	ErrPowerDet = PowerIn*(col[3][1]/Etot)
	Powerdet = "%6.3e" % PowerDet
	ErrPowerdet = "%6.3e" % ErrPowerDet
	outfile.write("Scattered power to the detector: "+str(Powerdet)+"+-"+str(ErrPowerdet)+" W"+"\n")
	outfile.write("\n")
	fractPowerSi = 100*col[0][0]/Etot
	ErrfractPowerSi = 100*col[0][1]/Etot
	fractPowerSil = "%.3f" % fractPowerSi
	ErrfractPowerSil = "%.3f" % ErrfractPowerSi
	outfile.write("Absorbed energy fraction in the silicon: "+str(fractPowerSil)+"+-"+str(ErrfractPowerSil)+" %"+"\n")
	fractPowerDet = 100*col[3][0]/Etot
	ErrfractPowerDet = 100*col[3][1]/Etot
	fractPowerDete = "%.3f" % fractPowerDet
	ErrfractPowerDete = "%.3f" % ErrfractPowerDet
	outfile.write("Absorbed energy fraction in the detector: "+str(fractPowerDete)+"+-"+str(ErrfractPowerDete)+" %"+"\n")
	outfile.write("\n")
	outfile.write('***********************************************'+"\n")
	outfile.write('**************************************************'+"\n")
	outfile.close()
	del outfile
	del col
	del riga


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
	try:
		os.remove("DepositedEnergyDetectorUp.dat")
	except:
		pass

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
		col[c] = col[c].split()[0:9]
		outfilePLOT.write(str(col[c][0])+"  "+str(col[c][1])+"  "+str(col[c][2])+"  "+str(col[c][3])+"  "+str(col[c][4])+"  "+str(col[c][5])+"  "
		+str(col[c][6])+"  "+str(col[c][7])+"\n")
	outfilePLOT.write(str(col[pl-1][1])+"  "+str(col[pl-1][1])+"  "+str(col[pl-1][2])+"  "+str(col[pl-1][3])+"  "+str(col[pl-1][4])+"  "
	+str(col[pl-1][5])+"  "+str(col[pl-1][6])+"  "+str(col[pl-1][7])+"\n")		
	del col
	outfilePLOT.close()
	del outfilePLOT


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
#        pass

#************************ End ******************************************************************************
