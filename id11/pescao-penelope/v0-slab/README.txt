				pescao - version 1
                                  by Eleonora Secco
  				        ESRF

1) PURPOSE OF THIS TOOL
2) DESCRIPTION OF THE FILES
3) HOW TO RECOMPILE
4) HOW TO USE IT
5) EXAMPLE
6) CONTACT

***********************************************************************************************************************************************************

1) PURPOSE OF THIS TOOL

PESCAO is a computer tool build for assesting the amount of absorbed and back-scattered power from the optical elements. The geometry consists of a silicon slab (optical element) and of a eminfinite plane of water (back-scattered particles counters). For attesting those quantities a full simulation package, named PENELOPE, and its main program penEasy are used. PENELOPE simulates the transport of electrons, positrons and photons in different materials and in complex geometries.
PESCAO is a computer tools written in Python that permits to run PENELOPE and penEasy automatically, to extract and displaying the results.

***********************************************************************************************************************************************************

2) DESCRIPTION OF THE FILES

- ./documentation/changeHistory.txt
  List of changes with respect to previous versions.

- ./documentation/sourceXXX.txt
  ./documentation/tallyXXX.txt
  Text files containing a description of the input data for each source and tally.

- ./fortranCode/penaux.f
  Auxiliary routines that initialize penEasy/PENELOPE.

- ./fortranCode/penEasy.f
  The main program source code.

- ./fortranCode/penpatch.f
  Subroutines that supersede their counterparts in the PENELOPE package. The physics in these routines is identical to the original, but provide more 
  information on the outcome of some interaction processes.

- ./fortranCode/penvox.f
  Geometry package for the simulation of voxelized geometries.

- ./fortranCode/penvr.f
  Implementation of various simple variance reduction techniques.

- ./fortranCode/sourceXXX.f
  Each source model is coded in a different file, e.g. sourcePhaseSpaceFile.f. They are described in the accompanying ./documentation/sourceXXX.txt files.

- ./fortranCode/tallyXXX.f
  Each tally is coded in a different file, e.g. tallySpatialDoseDistrib.f. They are described in the accompanying ./documentation/tallyXXX.txt files.

- ./fortranCode/timing.f

- ./XXXX.mat
  These are the material file of PENELOPE
- ./phantom.geo.tmp
  This is the file of the considered geometry (simple quadric geometry)

- ./penEasy.x
  Executable file that permits to run penEasy

- ./penEasy.in.XXXX.tmp
  Files di imput di penEasy that contain the source models, the fundamentals parameter for the simulation, tallies, ...

- ./tallyXXXX.gpl and HistoXXXX.py
  Files for plotting data

- ./penEasyRUN.py
  Computer tool for the assestment of transmitted, absorbed, irradiated radiation on the optics

- ./penEasyANALYSIS
  File bash for the results

***********************************************************************************************************************************************************


3) HOW TO RECOMPILE
 
1) You enter in fortranCode directory: cd fortranCode
2) You make you changes to the code
3) You compile: gfortran penEasy.f -o penEasy.x
4) You copy penEasy.x in your main directory: cp penEasy.x ..

***********************************************************************************************************************************************************


4) HOW TO USE IT
 
NOTE: Hereafter is assumed that the installation is carried out on a GNU/Linux system and that Python and its library "matplotlib" are installed in your system.

1) Create a new directory (e.g., mkdir ~/pescao, cd ~/pescao) in your working area where all the files for your job will be stored. Type the command:

python <installation_path>/pescao.py  

(e.g., python ~/GIT/pescao/pescaoVersion1/pescao.py )

and all the files that you need for running PESCAO will be copied automatically in your new directory and the program start to run. 

For further runnings you just type the command "python pescao.py".


2) After running the program, if you want to modify ~/tallyXXXX.gpl or HistoXXXX.py for improving the plots, you can show again the result, without
   running again the program, typing ~/penEasyANALYSIS.


***********************************************************************************************************************************************************
5) EXAMPLE - MONOCHROMATIC SOURCE or SPECTRUM SOURCE

Follow these steps:

a) Create a working directory and cd to it

b) Run the program:
   python <installation directoru>/pescao.py

c) You will find the following questions: 

- Number of particles [e.g. 10000] 
- Thickness of the silicon plane (cm) [e.g. 10cm]
- Distance between the detector plane and the silicon one (cm) [e.g. 10cm]
- Theta [deg] (grazing) [0.17 deg]
- Number of bins of the plots [e.g. 100]
- Type of source? [1] Monochromatic [2] From file
In the case of the monochromatic source you enter the energy beam in eV and in the case of the spectrum you have to enter a "file.dat" of two columns (energy in eV and the flux in ph/s/0.1%bw; there is an example of spectrum source in the file "Spectrum.dat" (999 ponits maximum))

d) After enetring your simulation parameters, you will obtain the results:

NUMERICAL RESULTS
- beam total incident energy [eV]
- absorbed energy in the silicon [eV]
- scattered energy to the detector [eV]
- absorbed power fraction in the silicon [%]
- absorbed power fraction in the detector [%]

GRAPHICAL RESULTS
- particle energy distribution in the silicon
- particle energy distribution in the detector
- superimposition of source, silicon and detector particle energy distributions
- silicon particle energy distribution vs (x,y,z)
- detector particle energy distribution vs (x,y,z)

*************************************************************************************************************************************************************
6) CONTACT

You can send your comments or questions to srio@esrf.eu
