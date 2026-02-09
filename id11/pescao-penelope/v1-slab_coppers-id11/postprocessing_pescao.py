# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 15:28:23 2026

@author: BRUMUND
"""

import numpy as np
import matplotlib.pyplot as plt

# %% Plot: 2D histogram absorption in depth
data = np.loadtxt('pescao0-HistoSilicon-id11_21.3keV-perp_Si.dat', skiprows=2)

en_tot = 3.05e9
P_tot = 1893
x_range = 0.2
z_range = 6

h, xedges, yedges, image = plt.hist2d(data[:,1], data[:,3], 100, density=False, weights=data[:,0]/en_tot*P_tot, 
                                      vmin = 0, vmax = 10, 
                                      range=[[-x_range/2, x_range/2], [-z_range, 0]])
plt.xlabel('transverse coordinate x [cm]')
plt.ylabel('depth coordinate z [cm]')
plt.title(r'$\mathrm{Power\ distribution\ in\ Silicon}$ (100 bins)')
plt.colorbar(label='Power (W)')
plt.show()

# %% Plots: Scattering ID11 on detector
data = np.loadtxt('pescao0-HistoDetectorUp-id11_21.3keV-th5.33.dat', skiprows=3)
# data = np.loadtxt('test_beam-DepositedEnergyDetectorUp.dat', skiprows=1, usecols=[1,2,3,4])


E_tot = 3.06e9
P_tot = 1143.718
scale = P_tot/E_tot
# =============================================================================
# Y: longitudinal
# =============================================================================
BIN = 20

y_angle = np.rad2deg(np.tan(data[:,2]/0.1))

plt.subplots()
n, bins, patches = plt.hist(data[:,2]*10, BIN, 
                         density=False, weights=data[:,0]*scale,
                         facecolor='red', range=(-40,80))



plt.xlabel('longitudinal coordinate y [mm]')
plt.ylabel('P [W]')
plt.title(r'$\mathrm{Power\ distribution\ in\ the\ upper\ Detector (scattering)}$ (4mm from Si)')
plt.show()


# =============================================================================
# X: transversal
# =============================================================================
BIN=20
plt.subplots()
n, bins, patches = plt.hist(data[:,1]*10, BIN, 
                         density=False, weights=data[:,0]*scale,
                         facecolor='red', range=(-20,20))



plt.xlabel('transverse coordinate x [mm]')
plt.ylabel('P [W]')
plt.title(r'$\mathrm{Power\ distribution\ in\ the\ upper\ Detector}$ (4mm from Si)')
plt.show()

# %% Plots: Pescao 2 (modified geometry) 
data_si = np.loadtxt('pescao1-HistoSilicon-id11_21.3keV-th5.33_2.dat', skiprows=0)
# data_si = np.loadtxt('test_beam-DepositedEnergySilicon.dat', skiprows=1, usecols=[1,2,3,4])
data_db = np.loadtxt('pescao1-HistoDetectorDown-id11_21.3keV-th5.33_2.dat', skiprows=0)

E_tot = 6.10e+09
P_tot = 1143.718
scale = P_tot/E_tot

# %%% Plot 1: silicon y-z
# Plot y-z plane (2-3) - silicon
plt.subplots()
# si
y_range = 100
z_range = 6
z_res = 10
y_res = 1
h, xedges, yedges, image = plt.hist2d(data_si[:,2]*10, data_si[:,3]*10, 
                                      bins = [int(y_range*y_res),int(z_range*z_res)], 
                                      vmin = 0, vmax = 10, 
                                      density=False, 
                                      weights=data_si[:,0]*scale*y_res*z_res, 
                                      range=[[-0.2*y_range, 0.8*y_range], [-z_range, 0]])
# plt.gca().set_aspect('equal')
plt.xlabel('longitudinal coordinate y [mm]')
plt.ylabel('depth coordinate z [mm]')
plt.title(r'$\mathrm{Power\ distribution\ in\ Silicon}$ (y-z-plane, transverse)')
plt.colorbar(label='power density (W/mm2)')
plt.show()
# h, xedges, yedges, image = plt.hist2d(data_db[:,2]*10, data_db[:,3]*10, bins = [int(y_range),2], density=False, weights=data_db[:,0]/en_tot*P_tot, 
#                                       vmin = 0, vmax = 10, 
#                                       range=[[-0.2*y_range, 0.8*y_range], [-62, -60]])

# %%% Plot 2: plot bottom detector 1D y
plt.subplots()
y_range = 400
y_res = 0.1
n, bins, patches = plt.hist(data_db[:,2]*10, int(y_range*y_res), 
                         density=False, weights=data_db[:,0]*scale,
                         facecolor='red', range=(-0.1*y_range,0.9*y_range))



plt.xlabel('longitudinal coordinate y [mm]')
plt.ylabel('power [W]')
plt.title(r'$\mathrm{Power\ distribution\ in\ the\ lower \ Detector (scattering)}$ (touching Si bottom (-6cm))')
plt.show()

# %%% Plot 3: plot bottom detector 2D x-y
# Plot y-z plane (2-3) - silicon
plt.subplots()
# si
x_range = 2
y_range = 60
x_res = 25
y_res = 1
h, xedges, yedges, image = plt.hist2d(data_si[:,2]*10, data_si[:,1]*10, 
                                      bins = [int(y_range*y_res),int(x_range*x_res)], 
                                      vmin = 0, vmax = 10, 
                                      density=False, 
                                      weights=data_si[:,0]*scale*y_res*x_res, 
                                      range=[[-0.2*y_range, 0.8*y_range], [-0.5*x_range, 0.5*x_range]])
# plt.gca().set_aspect('equal')
plt.xlabel('longitudinal coordinate y [mm]')
plt.ylabel('transverse coordinate x [mm]')
plt.title(r'$\mathrm{Power\ distribution\ in\ lower Detector}$ (x-y-plane)')
plt.colorbar(label='power density (W/mm2)')
plt.show()
# h, xedges, yedges, image = plt.hist2d(data_db[:,2]*10, data_db[:,3]*10, bins = [int(y_range),2], density=False, weights=data_db[:,0]/en_tot*P_tot, 
#                                       vmin = 0, vmax = 10, 
#                                       range=[[-0.2*y_range, 0.8*y_range], [-62, -60]])