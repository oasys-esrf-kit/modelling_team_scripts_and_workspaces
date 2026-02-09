#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 20:42:29 2020

@author: philipp
"""

#
# script to make the calculations (created by XOPPY:undulator_spectrum)
#
from XpyBru.Sources import Undulator, ebs
from XpyBru.XoppyTools import Filter, calc_1d_spectrum
from XpyBru.BruFEMts import generate_1d_apdl_table_file

from dabax.dabax_xraylib import DabaxXraylib
import xraylib
import numpy as np
import matplotlib.pyplot as plt

cpmu20 = Undulator(0.0205, 2000/20.05, K_max=2.60)
cpmu18 = Undulator(0.018, 2000/18, K_max=1.87)

def calc_all_ids_from_e(energy = 10000, h=1.2, v=1.2, d=32, de = 40):
    e1, _, sp_pow1, _ = calc_1d_spectrum(ebs, cpmu18, 
                                         k=cpmu18.find_K_value_from_E(energy)[0], 
                                         h=h, v=v, d=d, 
                                         emin=3000, emax=150000, de=de)
    e1, _, sp_pow2, _ = calc_1d_spectrum(ebs, cpmu20, 
                                         k=cpmu20.find_K_value_from_E(energy)[0], 
                                         h=h, v=v, d=d, 
                                         emin=3000, emax=150000, de=de)
    e = e1
    sp_pow = sp_pow1+sp_pow2
    
    return e, sp_pow, sp_pow1, sp_pow2

def calc_all_ids_from_k(k18 = 2, k20 = 2, h = 1.2, v = 1.2, d = 32, de = 40):
    e1, _, sp_pow1, _ = calc_1d_spectrum(ebs, cpmu18, 
                                         k=k18, 
                                         h=h, v=v, d=d, 
                                         emin=3000, emax=150000, de=de)
    e1, _, sp_pow2, _ = calc_1d_spectrum(ebs, cpmu20, 
                                         k=k20, 
                                         h=h, v=v, d=d, 
                                         emin=3000, emax=150000, de=de)
    e = e1
    sp_pow = sp_pow1+sp_pow2
    
    return e, sp_pow, sp_pow1, sp_pow2

title = 'EBS - ID11 Refurbishment, d=30m, H=2mm, V=1mm, Si absorption'


#%% Plot: Compare Energy absorption & Normal Linear Attenuation

h, v = 2.47, 1.37
fe_filter = Filter(['C'], [3.51], thicks=[0.3])

# plot total power absorbed vs. depth P(z)
fig_P, ax_P=plt.subplots()      
title = f'EBS - ID11 Refurbishment, d=32m, H={h:.2}mm, V={v:.2f}mm, Si absorption'
ax_P.grid(True)
ax_P.set_xlabel('$ z \, \mathrm{(mm)}$')
ax_P.set_ylabel('$ dP_{abs}/dz \, \mathrm{(W/mm)}$')
ax_P.set_title(title, fontsize='small')

si_block = Filter(['Si'], ['?'], thicks=[100])

e1, sp_pow, _, _ = calc_all_ids_from_e(21300, h, v, de = 150)

sp_pow, _, _ = fe_filter.apply_filter(e1, sp_pow, material_constants_library=DabaxXraylib(file_f1f2="f1f2_Windt.dat",file_CrossSec="CrossSec_NIST_MassEnergyAbsorption.dat"))

x, P_x, p_v, f_v, P_n = si_block.calc_volumetric_abs(e1,sp_pow, A=h*v, discr=51,
                                                     material_constants_library=xraylib)
ax_P.plot(x, np.gradient(P_x,x), label = 'new prog xraylib')
x, P_x, p_v, f_v, P_n = si_block.calc_volumetric_abs(e1,sp_pow, A=h*v, discr=51,
                                                     material_constants_library=DabaxXraylib(file_f1f2="f1f2_Windt.dat",file_CrossSec="CrossSec_NIST_MassEnergyAbsorption.dat"))
ax_P.plot(x, np.gradient(P_x,x), label = 'new prog DabaxX Mass Energy Abs.', color = 'tab:green')


# load and add pescao (10cm=100mm, 100bins is one bin per mm)
data = np.loadtxt('HistoSilicon-id11_21.3keV-perp_Si.dat', skiprows=2)
en_tot = 3.03e9
P_tot = 1893
n, bins, patches = ax_P.hist(-1*data[:,3]*10, 100, density=False, weights=data[:,0]/en_tot*P_tot,facecolor='red', alpha=0.5, label = 'Pescao-Penelope (Monte-Carlo)')
# n, bins = np.histogram(-1*data[:,3]*10, 100, density=False, weights=data[:,0]/en_tot*P_tot)
# ax_P.plot(bins[:-1], n, label = 'Pescao-Penelope (Monte-Carlo)')


ax_P.legend()