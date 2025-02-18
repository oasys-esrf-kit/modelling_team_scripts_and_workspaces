# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 18:25:02 2025

@author: BRUMUND
"""


import numpy
# from crystalpy.util.calc_xcrystal import calc_xcrystal_angular_scan, calc_xcrystal_energy_scan, calc_xcrystal_alphazachariasen_scan
from crystalpy.util.calc_xcrystal import calc_xcrystal_angular_scan
import matplotlib.pyplot as plt
from oasys.util.oasys_util import get_fwhm
from scipy.signal import find_peaks, peak_widths




fig, ax = plt.subplots()
lines = []
leg = []
fwhms = []

# =============================================================================
# Comparisons
# =============================================================================
es = [4200, 6000, 7000, 8000, 10886, 11218, 15000, 20000, 20000, 20000]
hkls = [[3,1,1],
       [3,3,3],
       [4,4,0],
       [4,4,4],
       [6,6,4],
       [8,4,4],
       [4,4,4],
       [3,1,1],
       [4,4,0],
       [4,4,4]]
powers = [1,1,1,1,1,1,1,1,1,1]
emin = 0e-6
emax = 150e-6
n = 3000
alphas = [0, 
          5, 
          4, 
          5, 
          0, 
          0, 
          5, 
          0, 
          4, 
          5]

# es = [4200, 6000, 7800, 15000, 19000]
# hkls = [[3,1,1],
#        [3,3,3],
#        [3,3,3],
#        [8,4,4],
#        [8,4,4]]
# powers = [1,1,1,1,1]
# emin = -5e-6
# emax = 150e-6
# n = 4000

# es = [4200, 6000, 7800, 15000, 19000]
# hkls = [[1,1,1],
#        [1,1,1],
#        [1,1,1],
#        [1,1,1],
#        [1,1,1]]
# powers = [1,1,1,1,1]
# emin = -5e-6
# emax = 150e-6
# n = 1000

# High power comparison
# es = [19000,19000, 19000]
# hkls = [[8,4,4],
#         [4,4,0],
#         [4,4,4]]
# # powers = [2, 4, 4]
# powers = [1,1,1]
# emin = -5e-6
# emax = 15e-6
# n = 1000

# =============================================================================
# Script
# =============================================================================
   
for e, hkl, power, al in zip(es, hkls, powers, alphas): 
    bunch_out_dict, diffraction_setup, deviations = calc_xcrystal_angular_scan(
        # material_constants_library_flag=self.material_constants_library_flag,
        crystal_name              = 'Si',
        thickness                 = 0.007,
        miller_h                  = hkl[0],
        miller_k                  = hkl[1],
        miller_l                  = hkl[2],
        asymmetry_angle           = numpy.deg2rad(al),
        energy                    = e,
        angle_deviation_min       = emin,
        angle_deviation_max       = emax,
        angle_deviation_points    = n,
        angle_center_flag         = 2,
        calculation_method        = 1, # 0=Zachariasen, 1=Guigay
        is_thick                  = 0,
        use_transfer_matrix       = 0,
        geometry_type_index       = 0,
        calculation_strategy_flag = 0, # 0=mpmath 1=numpy 2=numpy-truncated
                )
    
    tmp = numpy.zeros((bunch_out_dict["energies"].size,7))
    tmp[:, 0] = deviations / 1e-06
    tmp[:, 1] = 8000.0
    tmp[:, 2] = bunch_out_dict["phaseP"]
    tmp[:, 3] = bunch_out_dict["phaseS"]
    # tmp[:, 4] = circular polarization
    tmp[:, 5] = bunch_out_dict["intensityP"]
    tmp[:, 6] = bunch_out_dict["intensityS"]**power

    fwhm, quote, coordinates = get_fwhm(tmp[:,6], tmp[:,0], ret0=None)
    fwhms.append(fwhm)
    
    ln = ax.plot(tmp[:,0], tmp[:,6]*100)
    lines.append(ln[0])
    # leg.append(f'{e/1000:.1f} keV, Si({hkl[0]}{hkl[1]}{hkl[2]})')
    th = numpy.rad2deg(diffraction_setup.angleBragg(e))
    # leg.append(f'{e/1000:4.1f} keV, Si({hkl[0]}{hkl[1]}{hkl[2]}), fwhm={fwhm:.2f} urad')
    
    leg.append(f'{e/1000:4.1f} keV, Si({hkl[0]}{hkl[1]}{hkl[2]}), th={th:.1f} deg, fwhm={fwhm:.2f} urad')
    
    
for i in range(len(leg)):
    print(leg[i])
    
# from srxraylib.plot.gol import plot
# plot(tmp[:,0], tmp[:,6], tmp[:,0], tmp[:,5], xtitle="angle", legend=["S-pol","P-pol"])
# plot(tmp[:,0], tmp[:,6] xtitle="angle", legend=["S-pol","P-pol"])
ax.grid()
ax.legend(lines, leg)
ax.set_xlabel('angle $\Theta - \Theta_B \, (\mu rad)$')
ax.set_ylabel('reflectivity (%)')
ax.set_ylim([0,100])
# print(fwhm)