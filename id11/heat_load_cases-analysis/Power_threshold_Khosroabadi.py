# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 09:39:00 2025

@author: BRUMUND
"""

import numpy as np 
from scipy.optimize import newton
import matplotlib.pyplot as plt
from XpyBru import Monochromators

# Universal function from Khosroabadi 2024 (eq. 8)
def f_univ(P, Pd = 20):
    Tze = 127
    
    Tcu = 80
    # k = 2000
    k = 3000
    A = 100*50*2*1e-6
    Tb = Tcu + P/k/A
    
    C = 7e-5
    
    f = C * Tb**2 * (P*Pd)**0.5 - 2 * (1 + C * Tze * (P*Pd)**0.5) * Tb + 2 * Tze
    
    return f


fig, ax = plt.subplots()
# pd = np.logspace(np.log10(5), np.log10(70), 20)
pd = np.logspace(np.log10(10), np.log10(50), 50)
Ps = np.zeros(pd.shape)

for i, pdi in enumerate(pd):
    P = newton(f_univ, 150, args=(pdi,), disp = False)
    print(pdi, P)
    Ps[i] = P
    # if np.isfinite(P): ax.scatter(pdi, P, s = 4**2, c='tab:blue')
    
ax.plot(pd, Ps, label = 'Threshold power (Khosroabadi 2024)')    

ax.set_xlabel('Power density (normal to crystal surface) [W/mm^2]')
ax.set_ylabel('Power [W]')
ax.set_xlim([0,pd.max()+10])
# ax.set_ylim([0,Ps.max()+50])

ax.grid()

# =============================================================================
# # Add ID11 points
# =============================================================================
# th = np.array)
data = np.loadtxt(r'c:\Users\brumund\Work Folders\Documents\Projects\ID11-Upgrade\Mono-Loadcases.csv', 
                  delimiter = ',', 
                  skiprows = 2,
                  usecols = (0,11,12,13))
Es, H_V, p, P = data[:,0], np.sqrt(data[:,1]), data[:,2], data[:,3]

# th = Monochromators.bragg_angle(15000)
# pd11 = np.sin(th)*1000
# P11 = 1000 * 0.8**2
# ax.scatter(pd11, P11, label = f'ID11, E=15keV (th={np.rad2deg(th):.1f}deg)  p=1000 W/mm^2')
colors = [[0.7, 'tab:blue'],
          [0.8, 'tab:orange'],
          [1.0, 'tab:green'],
          [1.4, 'tab:red']]
colors = np.array(colors)
def get_color(H):
    find_arr = np.float32(colors[:,0])
    ind = np.argmin(np.abs(find_arr-H))
    return colors[ind,1]
    

for Ei, Hi, pi, Pi in zip(Es, H_V, p, P):
    if Ei == 15:
        ax.scatter(pi, Pi, label = f'E={Ei} keV, H=V={Hi:.1f} mm', 
               marker = 'o', color = get_color(Hi))
    elif Ei == 18:
        ax.scatter(pi, Pi, label = f'E={Ei} keV, H=V={Hi:.1f} mm', 
               marker = 's', color = get_color(Hi))
    elif Ei == 24:
        ax.scatter(pi, Pi, label = f'E={Ei} keV, H=V={Hi:.1f} mm', 
               marker = 'D', color = get_color(Hi))
    elif Ei == 31.3:
        ax.scatter(pi, Pi, label = f'E={Ei} keV, H=V={Hi:.1f} mm', 
               marker = 'v', color = get_color(Hi))
    elif Ei == 40.3:
        ax.scatter(pi, Pi, label = f'E={Ei:.1f} keV, H=V={Hi:.1f} mm', 
               marker = '^', color = get_color(Hi))
    else:
        print(Ei, 'keV data not plotted')
        pass
    
ax.legend(fontsize = 'small')
ax.set_title('EBS - ID11 - CPMU18 (2m) and CPMU20.5 (2m) combined')

