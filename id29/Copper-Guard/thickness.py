# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 17:24:01 2025

@author: BRUMUND
"""
import numpy as np
from srxraylib.plot.gol import plot
import contextlib
import matplotlib.pyplot as plt

def spectral_calc(H_off, gapH = 10e-3, gapV = 1, emax = 1.5e5, n = 2000, und = 'cpmu'):
    #
    # script to make the calculations (created by XOPPY:undulator_spectrum)
    #
    from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_spectrum
    with contextlib.redirect_stdout(None):
        if und == 'cpmu':
            energy, flux, spectral_power, cumulated_power = xoppy_calc_undulator_spectrum(
                ELECTRONENERGY=6.0,
                ELECTRONENERGYSPREAD=0.001,
                ELECTRONCURRENT=0.2,
                ELECTRONBEAMSIZEH=3.01836e-05,
                ELECTRONBEAMSIZEV=3.63641e-06,
                ELECTRONBEAMDIVERGENCEH=4.36821e-06,
                ELECTRONBEAMDIVERGENCEV=1.37498e-06,
                PERIODID=0.0164,
                NPERIODS=121,
                KV=1.325,
                KH=0.0,
                KPHASE=0.0,
                DISTANCE=31.5,
                GAPH=gapH * 1e-3,
                GAPV=gapV * 1e-3,
                GAPH_CENTER=H_off * 1e-3,
                GAPV_CENTER=0.0,
                PHOTONENERGYMIN=8500.0,
                PHOTONENERGYMAX=emax,
                PHOTONENERGYPOINTS=n,
                METHOD=2,
                USEEMITTANCES=1)
        elif und == 'ivu':
            energy, flux, spectral_power, cumulated_power = xoppy_calc_undulator_spectrum(
                ELECTRONENERGY=6.0,
                ELECTRONENERGYSPREAD=0.001,
                ELECTRONCURRENT=0.2,
                ELECTRONBEAMSIZEH=3.01836e-05,
                ELECTRONBEAMSIZEV=3.63641e-06,
                ELECTRONBEAMDIVERGENCEH=4.36821e-06,
                ELECTRONBEAMDIVERGENCEV=1.37498e-06,
                PERIODID=0.021,
                NPERIODS=95.238,
                KV=0.966,
                KH=0.0,
                KPHASE=0.0,
                DISTANCE=31.5,
                GAPH=gapH * 1e-3,
                GAPV=gapV * 1e-3,
                GAPH_CENTER=H_off * 1e-3,
                GAPV_CENTER=0.0,
                PHOTONENERGYMIN=8500.0,
                PHOTONENERGYMAX=emax,
                PHOTONENERGYPOINTS=n,
                METHOD=2,
                USEEMITTANCES=1)
    #
    # # example plot
    # #
    # if True:
    #     from srxraylib.plot.gol import plot
        
    #     plot(energy,flux,
    #         xtitle="Photon energy [eV]",ytitle="Flux [photons/s/o.1%bw]",title="Undulator Flux",
    #         xlog=False,ylog=False,show=False)
    #     plot(energy,spectral_power,
    #         xtitle="Photon energy [eV]",ytitle="Power [W/eV]",title="Undulator Spectral Power",
    #         xlog=False,ylog=False,show=False)
    #     plot(energy,cumulated_power,
    #         xtitle="Photon energy [eV]",ytitle="Cumulated Power [W]",title="Undulator Cumulated Power",
    #         xlog=False,ylog=False,show=True)
    # #
    # # end script
    # #

    return energy, spectral_power

def transmission_calc(H_off, t_cu = 0.1, gapH = 10e-3, gapV = 1, und = 'cpmu'):
    
    energy, spectral_power = spectral_calc(H_off, und = und)
    
    from xoppylib.power.xoppy_calc_power import xoppy_calc_power
    import xraylib
    
    if t_cu == 0: t_cu = 1e-9
    
    out_dictionary = xoppy_calc_power(
            energy,
            spectral_power,
            substance = ['C','Cu',],
            thick     = [0.3, t_cu,], # in mm (for filters)
            angle     = [3, 3,], # in mrad (for mirrors)
            dens      = [3.51, '?',],
            roughness = [0, 0,], # in A (for mirrors)
            flags     = [0, 0,], # 0=Filter, 1=Mirror
            nelements = 2,
            FILE_DUMP = 0,
            material_constants_library = xraylib,
            )
    
    
    # data to pass
    energy = out_dictionary["data"][0,:]
    spectral_power = out_dictionary["data"][-1,:]
    
    # Transmitted power through filters
    P_trans = np.trapz(spectral_power, energy)
    
    # Transmitted power density filters
    p_trans = P_trans / gapH / gapV
    
    
    return p_trans, P_trans
    
def transm_from_x_tcu(H_off, t_cu, und = 'cpmu'):
    P, p = [], []
    i = 0
    n = H_off.size
    
    for xi, ti in zip(H_off, t_cu):
        print(f'Running point {i:3.0f}/{n:3.0f}... xi={xi:5.3f}mm, ti={ti:5.3f}mm')
        p_trans, P_trans = transmission_calc(xi, ti, und = und)
        P.append(P_trans)
        p.append(p_trans)
        print(f'Calculated powers: p={p_trans:6.1f}W/mm^2, P={P_trans:6.1f}W')
        i += 1
    return np.array(P), np.array(p)

def d1(x, th = 1.12):
    t = np.deg2rad(th)
    return x*np.tan(t)

def d2(x, th = 1.12, L2 = 3):
    t = np.deg2rad(th)
    xlim = np.sin(t)*L2
    d2 = np.zeros(x.shape)
    d2 = np.where(x <= xlim, x/np.tan(t), L2*np.cos(t))
    return d2

def d3(x, th = 1.12, al=12, L2 = 3):
    t = np.deg2rad(th)
    a = np.deg2rad(al)
    xlim = np.sin(t)*L2
    d3 = np.zeros(x.shape)
    d3 = np.where(x <= xlim, 0, (x-xlim)/np.tan(t+a))
    return d3

def d(x, th = 1.12, al=12, L2 = 3):
    return d1(x,th) + d2(x, th, L2) + d3(x, th, al, L2)



if __name__ == '__main__':
    plt.close('all')
    
    # =============================================================================
    # Create geometry
    # =============================================================================
    th = 1.12
    t = np.deg2rad(th)
    x0 = np.sin(t)*60.5
    
    # Create x array up to edge of front guard and further
    x1 = np.linspace(0, np.sin(t)*3, 41)
    x2 = np.linspace(np.sin(t)*3, 1, 21)
    x2 = np.delete(x2, 0)
    x = np.insert(x2, 0, x1)
    

    dtot = d(x)
    
    # Transform x-coordinate (as rotation center on ML1 center)
    x_new = x + x0
    
    
    # Insert values for chamfer outside of front guard
    x_ins = np.linspace(np.sin(t)*45, np.sin(t)*60.5, 12)
    x_ins = np.delete(x_ins, -1)
    x_new = np.insert(x_new, 0, x_ins)
    dtot = d(x)
    d_ins = np.zeros(11)
    dtot = np.insert(dtot, 0, d_ins)


    plot(x_new, dtot, marker = 'x')


    # =============================================================================
    # Transmission calculations
    # =============================================================================
    # P, p = transm_from_x_tcu(x_new, dtot, und = 'ivu')
    # with open('transmission_data.npy', 'wb') as f:
    #     np.save(f, x_new)
    #     np.save(f, p)
    #     np.save(f, P)
    with open('transmission_data_cpmu.npy', 'rb') as f:
        x_new_c = np.load(f)
        p_cpmu = np.load(f)
        P_cpmu = np.load(f)
    with open('transmission_data_ivu.npy', 'rb') as f:
        x_new_i = np.load(f)
        p_ivu = np.load(f)
        P_ivu = np.load(f)
        
    # Plot absorbed power density at and image of guard
    
    
    fig, ax = plt.subplots(2,1, sharex = True)
    ax[0].plot(x_new_i, p_ivu, label = 'IVU21')
    ax[0].plot(x_new_c, p_cpmu, label = 'CPMU16.4')
    
    
    
    xlim = ax[0].get_xlim()
    dx = xlim[1] - xlim[0]
    
    
    ax[0].set_title('Power transmitted through Copper guard - 11 keV')
    ax[0].legend()
    
    image = plt.imread('29311075.PNG')
    ratio = image.shape[0]/image.shape[1]
    ax[1].imshow(image, alpha = 0.5, extent = (x0,x0 + 8,-8*ratio,0), aspect='auto')
    ax[1].set_xlim(xlim)
    ax[1].axis('off')
    ax[0].xaxis.set_tick_params(labelbottom=True)
    
    
    ax[0].set_xlabel('Horizontal x ($mm$)')
    ax[0].set_ylabel('Transm. power density p ($W/mm^2$)')
    ax[1].annotate("", xytext=(x0 + 0.1, -8*ratio), xy=(x0 + 0.1, -8*ratio+2),
            arrowprops=dict(arrowstyle="->", color = 'tab:red', lw = 3))
    ax[1].text(x0 + 0.15, -8*ratio+1, 'x-rays', color = 'tab:red')
    ax[1].text(x0 + 0.3, -8*ratio+3, 'copper guard', color = 'black')
    
    data_cpmu = np.zeros((x_new.size, 2))
    data_cpmu[:,0], data_cpmu[:,1] = x_new_c, p_cpmu
    data_ivu = np.zeros((x_new.size, 2))
    data_ivu[:,0], data_ivu[:,1] = x_new_i, p_ivu
    # np.savetxt('p_frontGuard_cpmu.csv', data_cpmu, delimiter=',')
    # np.savetxt('p_frontGuard_ivu.csv', data_ivu, delimiter=',')
    
    # =============================================================================
    # Spectral calculations
    # =============================================================================
    # # Check good energy discretization
    # n = [1000, 1500, 2000, 3000]
    # Pn = []
    # for ni in n:
    #     print('running ni={:.0f}'.format(ni))
    #     energy, sp_pow = spectral_calc(1.2, n = ni, emax = 1.5e5)
    #     Pi = np.trapz(sp_pow, energy)
    #     Pn.append(Pi)
        
    # # Check good max. energy discretization
    # emax = [1e5, 1.2e5, 1.4e5, 1.6e5, 1.8e5, 2e5, 2.5e5]
    # Pemax = []
    # for emaxi in emax:
    #     print('running emaxi={:.0f}'.format(emaxi))
    #     energy, sp_pow = spectral_calc(1.2, emax = emaxi, n = 2000)
    #     Pi = np.trapz(sp_pow, energy)
    #     Pemax.append(Pi)
    
    