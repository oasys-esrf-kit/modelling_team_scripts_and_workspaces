import matplotlib.pyplot as plt
import pandas as pd
import numpy
import xraylib
import scipy.constants as codata
from syned.util.json_tools import load_from_json_file
from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_spectrum
from xoppylib.power.xoppy_calc_power import xoppy_calc_power
from srxraylib.plot.gol import plot

def get_spectral_power(res_energy, und='CPMU19_2.5m', plot_spec=False):
    
    try:
        syned_json_file = 'ESRF_ID20_EBS_'+ und + '.json'
        syned_obj = load_from_json_file(syned_json_file)
        e = syned_obj.get_electron_beam()        
        u = syned_obj.get_magnetic_structure()
    except:
        print('Unable to load syned file', syned_json_file)    
    
    if und == 'CPMU19_2.5m':
        first_harm_max_energ = 15300
    elif und == 'CPMU20.5_2.5m':
        first_harm_max_energ = 14100
    elif und == 'CPMU20.5_4m':
        first_harm_max_energ = 14600
    else:
        RuntimeError(f'ERROR: script not implemente for undulator {und}')
    
    harmonic=1 if res_energy < first_harm_max_energ else 3

    (sigma_x, sig_div_x, sigma_y, sig_div_y) = e.get_sigmas_all()

    photon_energy,_,spectral_power,_ = xoppy_calc_undulator_spectrum(
            ELECTRONENERGY=e.energy(),
            ELECTRONENERGYSPREAD=e._energy_spread,
            ELECTRONCURRENT=e.current(),
            ELECTRONBEAMSIZEH=sigma_x,
            ELECTRONBEAMSIZEV=sigma_y,
            ELECTRONBEAMDIVERGENCEH=sig_div_x,
            ELECTRONBEAMDIVERGENCEV=sig_div_y,
            PERIODID=u.period_length(),
            NPERIODS=u.number_of_periods(),
            KV=u.get_K_from_photon_energy(res_energy, e.gamma(), harmonic=harmonic),
            KH=0.0,
            KPHASE=0.0,
            DISTANCE=31.0,
            GAPH=0.00084,
            GAPV=0.00084,
            GAPH_CENTER=0.0,
            GAPV_CENTER=0.0,
            PHOTONENERGYMIN=1000,
            PHOTONENERGYMAX=150000,
            PHOTONENERGYPOINTS=4000,
            METHOD=2,
            USEEMITTANCES=1)
    
    if plot_spec:
        plot(photon_energy,spectral_power, xtitle="Photon energy [eV]",
             ytitle="Power [W/eV]",title="Undulator Spectral Power",
             xlog=False,ylog=False,show=False)

    return photon_energy, spectral_power

def pow_transport_mirror(res_energy, photon_energy, spectral_power, atten_active=None):

    if res_energy < 10000:
        mirror_coating = 'Si'
        mc_dens = '?'
    elif 10000 <= res_energy < 12000:
        mirror_coating = 'SiC'
        mc_dens = '3.21'
    else:
        mirror_coating = 'Rh'
        mc_dens = '?' 

    #here we check if the the attenuator is active and over its on energy, example atten_active = 12000 eV   
    if atten_active and res_energy > atten_active:
        
        out_dictionary = xoppy_calc_power(
        photon_energy,
        spectral_power,
        substance = ['C','Cr',mirror_coating,],
        thick     = [0.05,0.0001,0.5,], # in mm (for filters)
        angle     = [3.1,3,3.1,], # in mrad (for mirrors)
        dens      = ['3.52','?',mc_dens,],
        roughness = [0,0,0,], # in A (for mirrors)
        flags     = [0,0,1,], # 0=Filter, 1=Mirror
        nelements = 3,
        FILE_DUMP = 0,
        material_constants_library = xraylib,
        )

    else:

        out_dictionary = xoppy_calc_power(
            photon_energy,
            spectral_power,
            substance = [mirror_coating,],        
            thick     = [0.5,], # in mm (for filters)
            angle     = [3.1,], # in mrad (for mirrors)
            dens      = [mc_dens,],        
            roughness = [0,], # in A (for mirrors)
            flags     = [1,], # 0=Filter, 1=Mirror
            nelements = 1,
            FILE_DUMP = 0,
            material_constants_library = xraylib,
            )
    pow_on_dcm = round(numpy.trapz(out_dictionary["data"][-1,:], x=photon_energy, axis=-1), 2)

    return pow_on_dcm

def run_calculations(res_energies, und='CPMU19_2.5m', save_file=True, plot_pow=True, atten_active=None):
    
    pow_on_dcm = []
    for res_energy in res_energies:
        photon_energy, spec_pow = get_spectral_power(res_energy, und=und, plot_spec=False)
        pow = pow_transport_mirror(res_energy, photon_energy, spec_pow, atten_active=atten_active)
        pow_on_dcm.append(pow)
    
    if save_file:
        df_pow = pd.DataFrame({'res_energy[eV]':res_energies, 'power_on_dcm[W]':pow_on_dcm})
        if atten_active:
            df_pow.to_csv(f'id20_pow_on_dcm_{und}_att.csv', index=False)
            print(f'File id20_pow_on_dcm_{und}_att.csv saved on disk')
        else:       
            df_pow.to_csv(f'id20_pow_on_dcm_{und}.csv', index=False)
            print(f'File id20_pow_on_dcm_{und}.csv saved on disk')

    if plot_pow:
        plot(res_energies,pow_on_dcm, xtitle="Photon energy [eV]",
             ytitle="Total power [W]", title=f"Power on DCM due {und}",
             xlog=False,ylog=False,show=False)

if __name__=="__main__":
    res_energies = numpy.linspace(4700, 20000, 200)    
    run_calculations(res_energies, und='CPMU19_2.5m', save_file=True, plot_pow=False, atten_active=12000)
    print('done for CPMU19_2.5m')

    res_energies = numpy.linspace(4100, 20000, 200)    
    run_calculations(res_energies, und='CPMU20.5_2.5m', save_file=True, plot_pow=False, atten_active=12000)
    print('done for CPMU20.5_2.5m')

    res_energies = numpy.linspace(4800, 20000, 200)    
    run_calculations(res_energies, und='CPMU20.5_4m', save_file=True, plot_pow=False, atten_active=12000)
    print('done for CPMU20.5_4m')