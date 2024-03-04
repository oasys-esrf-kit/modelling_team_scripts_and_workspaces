
import matplotlib.pyplot as plt
import pandas as pd
import numpy
import xraylib
import scipy.constants as codata
from syned.util.json_tools import load_from_json_file
from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_spectrum
from xoppylib.power.xoppy_calc_power import xoppy_calc_power
from srxraylib.plot.gol import plot

def get_spectral_power(res_energy, position='primary_slits', und='U17.6_4.8m', plot_spec=False):
    
    #loads the undulator characteristics from the syned JSON
    try:
        syned_json_file = 'ESRF_ID28_EBS_'+ und + '.json'
        syned_obj = load_from_json_file(syned_json_file)
        e = syned_obj.get_electron_beam()        
        u = syned_obj.get_magnetic_structure()
    except:
        print('Unable to load syned file', syned_json_file)    
    
    if und == 'U17.6_4.8m':
        harmonic = 1
        e_min = 3000
        e_max = 120000        
    elif und == 'U32_4.8m':
        e_min = 1000
        e_max = 100000
        if res_energy < 16000:
            harmonic = 3
        else:
            harmonic = 5
    elif und == 'IVU22':
        e_min = 1000
        e_max = 100000
        if res_energy <= 15000:
            harmonic = 1
        else:
            harmonic = 3
    elif und == 'IVU13.4':
        e_min = 1000
        e_max = 100000
        harmonic = 1    
    else:
        RuntimeError(f'ERROR: script not implemented for undulator {und}') 
    
    if position == 'primary_slits':
        distance = 27.0
        h_gap = 0.00233
        v_gap = 0.0046
        e_points = 1000
    elif position == 'Be_lenses':
        distance = 28.5
        h_gap = 0.00106
        v_gap = 0.00106
        e_points = 2000
    else:
        RuntimeError(f'ERROR: script not implemented for this position {position}')

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
            DISTANCE=distance,
            GAPH=h_gap,
            GAPV=v_gap,
            GAPH_CENTER=0.0,
            GAPV_CENTER=0.0,
            PHOTONENERGYMIN=e_min,
            PHOTONENERGYMAX=e_max,
            PHOTONENERGYPOINTS=e_points,
            METHOD=2,
            USEEMITTANCES=1)
    
    if plot_spec:
        plot(photon_energy,spectral_power, xtitle="Photon energy [eV]",
             ytitle="Power [W/eV]",title="Undulator Spectral Power",
             xlog=False,ylog=False,show=False)

    return photon_energy, spectral_power

def pow_transport(photon_energy, spectral_power, position='primary_slits'):     

    #here we check if the the attenuator is active and over its on energy, example atten_active = 12000 eV   
    if position == 'primary_slits':
        #total power is the power arriving to the primary slits, abs_power is the absorbed power by the diamond window
        out_dictionary = xoppy_calc_power(
            photon_energy,
            spectral_power,
            substance = ['C',],
            thick     = [0.3,], # in mm (for filters)
            angle     = [3,], # in mrad (for mirrors)
            dens      = ['3.52',],
            roughness = [0,], # in A (for mirrors)
            flags     = [0,], # 0=Filter, 1=Mirror
            nelements = 1,
            FILE_DUMP = 0,
            material_constants_library = xraylib,
            )
        total_power = round(numpy.trapz(out_dictionary["data"][-1,:], x=photon_energy, axis=-1)) 
        abs_power = round(numpy.trapz(out_dictionary["data"][-2,:], x=photon_energy, axis=-1))

    elif position == 'Be_lenses':
        #total power is the power after the primary slits, abs_power is the absorbed power by 1 mm3 of Be
        out_dictionary = xoppy_calc_power(
            photon_energy,
            spectral_power,
            substance = ['C','Be',],
            thick     = [0.3,1,], # in mm (for filters)
            angle     = [3,3,], # in mrad (for mirrors)
            dens      = ['3.52','?',],
            roughness = [0,0,], # in A (for mirrors)
            flags     = [0,0,], # 0=Filter, 1=Mirror
            nelements = 2,
            FILE_DUMP = 0,
            material_constants_library = xraylib,
            )
        total_power = round(numpy.trapz(out_dictionary["data"][7,:], x=photon_energy, axis=-1)) 
        abs_power = round(numpy.trapz(out_dictionary["data"][-2,:], x=photon_energy, axis=-1))

    else:
        RuntimeError(f'ERROR: script not implemented for this position {position}')    
    
    return total_power, abs_power

def run_calculations(res_energies, und='17.6_4.8m', save_file=True):
    
    pow_at_pp = []
    wind_abs_power = []
    pow_after_pp = []
    pow_abs_be = []
    for res_energy in res_energies:
        energy1, spec_pow1 = get_spectral_power(res_energy, und=und, position='primary_slits', plot_spec=False)
        pow1, abs1 = pow_transport(energy1, spec_pow1, position='primary_slits')
        pow_at_pp.append(pow1)
        wind_abs_power.append(abs1)
        energy2, spec_pow2 = get_spectral_power(res_energy, und=und, position='Be_lenses', plot_spec=False)
        pow2, abs2 = pow_transport(energy2, spec_pow2, position='Be_lenses')
        pow_after_pp.append(pow2)
        pow_abs_be.append(abs2)

    
    if save_file:
        df_pow = pd.DataFrame({'res_energy[eV]':res_energies, 'Pow_on_PS[W]':pow_at_pp, \
            'Abs_window[W]':wind_abs_power, 'Power_after_PP[W]': pow_after_pp, 'Power_abs_BE_cmm[W]': pow_abs_be})
        df_pow.to_csv(f'id28_pow_{und}.csv', index=False)
        print(f'File id28_pow_{und}.csv saved on disk')

    return df_pow            

if __name__=="__main__":    
    #energ1 = [12400, 15000, 21700, 23700]
    #run_calculations(energ1, und='U32_4.8m')
    #energ2 = [17800]
    #run_calculations(energ2, und='U17.6_4.8m')
    #energ3 = [12400, 15000, 17800, 21700, 23700]
    #run_calculations(energ3, und='IVU22')
    energ4 = [21700, 23700]
    run_calculations(energ4, und='IVU13.4')

    #pass