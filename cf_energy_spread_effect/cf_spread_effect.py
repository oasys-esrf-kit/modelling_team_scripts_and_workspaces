#
# Import section
#
import numpy as np
import scipy.constants as codata
from syned.util.json_tools import load_from_json_file
from wofryimpl.propagator.util.undulator_coherent_mode_decomposition_1d import UndulatorCoherentModeDecomposition1D   
import matplotlib.pyplot as plt
import pandas as pd
#from srxraylib.plot.gol import plot, plot_image


def get_cf_cmd_es(photon_energy, json_file, energy_spread=0, harmonic=1, scan_direction='H'):

    """ Function to get the Coherence fraction for a given photon energy, an
    ID-JSON file and energy_spread using Coherent Mode Descompositon """
    # gets the ID and Electron beam properties from the JSON file
    syned_obj = load_from_json_file(json_file)
    u = syned_obj.get_magnetic_structure()
    e = syned_obj.get_electron_beam()    

    #print(f'Undulator period lenght {u.period_length()} m')    

    lambdan = (codata.h / codata.electron_volt * codata.c / photon_energy)
    k = round(np.sqrt(2*(((lambdan*2*e.gamma()**2)/u.period_length())-1)), 6)

    if scan_direction == 'H':
        sigmaxx    = e.get_sigmas_horizontal()[0]
        sigmaxpxp  = e.get_sigmas_horizontal()[1]
    elif scan_direction == 'V':
        sigmaxx    = e.get_sigmas_vertical()[0]
        sigmaxpxp  = e.get_sigmas_vertical()[1]
    else:
        raise RuntimeError("ERROR: Unidentified scan direction")

    if energy_spread == 0:
        flag_e_sp = 0
    else:
        flag_e_sp = 1        

    print(f'Energy {photon_energy} corresponds a K value of {k}')

    coherent_mode_decomposition = UndulatorCoherentModeDecomposition1D(
        electron_energy=e.energy(),
        electron_current=e.current(),
        undulator_period=u.period_length(),
        undulator_nperiods=u.number_of_periods(),
        K=k,
        photon_energy=photon_energy,
        abscissas_interval=0.00025,
        number_of_points=1000,
        distance_to_screen=100,
        magnification_x_forward=100,
        magnification_x_backward=0.01,
        scan_direction=scan_direction,
        sigmaxx=sigmaxx,
        sigmaxpxp=sigmaxpxp,
        useGSMapproximation=False,
        e_energy_dispersion_flag=flag_e_sp,
        e_energy_dispersion_sigma_relative=energy_spread,
        e_energy_dispersion_interval_in_sigma_units=6,
        e_energy_dispersion_points=15)
    
    coherent_mode_decomposition.calculate()    
    cf_cmd = coherent_mode_decomposition.eigenvalues[0] / coherent_mode_decomposition.eigenvalues.sum()    
    print(f"Coherent Fraction from CMD is {round(cf_cmd, 6)}")
    return round(cf_cmd, 6)

def run_calculations(photon_energy, energies_spread, id='cpmu18', harmonic=1, scan_direction = 'H', save_file=False):
    cfs_cmd = []

    if id == 'cpmu18':
        json_file = 'EBS_ID18_CPMU18_on-energy.json'
    else:
        #to be impremented for other IDs and to consider other than the fisrt harmonic
        pass

    for energy_spread in energies_spread:
        cfs_cmd.append(get_cf_cmd_es(photon_energy, json_file, energy_spread=energy_spread, harmonic=harmonic, scan_direction=scan_direction))  
        
    df_cfs = pd.DataFrame({'Energy_spread':energies_spread,'CFs_CMD':cfs_cmd})
    
    if save_file:
        df_cfs.to_csv(f'{id}_{photon_energy}_cfs_results_{scan_direction}.csv', index=False)
        print(f'file {id}_{photon_energy}_cfs_results_{scan_direction}.csv has been saved to disk')

    return df_cfs

def plot_results(csv_files, plot_type='energy_spread'):

    """ Plots the CF from the created CSV files, please give the filenames as a list"""
    
    f_size = 14

    plt.ion()
    plt.figure()

    if plot_type == 'energy_spread':

        for csv_file in csv_files:
    
            if "H" in csv_file:
                label = 'Horizontal'
            elif "V" in csv_file:
                label = 'Vertical'
            else:
                raise RuntimeError("ERROR: Unidentified scan direction")
    
            df = pd.read_csv(csv_file, sep=',|\s+', comment = '#', engine='python')
            plt.plot(df['Energy_spread']*1e3, df['CFs_CMD'], 'o-', label=label)        
    
        plt.xlabel("Energy Spread (x10$^{-3}$)", fontsize=f_size)
        plt.ylabel("Coherence fraction", fontsize=f_size)    
        
        plt.ylim(0.0, 0.6)
        plt.title("EBS CPMU18 at 7 keV ", fontsize=f_size)

    elif plot_type == 'detuning':

        for csv_file in csv_files:
    
            if "H" in csv_file:
                label = 'Horizontal'
            elif "V" in csv_file:
                label = 'Vertical'
            else:
                raise RuntimeError("ERROR: Unidentified scan direction")
    
            df = pd.read_csv(csv_file, sep=',|\s+', comment = '#', engine='python')
            plt.plot(df['Photon_energy'], df['CFs_CMD'], 'o-', label=label) 
        
        plt.xlabel("Photon Energy (eV)", fontsize=f_size)
        plt.ylabel("Coherence fraction (a. u.)", fontsize=f_size)
        plt.title("EBS CPMU18 resonance at 7 keV ", fontsize=f_size)
    
    else:
        raise RuntimeError("ERROR: Unidentified plotting type")    

    #plt.xscale('log')
    plt.legend(fontsize=f_size)
    
    
    plt.grid(which='both', axis='both')
    plt.show()


def get_cf_cmd_det(photon_energy, json_file, detunig=False, energy_spread=0, resonance_energy=7000, harmonic=1, scan_direction='H'):

    """ Function to get the Coherence fraction for a given photon energy, an
    ID-JSON file and energy_spread using Coherent Mode Descompositon """
    # gets the ID and Electron beam properties from the JSON file
    syned_obj = load_from_json_file(json_file)
    u = syned_obj.get_magnetic_structure()
    e = syned_obj.get_electron_beam()    

    #print(f'Undulator period lenght {u.period_length()} m')    
    if detunig:
        lambdan = (codata.h / codata.electron_volt * codata.c / resonance_energy)
        k = round(np.sqrt(2*(((lambdan*2*e.gamma()**2)/u.period_length())-1)), 6)
    else:
        lambdan = (codata.h / codata.electron_volt * codata.c / resonance_energy)
        k = round(np.sqrt(2*(((lambdan*2*e.gamma()**2)/u.period_length())-1)), 6)
        
    if scan_direction == 'H':
        sigmaxx    = e.get_sigmas_horizontal()[0]
        sigmaxpxp  = e.get_sigmas_horizontal()[1]
    elif scan_direction == 'V':
        sigmaxx    = e.get_sigmas_vertical()[0]
        sigmaxpxp  = e.get_sigmas_vertical()[1]
    else:
        raise RuntimeError("ERROR: Unidentified scan direction")

    if energy_spread == 0:
        flag_e_sp = 0
    else:
        flag_e_sp = 1        

    print(f'Energy {photon_energy} corresponds a K value of {k}')

    coherent_mode_decomposition = UndulatorCoherentModeDecomposition1D(
        electron_energy=e.energy(),
        electron_current=e.current(),
        undulator_period=u.period_length(),
        undulator_nperiods=u.number_of_periods(),
        K=k,
        photon_energy=photon_energy,
        abscissas_interval=0.00025,
        number_of_points=1000,
        distance_to_screen=100,
        magnification_x_forward=100,
        magnification_x_backward=0.01,
        scan_direction=scan_direction,
        sigmaxx=sigmaxx,
        sigmaxpxp=sigmaxpxp,
        useGSMapproximation=False,
        e_energy_dispersion_flag=flag_e_sp,
        e_energy_dispersion_sigma_relative=energy_spread,
        e_energy_dispersion_interval_in_sigma_units=6,
        e_energy_dispersion_points=15)
    
    coherent_mode_decomposition.calculate()    
    cf_cmd = coherent_mode_decomposition.eigenvalues[0] / coherent_mode_decomposition.eigenvalues.sum()    
    print(f"Coherent Fraction from CMD is {round(cf_cmd, 6)}")
    return round(cf_cmd, 6)


def run_calc_detuning(photon_energies, energy_spread=0.001, id='cpmu18', detunig=False, harmonic=1, resonance_energy=7000, scan_direction = 'H', save_file=False):
    cfs_cmd = []

    if id == 'cpmu18':
        json_file = 'EBS_ID18_CPMU18_on-energy.json'
    else:
        #to be impremented for other IDs and to consider other than the fisrt harmonic
        pass

    for photon_energy in photon_energies:
        cfs_cmd.append(get_cf_cmd_det(photon_energy, json_file, detunig=detunig, energy_spread=energy_spread, resonance_energy=resonance_energy, harmonic=harmonic, scan_direction=scan_direction))  
        
    df_cfs = pd.DataFrame({'Photon_energy':photon_energies,'CFs_CMD':cfs_cmd})
    
    if save_file:
        df_cfs.to_csv(f'{id}_{resonance_energy*1e-3}keV_cfs_results_e-spread_{energy_spread}_{scan_direction}.csv', index=False)
        print(f'file {id}_{resonance_energy*1e-3}keV_cfs_results_e-spread_{energy_spread}_{scan_direction}.csv has been saved to disk')

    return df_cfs

def compare_prof(profile_files):

    f_size = 14

    plt.ion()
    plt.figure()

    for profile_file in profile_files:
    
            if "h" in profile_file:
                legend = 'Horizontal beam profiles at 30 m'
                h_axis = 'Horizontal axis (mm)'

                if '7000' in profile_file:
                    label = '7000 eV'
                elif '7020' in profile_file:
                    label = '7020 eV'
                else:
                    raise RuntimeError("ERROR: Not energy found")

            elif "v" in profile_file:
                legend = 'Vertical beam profiles at 30 m'
                h_axis = 'Vertical axis (mm)'

                if '7000' in profile_file:
                    label = '7000 eV'
                elif '7020' in profile_file:
                    label = '7020 eV'
                else:
                    raise RuntimeError("ERROR: Not energy found")
            else:
                raise RuntimeError("ERROR: Unidentified scan direction")
    
            df = pd.read_csv(profile_file, sep=',|\s+', comment = '#', engine='python')
            plt.plot(df.iloc[:, 0], df.iloc[:, 1]/df.iloc[:, 1].max(), 'o-', label=label)
    
    plt.xlabel(h_axis, fontsize=f_size)
    plt.ylabel("Normalized photon flux", fontsize=f_size)
    plt.title(legend, fontsize=f_size)

    plt.legend(fontsize=f_size)
    
    
    plt.grid(which='both', axis='both')
    plt.show()

    

if __name__ == "__main__":
    #energies_spread = np.arange(0, 0.0065, 0.0005)
    #run_calculations(7000, energies_spread, id='cpmu18', harmonic=1, scan_direction = 'V', save_file = True)
    #plot_results(['cpmu18_7000_cfs_results_H.csv', 'cpmu18_7000_cfs_results_V.csv'])
    pass
    #photon_energies = np.arange(6800, 7150, 10)
    #run_calc_detuning(photon_energies, energy_spread=0.001, id='cpmu18', detunig=True, resonance_energy=7000, scan_direction='H', save_file=True)
    #run_calc_detuning(photon_energies, energy_spread=0.001, id='cpmu18', detunig=True, resonance_energy=7000, scan_direction='V', save_file=True)