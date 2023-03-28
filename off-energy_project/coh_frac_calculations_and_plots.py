import numpy as np
import scipy.constants as codata
from syned.util.json_tools import load_from_json_file
from wofryimpl.propagator.util.undulator_coherent_mode_decomposition_1d import UndulatorCoherentModeDecomposition1D   
import matplotlib.pyplot as plt
import pandas as pd

def gammas_rad(json_file):
    syned_obj = load_from_json_file(json_file)
    u = syned_obj.get_magnetic_structure()
    e = syned_obj.get_electron_beam()
    return u.get_sigmas_radiation(e.gamma())

def photon_source(photon_energy, json_file):

    """ Function to get the radiation beam size for a given energy
    and an ID - JSON file """

    syned_obj = load_from_json_file(json_file)
    u_length = syned_obj.get_magnetic_structure().length()
    # calculate sizes of the photon undulator beam
    # see formulas 25 & 30 in Elleaume (Onuki & Elleaume)
    lambdan = (codata.h / codata.electron_volt * codata.c / photon_energy)
    photon_size= (2.740 / 4 / np.pi) * np.sqrt(lambdan * u_length)
    photon_div = 0.69 * np.sqrt(lambdan / u_length)    

    return photon_size, photon_div

def conv_source(photon_energy, json_file):
    """ Function to convolve the electron beam size for a given energy and
    an ID - JSON file """
    syned_obj = load_from_json_file(json_file)
    e_s = syned_obj.get_electron_beam().get_sigmas_all()
    #pure photon
    p_s = photon_source(photon_energy, json_file)	
    #conv source
    conv_sx = np.sqrt((e_s[0])**2 + (p_s[0])**2)
    conv_divx = np.sqrt((e_s[1])**2 + (p_s[1])**2)

    conv_sy = np.sqrt((e_s[2])**2 + (p_s[0])**2)
    conv_divy = np.sqrt((e_s[3])**2 + (p_s[1])**2)   

    return conv_sx, conv_divx, conv_sy, conv_divy

def get_cf_from_alg(photon_energy, json_file):

    """ Function to get the Coherence fraction for a given photon energy
    and an ID-JSON file using the classical algebraic expression """

    p_s = photon_source(photon_energy, json_file)
    conv_s = conv_source(photon_energy, json_file)

    cf_alg = round((p_s[0] * p_s[1]) / (conv_s[0] * conv_s[1]), 6)

    print(f'Algebraic CF_for on-energy is {cf_alg}')

    return cf_alg

def get_cf_from_cmd(photon_energy, json_file):

    """ Function to get the Coherence fraction for a given photon energy
    and an ID-JSON file using Coherent Mode Descompositon """

    syned_obj = load_from_json_file(json_file)
    u = syned_obj.get_magnetic_structure()
    e = syned_obj.get_electron_beam()    

    #print(f'Undulator period lenght {u.period_length()} m')

    lambdan = (codata.h / codata.electron_volt * codata.c / photon_energy)
    k = round(np.sqrt(2*(((lambdan*2*e.gamma()**2)/u.period_length())-1)), 6)

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
        scan_direction='H',
        sigmaxx=e.get_sigmas_horizontal()[0],
        sigmaxpxp=e.get_sigmas_horizontal()[1],
        useGSMapproximation=False,)
        # make calculation
    coherent_mode_decomposition.calculate()    
    cf_cmd = coherent_mode_decomposition.eigenvalues[0] / coherent_mode_decomposition.eigenvalues.sum()    
    print(f"Coherent Fraction from CMD is {round(cf_cmd, 6)}")
    return round(cf_cmd, 6)

def plot_results(csv_file):
    f_size = 14

    df = pd.read_csv(csv_file, sep=',|\s+', comment = '#', engine='python')

    plt.ion()
    plt.figure()

    plt.plot(df['Photon_energy'], df['CF_Alg_On'], label='CF$^{Alg}$')
    plt.plot(df['Photon_energy'], df['CF_CMD_On'], label='CF$^{CMD}$')

    plt.xlabel("Photon energy [eV]",fontsize=f_size)
    plt.ylabel("Coherence fraction (a.u.)",fontsize=f_size)

    plt.legend(fontsize=f_size)
    plt.title("On-energy", fontsize=f_size)

    #plt.xscale('log')
    #plt.xlim(0.1, 10)
    plt.grid(which='both', axis='both')
    plt.show()

    plt.ion()
    plt.figure()

    plt.plot(df['Photon_energy'], df['CF_Alg_Off'], label='CF$^{Alg}$')
    plt.plot(df['Photon_energy'], df['CF_CMD_Off'], label='CF$^{CMD}$')

    plt.xlabel("Photon energy [eV]",fontsize=f_size)
    plt.ylabel("Coherence fraction (a.u.)",fontsize=f_size)

    plt.legend(fontsize=f_size)
    plt.title("Off-energy", fontsize=f_size)

    #plt.xscale('log')
    #plt.xlim(0.1, 10)
    plt.grid(which='both', axis='both')
    plt.show()

    plt.ion()
    plt.figure()


    plt.plot(df['Photon_energy'], df['CF_Alg_Off']/df['CF_Alg_On']*100-100, label='Algebraic')
    plt.plot(df['Photon_energy'], df['CF_CMD_Off']/df['CF_CMD_On']*100-100, label='CMD')

    plt.xlabel("Photon energy [eV]",fontsize=f_size)
    plt.ylabel("CF increment [%]",fontsize=f_size)

    plt.legend(fontsize=f_size)

    #plt.xscale('log')
    #plt.xlim(0.1, 10)
    plt.grid(which='both', axis='both')
    plt.show()

def compare_plots_occupation(oc_file1, oc_file2):
    f_size = 14

    df1 = pd.read_csv(oc_file1, sep=',|\s+', comment = '#', header=None, engine='python')
    df2 = pd.read_csv(oc_file2, sep=',|\s+', comment = '#', header=None, engine='python')

    x1 = np.array(df1.iloc[:, 0])
    #print(x1)
    y1 = np.array(df1.iloc[:, 1])

    x2 = np.array(df2.iloc[:, 0])
    y2 = np.array(df2.iloc[:, 1])

    plt.ion()
    plt.figure()


    plt.plot(x1, y1, 'ro-', label='On-energy')
    plt.plot(x2, y2, 'bo-', label='Off-energy')    

    plt.xlabel("Mode Index", fontsize=f_size)
    plt.ylabel("Occupation", fontsize=f_size)

    plt.legend(fontsize=f_size)
    plt.title("7 keV photon energy", fontsize=f_size)

    #plt.xscale('log')
    #plt.xlim(0.1, 10)
    plt.grid(which='both', axis='both')
    plt.show()

def compare_plots_size(oc_file1, oc_file2):
    f_size = 14

    df1 = pd.read_csv(oc_file1, sep=',|\s+', comment = '#', header=None, engine='python')
    df2 = pd.read_csv(oc_file2, sep=',|\s+', comment = '#', header=None, engine='python')

    x1 = np.array(df1.iloc[:, 0])
    #print(x1)
    y1 = np.array(df1.iloc[:, 1])

    x2 = np.array(df2.iloc[:, 0])
    y2 = np.array(df2.iloc[:, 1])

    plt.ion()
    plt.figure()


    plt.plot(x1, y1, 'r', label='On-energy')
    plt.plot(x2, y2, 'b', label='Off-energy')    

    plt.xlabel("Horizontal axis [$\mu$m]", fontsize=f_size)
    plt.ylabel("Intensity", fontsize=f_size)

    plt.legend(fontsize=f_size)
    plt.title("7 keV photon energy", fontsize=f_size)

    #plt.xscale('log')
    #plt.xlim(0.1, 10)
    plt.grid(which='both', axis='both')
    plt.show()


def run_calculations():

    photon_energies = np.arange(7000, 19000, 500)

    cf_alg_on = []
    cf_alg_off = []
    cf_cmd_on = []
    cf_cmd_off = []

    for e in photon_energies:
        cf_alg_on.append(get_cf_from_alg(e, 'EBS_ID18_CPMU18_on-energy.json'))
        cf_alg_off.append(get_cf_from_alg(e, 'EBS_ID18_CPMU18_off-energy.json'))
        cf_cmd_on.append(get_cf_from_cmd(e, 'EBS_ID18_CPMU18_on-energy.json'))
        cf_cmd_off.append(get_cf_from_cmd(e, 'EBS_ID18_CPMU18_off-energy.json'))

    df_cfs = pd.DataFrame({'Photon_energy':photon_energies,'CF_Alg_On':cf_alg_on, 'CF_CMD_On':cf_cmd_on, 'CF_Alg_Off':cf_alg_off, 'CF_CMD_Off':cf_cmd_off})
    df_cfs.to_csv('cpmu18_cf_results.csv', index=False)
    print('file cpmu18_cf_results.csv has been saved to disk')


if __name__ == "__main__":
    pass
    
    
    

        


    
