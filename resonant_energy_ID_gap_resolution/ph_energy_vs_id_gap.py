import os, sys
import numpy
import scipy.constants as codata
import pandas as pd
import matplotlib.pyplot as plt

def get_id(section, id_name):

    df = pd.read_csv('jsrund.csv', skiprows=1, header=1)

    id_df = df[(df["Straight Section"] == section) & (df["Name"] == id_name)]

    if id_df.empty == True:
        raise RuntimeError("WARNING: There is no ID in the given section")
    
    return id_df

def gap_to_res_energy(id_gap_mm, section=19, id_name='U32a', harmonic=1):
    # CPMU16_4, section = 3, id_name='CPMU16-4b'
    # CPMU18, section = 6, id_name='CPMU18a' 
    # U27, section = 24, id_name='U27a'
    # U32, section = 19, id_name='U32c'
    
    electron_energy_in_GeV = 6.0
    id_df = get_id(section, id_name)

    period_length_mm = float(id_df['Period [mm]'])
    a0 = float(id_df['A(0)'])
    a1 = float(id_df['A(1)'])

    if "CPMU" in id_name:
        a2 = float(id_df['A(2)'])
        a3 = float(id_df['A(3)'])
        a4 = float(id_df['A(4)'])
        a5 = float(id_df['A(5)'])
        B  = a0 * numpy.exp(-numpy.pi * a3 * 1 * id_gap_mm / period_length_mm)
        B += a1 * numpy.exp(-numpy.pi * a4 * 2 * id_gap_mm / period_length_mm)
        B += a2 * numpy.exp(-numpy.pi * a5 * 3 * id_gap_mm / period_length_mm)
        
    elif "HU" in id_name:
        raise RuntimeError("ERROR: Helical/Apple undulators not implemented in this app (wrong results")
         
    else:
        B = a0 * numpy.exp(-numpy.pi * a1 * id_gap_mm / period_length_mm)

    K = B * (period_length_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)

    gamma = 1e9 * electron_energy_in_GeV / (codata.m_e *  codata.c**2 / codata.e)

    resonance_wavelength =  ((period_length_mm * 1e-3 / (2.0 * gamma **2)) * (1 + K**2 / 2.0))/harmonic
    
    resonance_energy =  codata.c * codata.h / codata.e / resonance_wavelength
    
    print(f'For a gap of {id_gap_mm} mm, for the {id_name} the Resonance energy is {resonance_energy} eV')

    return resonance_energy

def run_calcualtions(section=3, id_name='CPMU16-4b', save_file=True):
    # CPMU16_4, section = 3, id_name='CPMU16-4b'
    # CPMU18, section = 6, id_name='CPMU18a' 
    # U27, section = 24, id_name='U27a'
    # U32, section = 19, id_name='U32c'
    
    id_df = get_id(section, id_name)

    min_gap = float(id_df['Gap min. [mm]'])

    id_gaps = numpy.arange(min_gap, min_gap + 0.101, 0.001) #Scan is done higher

    res_energy1 = []
    res_energy3 = []
    res_energy5 = []

    for id_gap_mm in id_gaps:

        res_energy1.append(gap_to_res_energy(id_gap_mm, section=section, id_name=id_name, harmonic=1))
        res_energy3.append(gap_to_res_energy(id_gap_mm, section=section, id_name=id_name, harmonic=3))
        res_energy5.append(gap_to_res_energy(id_gap_mm, section=section, id_name=id_name, harmonic=5))
    
    df_res_e = pd.DataFrame({'ID_gap':id_gaps, 'harm1':res_energy1, 'harm3':res_energy3, 'harm5':res_energy5})
    
    if save_file:

        df_res_e.to_csv(f'Resonant_energy_results_{id_name}.csv')    
        print(f'File Resonant_energy_results_{id_name}.csv ahs been saved to disk')

    return df_res_e

def plot_results(csv_file, plot_harmonic=1):
    f_size = 14
    df = pd.read_csv(csv_file, sep=',|\s+', comment = '#', engine='python')

    diff1 = []
    diff3 = []
    diff5 = []

    for i in range(len(df['harm1']) - 1):
        diff1.append(df['harm1'][i+1] - df['harm1'][i])
        diff3.append(df['harm3'][i+1] - df['harm3'][i])
        diff5.append(df['harm5'][i+1] - df['harm5'][i])
    
    print(f'For the 1st Harmonic the average of difference of energy for 1 micron is {round(numpy.average(diff1), 5)} eV')
    print(f'For the 3rd Harmonic the average of difference of energy for 1 micron is {round(numpy.average(diff3), 5)} eV')
    print(f'For the 5th Harmonic the average of difference of energy for 1 micron is {round(numpy.average(diff5), 5)} eV')

    plt.ion()
    plt.figure()    

    if plot_harmonic == 1:
        plt.plot(df['ID_gap'], df['harm1'], 'bo-', label='1$^{st}$ Harmonic')
    elif plot_harmonic == 3:
        plt.plot(df['ID_gap'], df['harm3'], 'go-', label='3$^{rd}$ Harmonic')
    elif plot_harmonic == 5:
        plt.plot(df['ID_gap'], df['harm5'], 'ko-', label='5$^{th}$ Harmonic')
    else:
        raise RuntimeError("ERROR: This harmonic has not been yet considered")

    plt.xlabel("ID gap [mm]",fontsize=f_size)
    plt.ylabel("Resonance energy (eV)",fontsize=f_size)

    plt.legend(fontsize=f_size)
    #plt.title("On-energy", fontsize=f_size)

    #plt.yscale('log')
    #plt.xscale('log')
    #plt.ylim(0.1, 10)
    plt.grid(which='both', axis='both')
    plt.show()

if __name__ == "__main__":
    pass

