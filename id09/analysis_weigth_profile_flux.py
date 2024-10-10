import h5py
import pandas as pd
import numpy
import matplotlib.pyplot as plt
from srxraylib.util.h5_simple_writer import H5SimpleWriter
import scipy.constants as codata


def get_flux_spectrum(spectrum_file, new_energies):

    """ Function to interpolate the new energies and obtain the corresponding flux """
    
    df_spectrum = pd.read_csv(spectrum_file, sep=',', comment = '#', engine='python') 
    
    if 'power' in spectrum_file:

        energy = df_spectrum['Photon Energy [eV]']
        power = df_spectrum['Spectral power [W/eV]']
        flux = power / (codata.e * 1e3)
    
    else:
        energy = df_spectrum['Energy [eV]']
        flux = df_spectrum['Flux [Phot/sec/0.1%bw]']


    new_flux = numpy.interp(new_energies, energy, flux)

    return new_flux

def weigth_profile(file, spectrum_file, save_file = False): 

    """ Get all the profiels for the different energies and weight them using 
      the photon flux  """   

    h5_file = h5py.File(file, 'r')
    
    profiles = numpy.copy(h5_file['Run/Profiles/beam'])    
    photon_energies = numpy.copy(h5_file['Run/Profiles/photon_energy'])    
    x = numpy.copy(h5_file['Run/Profiles/x'])
    y = numpy.copy(h5_file['Run/Profiles/y'])

    h5_file.close()
    
    weigths = get_flux_spectrum(spectrum_file, photon_energies) / max(get_flux_spectrum(spectrum_file, photon_energies))
    
    weigth_profile = numpy.zeros_like(profiles[0, :, :])

    for i in range(len(weigths)):
        weigth_profile += profiles[i, :, :] * weigths[i]

    if save_file:
        h5w = H5SimpleWriter.initialize_file(f"weight_profile_peak.h5",creator="h5_basic_writer.py")

        h5w.add_image(weigth_profile.T, x, y,
                      image_name="weight_profile",
                      title_x="h [um]",title_y="v [um]")
        print('weight_profile_peak.h5 has been save to disk')
    
    return x, y, weigth_profile

def flux_at_element(h5files, flux_labels, full_spectrum, plot=True, f_size=12):

    """ Get the flux at a given element using the transmitivity data and the photon
     flux """
     
    photon_energies = []
    fluxes_in_element = []
    total_flux = []    

    for i, h5file in enumerate(h5files):

        with h5py.File(h5file, 'r') as h5_file:
    
            photon_energies.append(numpy.copy(h5_file['Run/Transmitivity/photon_energy']))

            flux_in_element = get_flux_spectrum(full_spectrum, photon_energies[i]) *\
                              numpy.copy(h5_file['Run/Transmitivity/transmitivity'])
            
            total_flux.append(numpy.trapz(flux_in_element*1e3 / photon_energies[i],
                                          x=photon_energies[i], axis=-1))
            
            print("For the element {}: Total flux of {:0.2e} Photons/s".format(flux_labels[i], total_flux[i]) )

            fluxes_in_element.append(flux_in_element)
             
        
    if plot:
        for i, photon_energy in enumerate(photon_energies):
            plt.plot(photon_energy, fluxes_in_element[i], label=flux_labels[i])
        plt.xticks(fontsize=f_size)
        plt.yticks(fontsize=f_size)        
        plt.xlabel("Photon energy [eV]", fontsize=f_size)
        plt.ylabel("Flux [Photon/s/0.1%bw]", fontsize=f_size)
        plt.legend(fontsize=f_size)
        plt.show()


def plot_profile(x, y, profile, plot=True):

    #TODO add a fitting rutine to get directly the FWHMs
    # For now we can obtain them by opening the weight_profile_peak.h5 with the
    # ASYS HDF5 reader 

    if plot:
        plt.pcolormesh(x,y,profile, cmap=plt.cm.viridis)
        plt.colorbar()
        plt.ylabel("Vertical [um]",fontsize=12)
        plt.xlabel("Horizontal [um]",fontsize=12)
        plt.show()

if __name__=="__main__":
    #pass

    
    #x, y, profile = weigth_profile('results/peak_run.hdf5', 'u17_spectrum_15keV_peak.csv')    

    #flux_at_element(['results/opened_chopper1.hdf5','results/closed_chopper.hdf5'],
    #                ['Full open chopper','Actual aperture chopper'], 'u17_power_spectrum_15keV_peak.csv')

    flux_at_element(['results/sample_ideal_lens.hdf5', 'results/sample_50_200um_diamond.hdf5',
                    'results/sample_12_50um_diamond.hdf5'], ['Ideal lens',
                    '50_Diamond_lens_r=200 $\mu$m', '12_Diamond_lens_r=50 $\mu$m'],
                    'u17_power_spectrum_15keV_peak.csv')

    #x, y, profile = weigth_profile('results/opened_chopper1.hdf5', 'u17_spectrum_15keV_peak.csv', save_file=True)

    #x, y, profile = weigth_profile('results/sample_ideal_lens.hdf5', 'u17_spectrum_15keV_peak.csv', save_file=True)

    #x, y, profile = weigth_profile('results/sample_50_200um_diamond.hdf5', 'u17_spectrum_15keV_peak.csv', save_file=True)

    #x, y, profile = weigth_profile('results/sample_12_50um_diamond.hdf5', 'u17_spectrum_15keV_peak.csv', save_file=True)

    #plot_profile(x, y, profile)



    




 
