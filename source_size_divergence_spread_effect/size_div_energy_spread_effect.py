import numpy
from scipy.special import erf
import scipy.constants as codata
from syned.util.json_tools import load_from_json_file
import matplotlib.pyplot as plt
import pandas as pd
from srxraylib.plot.gol import plot as manolo_plot
from silx.math.fit import FitManager
from silx.math.fit import fittheories
from silx.math.fit.functions import sum_lorentz, sum_gauss
import h5py
import re


def ev_to_m(photon_energy):
    """ very short function just to get the wavelength in meters from the photon
    energy in eV using the formula: lambda = hc/E"""
    # scipy plack constant is in J/s, which is transformed to eV/s 
    return (codata.h / codata.electron_volt * codata.c / photon_energy)

def delta_res_energy(delta, json_file, harmonic=1):
    
    syned_obj = load_from_json_file(json_file)
    e = syned_obj.get_electron_beam()
    u = syned_obj.get_magnetic_structure()

    electron_energy = e.energy()
    gamma = e.gamma()
    resonance_wavelength =  ((u.period_length() / (2.0 * gamma **2)) * (1 + u.K()**2 / 2.0))/harmonic
    print(f'Resonance_wavelenght {resonance_wavelength}')
    resonance_energy =  codata.c * codata.h / codata.e / resonance_wavelength

    new_res_energy = (delta * resonance_energy) / (u.number_of_periods() * harmonic) + resonance_energy   

    print(f"Resonant energy is: {resonance_energy}")
    print(f"New resonant energy is: {new_res_energy}")

    return new_res_energy

def photon_source(photon_energy, syned_json_file, aprox='fit'):
    """ Function to get the radiation beam size for a given energy
    and undulator length for Onuki&Elleaume eq.25,30 and ('fit') or Gaussian aprox ('Gauss') """

    syned_obj = load_from_json_file(syned_json_file)
    u = syned_obj.get_magnetic_structure()
    u_length = u.period_length()*u.number_of_periods()   
    
    if aprox == 'fit':
        #see formulas 25 & 30 in Elleaume (Onuki & Elleaume)        
        photon_size = 2.740 / (4e0 * numpy.pi) * numpy.sqrt(ev_to_m(photon_energy) * u_length)
        photon_div  = 0.69 * numpy.sqrt(ev_to_m(photon_energy) /u_length)
                                    
    elif  aprox == 'Gauss':
        photon_size = (1/(4 * numpy.pi)) * numpy.sqrt(2 * ev_to_m(photon_energy) * u_length)
        photon_div  = numpy.sqrt(ev_to_m(photon_energy) / (2 * u_length))
        
    else:
        raise RuntimeError("ERROR: Unidentified aproximation")

    return photon_size, photon_div

def norm_energ_spr(energy_spread, json_file, harmonic=1):
    """ Tanaka & Kitamura 2009 Normalized energy spread
    equation (13)""" 

    syned_obj = load_from_json_file(json_file)
    u = syned_obj.get_magnetic_structure()
    return 2 * numpy.pi * harmonic * u.number_of_periods() * energy_spread

def q_a(x):
    """ Tanaka & Kitamura 2009 equation (17), hard modified to give
    1 in the limit close to zero"""

    if x > 1e-5:
        f_1 = -1 + numpy.exp(-2*numpy.square(x)) + numpy.sqrt(2*numpy.pi) * x * erf(numpy.sqrt(2)*x)
        value = numpy.sqrt(2*numpy.square(x)/f_1)
    elif x < 1e-5 and x >= 0:
        value = 1.0
    else:
        RuntimeError('ERROR: Please provide a valid input')
    return value

def q_s(x, factor=1):
    """ Tanaka & Kitamura 2009 equation (24)"""    
    return 2 * numpy.power(q_a(x/4), 2/3) * factor

def source_size_typical(photon_energy, syned_json_file, aprox='fit', emittance=True):
    """ Typical function to get the source size by convolution of electron beam
     and photon source can use Gaussiand or 'fit' Onuki&Elleaume"""
    if emittance:
        syned_obj = load_from_json_file(syned_json_file)
        e = syned_obj.get_electron_beam()
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = e.get_sigmas_all()
    else:
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = (0.0, 0.0, 0.0, 0.0)

    (photon_size, photon_div) = photon_source(photon_energy, syned_json_file, aprox=aprox)    

    source_h = numpy.sqrt( sigma_x**2 + photon_size**2 )
    source_v = numpy.sqrt( sigma_y**2 + photon_size**2 )
    source_hp = numpy.sqrt(sig_div_x**2 + photon_div**2 )
    source_vp = numpy.sqrt(sig_div_y**2 + photon_div**2 )
        
    return source_h, source_v, source_hp, source_vp

def source_size_e_spread(photon_energy, syned_json_file, energy_spread, aprox='Gauss', harmonic=1, emittance=True):
    """ Get source size by convolution with correction Tanaka&Kitamura(2009)
     equations 27 and 28 """
    if emittance:
        syned_obj = load_from_json_file(syned_json_file)
        e = syned_obj.get_electron_beam()
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = e.get_sigmas_all()
    else:
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = (0.0, 0.0, 0.0, 0.0)

    #Here we calculate the photon size and div for a given Gauss or 'fit' Onuki & Elleaume
    (photon_size, photon_div) = photon_source(photon_energy, syned_json_file, aprox=aprox)  
    #get the normilized energy spread
    nor_e_spr = norm_energ_spr(energy_spread, syned_json_file, harmonic=harmonic)

    source_h = numpy.sqrt( sigma_x**2 + (photon_size**2 * q_s(nor_e_spr)**2))
    source_v = numpy.sqrt( sigma_y**2 + (photon_size**2 * q_s(nor_e_spr)**2))
    source_hp = numpy.sqrt(sig_div_x**2 + (photon_div**2 * q_a(nor_e_spr)**2))
    source_vp = numpy.sqrt(sig_div_y**2 + (photon_div**2 * q_a(nor_e_spr)**2))
        
    return source_h, source_v, source_hp, source_vp

def diff_res_energy(electron_energy, json_file, harmonic=1):

    #electron energy in GeV
    
    syned_obj = load_from_json_file(json_file)
    e = syned_obj.get_electron_beam()
    u = syned_obj.get_magnetic_structure()

    electron_energy = electron_energy
    gamma = 1e9 * electron_energy / (codata.m_e *  codata.c**2 / codata.e)
    resonance_wavelength =  ((u.period_length() / (2.0 * gamma **2)) * (1 + u.K()**2 / 2.0))/harmonic
    print(f'Resonance_wavelenght {resonance_wavelength}')
    resonance_energy =  codata.c * codata.h / codata.e / resonance_wavelength
    print(f'For {electron_energy} GeV the resonance energy is {resonance_energy}')    

    return resonance_energy


def source_size_e_spread_det(json_file, energy_spread, sigma_units=6.0, num_points=11, aprox='fit', harmonic=1):

    syned_obj = load_from_json_file(json_file)
    e = syned_obj.get_electron_beam()
    electron_energy_center = e.energy()
    (sigma_x, sig_div_x, sigma_y, sig_div_y) = e.get_sigmas_all()

    sigma_delta = energy_spread * electron_energy_center
    interval = sigma_units * sigma_delta
    electron_energies = numpy.linspace(electron_energy_center - 0.5 * interval,
                                      electron_energy_center + 0.5 * interval,
                                      num_points)
    weights = numpy.exp(-(electron_energies - electron_energy_center) ** 2 / 2 / sigma_delta ** 2)
    weights /= weights.sum()
    print("electron energies: ", electron_energies)
    print("weigths: ", weights)

    for i, electron_energy_i in enumerate(electron_energies):
        electron_energy = electron_energy_i

        res_photon_energy = diff_res_energy(electron_energy, json_file, harmonic=harmonic)
        (photon_size, photon_div) = photon_source(res_photon_energy, json_file, aprox=aprox)
        
        source_h = numpy.sqrt( sigma_x**2 + photon_size**2 )
        source_v = numpy.sqrt( sigma_y**2 + photon_size**2 )
        source_hp = numpy.sqrt(sig_div_x**2 + photon_div**2 )
        source_vp = numpy.sqrt(sig_div_y**2 + photon_div**2 )

        if i == 0:
            source_h *= weights[i]
            source_v *= weights[i] 
            source_hp *= weights[i]
            source_vp *= weights[i]
        else:
            source_h += source_h * weights[i]
            source_v += source_v * weights[i] 
            source_hp += source_hp * weights[i]
            source_vp += source_vp * weights[i]
        
    return source_h, source_v, source_hp, source_vp

def run_cal_tanaka(syned_json_file, energy_spread, source_dim='hor_div', save_file=True, emittance=True):

    if source_dim == 'hor_size':
        df_spectra = pd.read_csv('spectra_results/SPECTRA_EBS_cpmu18_hor_size.txt', sep='\t', comment = '#', engine='python')
        indx_cal = 0

    elif source_dim == 'ver_size':
        df_spectra = pd.read_csv('spectra_results/SPECTRA_EBS_cpmu18_ver_size.txt', sep='\t', comment = '#', engine='python')
        indx_cal = 1
    
    elif source_dim == 'hor_div':
        df_spectra = pd.read_csv('spectra_results/SPECTRA_EBS_cpmu18_hor_div.txt', sep='\t', comment = '#', engine='python')
        indx_cal = 2

    elif source_dim == 'ver_div':
        df_spectra = pd.read_csv('spectra_results/SPECTRA_EBS_cpmu18_ver_div.txt', sep='\t', comment = '#', engine='python')
        indx_cal = 3               

    else:
        raise RuntimeError("ERROR: Unidentified source dimension")
    
    list_energy = []
    for i in range(int(df_spectra.shape[1]/2)):
        list_energy.append(df_spectra.iloc[:, i*2])
    harm_energies = numpy.array(list_energy)
    list_dim = []

   
    name_output_file = f'Tanaka_{source_dim}'        
    for i in range(len(list_energy)):
        dim = []
        for energy in list_energy[i]:
            dim.append(source_size_e_spread(energy, syned_json_file, energy_spread, harmonic=((i*2)+1), emittance=emittance)[indx_cal])
        list_dim.append(dim)    

    harm_dim = numpy.array(list_dim)

    if save_file:
        for i in range(len(list_energy)):
            if i == 0:
                df_tanaka = pd.DataFrame({f'Photon Energy {(i*2)+1}':list_energy[i], f'{source_dim} {(i*2)+1}':list_dim[i]})
            else:
                df_tanaka[f'Photon Energy {(i*2)+1}'] = list_energy[i]
                df_tanaka[f'{source_dim} {(i*2)+1}'] = list_dim[i]
        
        if emittance:
            df_tanaka.to_csv(f'{name_output_file}_spread_{energy_spread}.csv', index=False)
            print(f'file {name_output_file}_spread_{energy_spread}.csv has been saved to disk')
        else:
            df_tanaka.to_csv(f'{name_output_file}_spread_{energy_spread}_zero_emit.csv', index=False)
            print(f'file {name_output_file}_spread_{energy_spread}_zero_emit.csv has been saved to disk')        
    
    return harm_energies, harm_dim

def run_cal_shadow(syned_json_file, source_dim='hor_div', save_file=True, emittance=True):

    if source_dim == 'hor_size':
        df_spectra = pd.read_csv('spectra_results/SPECTRA_EBS_cpmu18_hor_size.txt', sep='\t', comment = '#', engine='python')
        indx_cal = 0

    elif source_dim == 'ver_size':
        df_spectra = pd.read_csv('spectra_results/SPECTRA_EBS_cpmu18_ver_size.txt', sep='\t', comment = '#', engine='python')
        indx_cal = 1
    
    elif source_dim == 'hor_div':
        df_spectra = pd.read_csv('spectra_results/SPECTRA_EBS_cpmu18_hor_div.txt', sep='\t', comment = '#', engine='python')
        dimension = source_dim
        indx_cal = 2

    elif source_dim == 'ver_div':
        df_spectra = pd.read_csv('spectra_results/SPECTRA_EBS_cpmu18_ver_div.txt', sep='\t', comment = '#', engine='python')
        indx_cal = 3             

    else:
        raise RuntimeError("ERROR: Unidentified source dimension")
    
    list_energy = []
    for i in range(int(df_spectra.shape[1]/2)):
        list_energy.append(df_spectra.iloc[:, i*2])
    harm_energies = numpy.array(list_energy)
    list_dim = []        
    for i in range(len(list_energy)):
        dim = []
        for energy in list_energy[i]:
            dim.append(source_size_typical(energy, syned_json_file, aprox='fit', emittance=emittance)[indx_cal])
        list_dim.append(dim)
    harm_dim = numpy.array(list_dim)

    if save_file:
        for i in range(len(list_energy)):
            if i == 0:
                df_shadow = pd.DataFrame({f'Photon Energy {(i*2)+1}':list_energy[i], f'{source_dim} {(i*2)+1}':list_dim[i]})
            else:
                df_shadow[f'Photon Energy {(i*2)+1}'] = list_energy[i]
                df_shadow[f'{source_dim} {(i*2)+1}'] = list_dim[i]
        if emittance:
            df_shadow.to_csv(f'Shadow_{source_dim}.csv', index=False)
            print(f'file Shadow_{source_dim}.csv has been saved to disk')
        else:
            df_shadow.to_csv(f'Shadow_{source_dim}_zero_emit.csv', index=False)
            print(f'file Shadow_{source_dim}_zero_emit.csv has been saved to disk')
            
    return harm_energies, harm_dim

def plot_comparison(files, variable_name = 'Horizontal divergence', units='urad'):
    """ Plot comparison between the different methods, files is a list """
    f_size = 12
    for file_csv in files:
        if "Shadow" in file_csv:
            df = pd.read_csv(file_csv, sep=',', comment = '#', engine='python')
            label = 'Shadow'
            
        elif "Tanaka" in file_csv:
            df = pd.read_csv(file_csv, sep=',', comment = '#', engine='python')
            label = 'Tanaka & Kitamura (2009)'

        elif "spectra" in file_csv:
            df = pd.read_csv(file_csv,  sep='\t', comment = '#', engine='python')
            label='SPECTRA: Gaussian approx'

        elif "sirepo" in file_csv:
            df = pd.read_csv(file_csv,  sep='\t', comment = '#', engine='python')
            label = 'Sirepo'

        else:
           raise RuntimeError("ERROR: Unidentified calculation approach - plot aborted")
        
        photon_energy = []
        source_variable = []
        for i in range(int(df.shape[1]/2)):
            if "sirepo" in file_csv:
                photon_energy += [j * 1e3 for j in df.iloc[:, i*2].to_list()] + [numpy.nan]
            else:
                photon_energy += df.iloc[:, i*2].to_list() + [numpy.nan]

            source_variable += [j * 1e6 for j in df.iloc[:, (i*2)+1].to_list()] + [numpy.nan]
        plt.plot(photon_energy, source_variable, label=label)
    
    plt.xlabel("Photon energy [eV]", fontsize= f_size)    
    plt.ylabel(f"{variable_name} ({units})", fontsize= f_size)

    plt.xlim(4e3, 21e4)
    #plt.ylim(1e12, 5e15)
    plt.xticks(fontsize= f_size)
    plt.yticks(fontsize= f_size)
    #plt.yscale("log")
    plt.grid(which='both', axis='y')
    plt.title("EBS ID06 CMPU18, On Resonance, $\epsilon$=0, $\delta$=0", fontsize=f_size)     

    plt.legend(fontsize = f_size)
       
    plt.show() 
            
def get_sigmas_from_spectra_wigner(h5file, source_dim='hor_dim', fit_meth='mid', plot_fit=False):
    """Function to get the RMS from fitting aprox a Gaussian to the central
    image profile in vertical and horizontal line from """

    #read the h5file
    f = h5py.File(h5file, 'r')
    
    try:
        scan = f['Scan']               
    except:
        print('Failing on getting the harmonics on Scan group')
    
    for i, key in enumerate(scan.keys()):
                
        harmonic = scan.get(key)
        flux_dens = numpy.array(harmonic['flux_dens'])
        photon_energies = numpy.array(harmonic['energy'])
        x = numpy.array(harmonic['x'])
        y = numpy.array(harmonic['y'])  
        harm_sigmas = []                 
        for j, photon_energy in enumerate(photon_energies):
            e_flux_dens = flux_dens[j, :, :]
            #get the central profiles
            if fit_meth == 'max':
                max_indx = numpy.unravel_index(e_flux_dens.argmax(), e_flux_dens.shape)
                max_indx_hor = max_indx[0]
                max_indx_ver = max_indx[0]

            elif fit_meth == 'mid':
                max_indx_hor = int(len(x)/2)
                max_indx_ver = int(len(y)/2)
            else:
                RuntimeError("ERROR: Method to get the central profile is not implemented, options: 'max' or 'mid'")

            if 'hor' in source_dim:
                axe = x + 1e-8 #for some reason it was falling all the time for the fitting, by modifying the number it worked
                f_axe = e_flux_dens[max_indx_hor,:]
                label_axe = 'horizontal'
            elif 'ver' in source_dim:
                axe = y + 1e-8
                f_axe = e_flux_dens[:, max_indx_ver]
                label_axe = 'vertical'
            else:
                raise RuntimeError("ERROR: Unidentified dim to analize")
            
            #Fitting process for each slit position            
            try:    
                fit = FitManager()
                fit.setdata(x=axe, y=f_axe)
                fit.loadtheories(fittheories)
                fit.settheory('Gaussians') #'Gaussians', 'Lorentz'
                fit.setbackground('No Background')
                fit.estimate()
                fit.runfit()    
                height, pp, fwhm = (param['fitresult'] for param in fit.fit_results)
                fit_axe = sum_gauss(numpy.array(axe), *[height, pp, fwhm])
                harm_sigmas.append(fwhm * 1e3 / (2 * numpy.sqrt(2 * numpy.log(2))))
                if plot_fit:
                    #Warning plot fit for every pow dens 2D map
                    manolo_plot(axe, f_axe, axe, fit_axe, legend=['{} profile'.format(label_axe), 'Gaussian fit'],
                                title='{}, photon energy:{}, FWHM:{}'.format(key, photon_energy, fwhm))
            except:
                print(f'Fitting failed for {key}:{photon_energy} indx={j}, trying fot cutted profile around max peak to fit')
                range_cut = 100
                print(f'Here is x:{axe} and here is the type {type(axe)}')
                print(f'Here is y:{f_axe} and here is the type {type(axe)}')
                #define a range around the peak:
                peak_ind = numpy.where(f_axe == numpy.amax(f_axe))[0][0]
                #print("peak index: ", peak_ind)        
                cut_axe = axe[peak_ind - int(range_cut/2) : peak_ind + int(range_cut/2)]
                cut_f_axe = f_axe[peak_ind - int(range_cut/2) : peak_ind + int(range_cut/2)]
                manolo_plot(cut_axe, cut_f_axe, legend=['Cutted {} profile'.format(label_axe)],
                                title='{}, photon energy:{}, FWHM:{}'.format(key, photon_energy, fwhm))               
                #Fitting to the cutted profile
                fit = FitManager()
                fit.setdata(x=cut_axe, y=cut_f_axe)
                fit.loadtheories(fittheories)
                fit.settheory('Gaussians') #'Gaussians', 'Lorentz'
                fit.setbackground('No Background')
                fit.estimate()
                fit.runfit()    
                height, pp, fwhm = (param['fitresult'] for param in fit.fit_results)
                cut_fit_axe = sum_gauss(numpy.array(cut_axe), *[height, pp, fwhm])
                harm_sigmas.append(fwhm * 1e3 / (2 * numpy.sqrt(2 * numpy.log(2))))
                if plot_fit:
                    #Warning plot fit for every pow dens 2D map
                    manolo_plot(cut_axe, cut_f_axe, cut_axe, cut_fit_axe,
                                legend=['Cutted {} profile'.format(label_axe), 'Gaussian fit'],
                                title='Photon energy {} eV with a FWHM of {}'.format(photon_energy, fwhm))            

        if i == 0:
            df_wigner = pd.DataFrame({f'{key}':photon_energies, f'{source_dim}_1':harm_sigmas})
        else:
            
            df_wigner[f'{key}'] = photon_energies
            df_wigner['{}_{}'.format(source_dim, re.findall(r'\d+', key)[0])] = harm_sigmas
        print('Done with {} going to the next...'.format(key))    
            
    df_wigner.to_csv(f'Spectra_Wigner_{source_dim}.csv', index=False)
    print(f'file Spectra_Wigner_{source_dim}.csv has been saved to disk')        
        
    f.close()

    #return df_wigner

def test_fit(h5file, harm_num=1, e_indx=0):
    #read the h5file
    f = h5py.File(h5file, 'r')
    
    try:
        data = f[f'Scan/harmonic 0{harm_num}']               
    except:
        print('Failing on getting the harmonics on Scan group')

    flux_dens = numpy.array(data['flux_dens'])
    photon_energies = numpy.array(data['energy'])    
    x = numpy.array(data['x'])
    f_e = numpy.copy(flux_dens[e_indx,:,:])

    print(f_e)

    max_indx_hor = int(len(x)/2)

    f_x = f_e[max_indx_hor,:]

    print(type(f_x))

    manolo_plot(x, f_x)
    
    try:
        print(f'This fit is for Energy:{photon_energies[e_indx]} taking the line with indx{max_indx_hor}')
            
        fit = FitManager()   
        fit.setdata(x=numpy.round(x,3), y=f_x)
        fit.loadtheories(fittheories)
        fit.settheory('Gaussians') #'Gaussians', 'Lorentz'
        fit.setbackground('No Background')
        fit.estimate()
        fit.runfit()    
        height, pp, fwhm = (param['fitresult'] for param in fit.fit_results)
        fit_y = sum_gauss(x, *[height, pp, fwhm])    
        
        manolo_plot(x, f_x, x, fit_y, legend=['Original', 'Gaussian fit'],
                    title='Photon energy {} eV with a FWHM of {}'.format(photon_energies[e_indx], fwhm))
    except:
        print('fitting failed')   
    f.close()
    return x, f_x 

def test_fitting(x, y):        
    #print(f'Here is x:{x} and here is the type {type(x)}')
    #print(f'Here is y:{y} and here is the type {type(y)}')
    fit = FitManager()
    fit.setdata(x=x, y=y)
    fit.loadtheories(fittheories)
    fit.settheory('Gaussians') #'Gaussians', 'Lorentz'
    fit.setbackground('No Background')
    fit.estimate()
    fit.runfit()    
    height, pp, fwhm = (param['fitresult'] for param in fit.fit_results)
    fit_y = sum_gauss(x, *[height, pp, fwhm])
    
    
    manolo_plot(x, y, x, fit_y, legend=['Original', 'Gaussian fit'])

if __name__=='__main__':
    #h_df = get_sigmas_from_spectra_wigner('spectra_results/Spectra_charac_at_source_point_ESRF_ID06_EBS_CPMU18_1_full.hdf5', source_dim='hor_dim', plot_fit=True)
    #v_df = get_sigmas_from_spectra_wigner('spectra_results/Spectra_charac_at_source_point_ESRF_ID06_EBS_CPMU18_1_full.hdf5', source_dim='ver_dim', plot_fit=True)
    #run_cal_shadow('ESRF_ID06_EBS_CPMU18_1.json', source_dim='hor_size', save_file=True, emittance=False)
    #run_cal_shadow('ESRF_ID06_EBS_CPMU18_1.json', source_dim='ver_size', save_file=True, emittance=False)

    #run_cal_tanaka('ESRF_ID06_EBS_CPMU18_1.json', 0.0, source_dim='hor_size', save_file=True, emittance=False)
    #run_cal_tanaka('ESRF_ID06_EBS_CPMU18_1.json', 0.0, source_dim='ver_size', save_file=True, emittance=False)
    #run_cal_tanaka('ESRF_ID06_EBS_CPMU18_1.json', 0.0, source_dim='hor_div', save_file=True, emittance=False)
    #run_cal_tanaka('ESRF_ID06_EBS_CPMU18_1.json', 0.0, source_dim='ver_div', save_file=True, emittance=False)
    dim = 'hor_size' #hor_div
    plot_comparison([f'Shadow_{dim}_zero_emit.csv',f'Tanaka_{dim}_spread_0.0_zero_emit.csv',f'sirepo_results/Sirepo_EBS_cpmu18_{dim}_zero_emit_zero_disp.dat', f'spectra_results/SPECTRA_EBS_cpmu18_{dim}_zero_emit_zero_disp.txt'], variable_name = 'Photon source size', units='um')