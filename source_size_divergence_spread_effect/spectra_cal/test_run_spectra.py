import spectra.python.spectra as spectra
import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy
import scipy.constants as codata
from syned.util.json_tools import load_from_json_file
import h5py

def e_to_k(photon_energy, syned_json_file, harmonic=1):
    
    m2ev = codata.c * codata.h / codata.e
    syned_obj = load_from_json_file(syned_json_file)
    e = syned_obj.get_electron_beam()
    u = syned_obj.get_magnetic_structure()
    
    wavelength = harmonic * m2ev / photon_energy
    
    K = round(numpy.sqrt(2*(((wavelength*2*e.gamma()**2)/u.period_length())-1)), 6)
    
    print(f'Energy {photon_energy} corresponds a K value of {K}')
    
    return K

def run_spectra(photon_energy, syned_json_file, harmonic=1, plot=False):

    #for now it works only for 
    f_size = 12
    # open an input file
    f = open("C://Work_JR//spectra_win//config_set_undulator_at_point_source_spatial_profile.json")
    prm = json.load(f)
    
    syned_obj = load_from_json_file(syned_json_file)
    u = syned_obj.get_magnetic_structure()
    und_name = syned_obj.get_name()
    #change the undulator magnetic period lenght
    prm["Light Source"]["&lambda;<sub>u</sub> (mm)"] = u.period_length() * 1e3 #mm
    prm["Light Source"]["Device Length (m)"] =  u.length()
    #prm["Light Source"]["K value"] = u.K_vertical()
    prm["Light Source"]["K value"] = e_to_k(photon_energy, syned_json_file, harmonic=harmonic)    
    prm["Configurations"]["Target Harmonic"] = harmonic
    
    prmstr = json.dumps(prm)
    
    # call solver with the input string (JSON format)
    solver = spectra.Solver(prmstr)
    
    # check if the paramete load is OK
    isready = solver.IsReady()
    if isready == False:
        print("Parameter load failed.")
        sys.exit()
    
    # start calculation
    solver.Run()
    
    data = solver.GetData()
    xyarray = data["variables"]
    flux_dens = data["data"]
        
    x = xyarray[0]
    y = xyarray[1]
        
    nx = len(x)
    ny = len(y)
    zdata = numpy.array(flux_dens).reshape([ny, nx])
    
    if plot:
        # get results
        # get titles and units
        captions = solver.GetCaptions()
        
        titles = captions["titles"]
        units = captions["units"]
               
        #plt.ion()
        #plt.figure()
        plt.pcolormesh(x, y, zdata, cmap=plt.cm.viridis, shading='auto')	
        plt.colorbar().ax.tick_params(axis='y', labelsize=f_size)
        plt.title("{} Photon Energy {} eV".format(\
                    und_name, photon_energy), fontsize=f_size)
        
        plt.xlabel(titles[0], fontsize=f_size)
        plt.ylabel(titles[1], fontsize=f_size)
        plt.xticks(fontsize=f_size)
        plt.yticks(fontsize=f_size)
        
        plt.show()
        
    return x, y, zdata
    
def test_save_h5(photon_energies, syned_json_file, harmonic=1, spectra_data= 'charac_at_source_point', save_file=True):

        
    syned_obj = load_from_json_file(syned_json_file)
    
    flux_densities = []
    
    for photon_energy in photon_energies:
        
        x, y, zdata = run_spectra(photon_energy, syned_json_file, harmonic=harmonic)
        flux_densities.append(zdata)
        #tz_data = numpy.array(zdata).T
        #flux_densities.append(tz_data)
    
    #flux_densities = numpy.array(flux_densities)
    #array_energies = numpy.array(photon_energies)
    #array_energies = array_energies.astype(numpy.float)   
    
    if save_file:
    
        HDF5_FILE = 'test_spectra_{}_{}.hdf5'.format(spectra_data, syned_obj.get_name())
        
        f = h5py.File(HDF5_FILE, 'a')
        
        nxentry = f.create_group('Scan')
        nxentry.attrs['NX_class'] = 'NXentry'
        
        nxdata = nxentry.create_group('harmonic {}'.format(harmonic))
        nxdata.attrs['NX_class'] = 'NXdata'
        nxdata.attrs['signal'] = '%s'%("flux_dens")
        nxdata.attrs['axes'] = ['energy', 'y', 'x']
            
        # Image data
        fs = nxdata.create_dataset("flux_dens", data=flux_densities)
        fs.attrs['long_name'] = 'Flux density'
    
        es = nxdata.create_dataset('energy', data=photon_energies)
        es.attrs['units'] = 'eV'
        es.attrs['long_name'] = 'Photon Energy (eV)'    
    
        # X axis data
        xs = nxdata.create_dataset('x', data=x)
        xs.attrs['units'] = 'mm'
        xs.attrs['long_name'] = 'horizontal (mm)'    
    
        # Y axis data
        ys = nxdata.create_dataset('y', data=y)
        ys.attrs['units'] = 'mm'
        ys.attrs['long_name'] = 'vertical (mm)'    
    
        f.close()
        
        print ('File Spectra_{}_{}.hdf5 save to the disk'.format(spectra_data, syned_obj.get_name())) 
    else:
        pass
    
    
    
def save_full_h5(ref_photon_energy_file, syned_json_file, spectra_data= 'charac_at_source_point'):

    df_e_ref = pd.read_csv(ref_photon_energy_file, sep='\t', comment = '#', engine='python')
    
    syned_obj = load_from_json_file(syned_json_file)
    
    #get the energies
    list_energy = []
    for i in range(int(df_e_ref.shape[1]/2)):
        list_energy.append(df_e_ref.iloc[:, i*2])    
    
    HDF5_FILE = 'Spectra_{}_{}_full.hdf5'.format(spectra_data, syned_obj.get_name())
        
    f = h5py.File(HDF5_FILE, 'a')
        
    nxentry = f.create_group('Scan')
    nxentry.attrs['NX_class'] = 'NXentry'
    
    
    for i in range(len(list_energy)):
        flux_densities = []
        for photon_energy in list_energy[i]:
            print('Calculating SPECTRA for photon energy {} and harmonic {}'.format(photon_energy, "%02d"%(i*2+1)))
            x, y, zdata = run_spectra(photon_energy, syned_json_file, harmonic=((i*2)+1))
            flux_densities.append(zdata)
        #dump in the h5file
        print('Done for the harmonic {} dumping in the h5 file...'.format("%02d"%(i*2+1)))
        
        nxdata = nxentry.create_group('harmonic {}'.format("%02d"%(i*2+1)))
        nxdata.attrs['NX_class'] = 'NXdata'
        nxdata.attrs['signal'] = '%s'%("flux_dens")
        nxdata.attrs['axes'] = ['energy', 'y', 'x']
            
        # Image data
        fs = nxdata.create_dataset("flux_dens", data=flux_densities)
        fs.attrs['long_name'] = 'Flux density'
    
        es = nxdata.create_dataset('energy', data=list_energy[i])
        es.attrs['units'] = 'eV'
        es.attrs['long_name'] = 'Photon Energy (eV)'    
    
        # X axis data
        xs = nxdata.create_dataset('x', data=x)
        xs.attrs['units'] = 'mm'
        xs.attrs['long_name'] = 'horizontal (mm)'    
    
        # Y axis data
        ys = nxdata.create_dataset('y', data=y)
        ys.attrs['units'] = 'mm'
        ys.attrs['long_name'] = 'vertical (mm)'
    
    f.close()
    
    print ('File {} save to the disk'.format(HDF5_FILE))   
    
if __name__=='__main__':
    #test_save_h5([9.07126e+3, 1.86203e+4], "ESRF_ID06_EBS_CPMU18_1.json", harmonic=1, save_file=True)
    save_full_h5('spectra_results/SPECTRA_EBS_cpmu18_hor_size.txt', 'ESRF_ID06_EBS_CPMU18_1.json', spectra_data= 'charac_at_source_point')
