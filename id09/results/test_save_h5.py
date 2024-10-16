import matplotlib.pyplot as plt
import numpy
import h5py

#0 for the fist run, 1 for the next
ii = 0

energy_points = 10


if ii == 0:
    print("run zero: ", ii)    
    photon_energies = []
    print(photon_energies)
    data_xy = []
    
else:       

    #getting the beam

    beam = in_object_1._beam
    photon_energies.append(round(beam.getshcol(11)[0], 2))
    print(f'Energy step {len(photon_energies)}') 
    print(f'Foton energies {photon_energies}') 

    #getting the histogram and dimensions
    tkt = beam.histo2(1, 3, ref=23, xrange=[-20e-6, 20e-6], yrange=[-20e-6, 20e-6], nbins=101, nolost=1)

    x = tkt['bin_h_center']
    y = tkt['bin_v_center']
    histogram = tkt['histogram']

    data_xy.append(histogram.T)

if len(photon_energies) == energy_points:
    f = h5py.File('/home/esrf/reyesher/OASYS/ID09/new_study_2024/results/test2.hdf5', 'a')

    nxentry = f.create_group('Run')
    nxentry.attrs['NX_class'] = 'NXentry'

    nxdata = nxentry.create_group('data1')
    nxdata.attrs['NX_class'] = 'NXdata'
    nxdata.attrs['signal'] = 'beam'
    nxdata.attrs['axes'] = ['photon_energy','y','x']
    
    #image data
    bm = nxdata.create_dataset('beam', data=numpy.array(data_xy))
    bm.attrs['long_name'] = 'beam profile'

    pe = nxdata.create_dataset('photon energy', data=numpy.array(photon_energies))
    pe.attrs['units'] = 'eV'
    pe.attrs['long_name'] = 'Photon energy (eV)'
    # X axis data
    xs = nxdata.create_dataset('x', data=x)
    xs.attrs['units'] = 'm'
    xs.attrs['long_name'] = 'Horizontal (m)'
    # Y axis data
    ys = nxdata.create_dataset('y', data=y)
    xs.attrs['units'] = 'm'
    xs.attrs['long_name'] = 'Vertical (m)'

    f.close()
    
    print('File test1.hdf5 save to disk')    

else:
    pass
 
