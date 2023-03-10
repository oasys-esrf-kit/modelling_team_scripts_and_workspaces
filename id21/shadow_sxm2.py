import Shadow
import numpy
from srxraylib.util.h5_simple_writer import H5SimpleWriter
import scipy.constants as codata

def ev_to_m(photon_energy):
    """ very short function just to get the wavelength in meters from the photon
    energy in eV using the formula: lambda = hc/E"""
    # scipy plack constant is in J/s, which is transformed to eV/s 
    return (codata.h / codata.electron_volt * codata.c / photon_energy)

def photon_source(photon_energy):
    """ Function to convolve the electron beam size with the ID
    for now is only for U42 """
    #electron properties
    e_size_h = 3.043e-05
    e_size_v = 3.645e-06
    e_div_h = 4.4e-06
    e_div_v = 1.37e-06
    #u42 length
    u42_l = 1.62
    #photon source
    photon_divx = numpy.sqrt((e_div_h)**2 + \
                 (0.69 * numpy.sqrt(ev_to_m(photon_energy) / u42_l))**2)    
    photon_divy = numpy.sqrt((e_div_v)**2 + \
                 (0.69 * numpy.sqrt(ev_to_m(photon_energy) / u42_l))**2)
    photon_sx = numpy.sqrt((e_size_h)**2 + ((2.704 / (4 * numpy.pi)) \
                * numpy.sqrt(ev_to_m(photon_energy) * u42_l))**2)
    photon_sy = numpy.sqrt((e_size_v)**2 + ((2.704 / (4 * numpy.pi)) \
                * numpy.sqrt(ev_to_m(photon_energy) * u42_l))**2)

    return photon_divx, photon_divy, photon_sx, photon_sy    

def dcmadeg(photon_energy):
    """ Short function to get the crystal angle in degress in function 
    for now is only for Si(111) """
    # Si (111) D spacing [m]
    dspa = 3.135416e-10

    return numpy.around((numpy.arcsin(ev_to_m(photon_energy)/(2 * dspa)) * 180 / numpy.pi), 10)

def foc_dist(working_distance, kbm='vfm'):
    """ This function gets the focal distance in meters for a given working
    distance (m) and a KB mirror: 'vfm' or 'hfm'. For the present ID21 SXM-II
    """
    ver_p = 50.1 #distance of vfm to source
    hor_p = 50.19 #distance of hfm to source
    hfm_l = 0.06 #hfm length

    if kbm == 'vfm':
        q = hor_p - ver_p + hfm_l/2 + working_distance
    elif kbm == 'hfm':
        q = hfm_l/2 + working_distance
    else:
        raise RuntimeError("ERROR: Unidentified mirror, please enter 'vfm' or hfm'")
    
    return numpy.around(q, decimals=2)

def diff_fwhm(photon_energy, working_distance, kbm='vfm'):
    """ This function provides a diffraction contribution term
    (in meters) to include it in the simulations results """

    #both mirrors grazing angle
    m_angle = 0.006

    if kbm == 'vfm':
        mirror_active_length = 0.090
        b_i = mirror_active_length * numpy.sin(m_angle) #beam intercept by mirror (b_i)      
        n_a = (b_i/2)/foc_dist(working_distance, kbm='vfm')  #numerical aperture (n_a)
        
    elif kbm == 'hfm':
        mirror_active_length = 0.050
        b_i = mirror_active_length * numpy.sin(m_angle)     
        n_a = (b_i/2)/foc_dist(working_distance, kbm='hfm')
        
    else:
        raise RuntimeError("ERROR: Unidentified mirror, please enter 'vfm' or hfm'")

    #difracction FWHM contribution (diff) get the wavelength first in meters
    diff = 0.44 * (ev_to_m(photon_energy) / n_a)
    
    return diff        

def run_shadow(photon_energy, working_distance, n_rays=1e6):
      
    """This function run shadow for the full beamline"""
    
    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0    
    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()
    oe5 = Shadow.OE()
    oe6 = Shadow.OE()
    oe7 = Shadow.OE()
    oe8 = Shadow.OE()
    oe9 = Shadow.OE()
    oe10 = Shadow.OE()
    oe11 = Shadow.OE()    
    #
    # Define variables. See meaning of variables in: 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #
    # Source    
    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.ISTAR1 = 0
    oe0.NPOINT = n_rays
    oe0.PH1 = photon_energy - 2.5 #small bandwidth
    oe0.PH2 = photon_energy + 2.5 #small bandwidth
    oe0.SIGDIX = photon_source(photon_energy)[0]
    oe0.SIGDIZ = photon_source(photon_energy)[1]
    oe0.SIGMAX = photon_source(photon_energy)[2]
    oe0.SIGMAZ = photon_source(photon_energy)[3]
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0
    # Primary Slits
    oe1.DUMMY = 100.0
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.N_SCREEN = 1
    oe1.RX_SLIT = numpy.array([0.00175, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.RZ_SLIT = numpy.array([0.0015, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 27.2
    #M0_1
    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FHIT_C = 1
    if 2500 <= photon_energy <= 7000:
        oe2.FILE_REFL = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/Ni.dat'
        oe2.FILE_RIP = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/M0_1_Ni_shadow.dat'
    elif 7000 < photon_energy <= 10000:
        oe2.FILE_REFL = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/Pt.dat'
        oe2.FILE_RIP = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/M0_1_Pt_shadow.dat'
    else:
        raise RuntimeError(f"ERROR: Photon energy {photon_energy} eV out of range")
    oe2.FWRITE = 3
    oe2.F_G_S = 2
    oe2.F_REFLEC = 1
    oe2.F_RIPPLE = 1
    oe2.RLEN1 = 0.11
    oe2.RLEN2 = 0.11
    oe2.RWIDX1 = 0.04
    oe2.RWIDX2 = 0.04
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 89.5989295434
    oe2.T_REFLECTION = 89.5989295434
    oe2.T_SOURCE = 2.7
    #M0_2
    oe3.ALPHA = 180.0
    oe3.DUMMY = 100.0
    oe3.FHIT_C = 1
    if 2500 <= photon_energy <= 7000:
        oe3.FILE_REFL = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/Ni.dat'
        oe3.FILE_RIP = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/M0_2_Ni_shadow.dat'
    elif 7000 < photon_energy <= 10000:
        oe3.FILE_REFL = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/Pt.dat'
        oe3.FILE_RIP = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/M0_2_Pt_shadow.dat'
    else:
        raise RuntimeError(f"ERROR: Photon energy {photon_energy} eV out of range")        
    oe3.FWRITE = 3
    oe3.F_G_S = 2
    oe3.F_MOVE = 1
    oe3.F_REFLEC = 1
    oe3.F_RIPPLE = 1
    oe3.OFFY = -0.212
    oe3.RLEN1 = 0.484125
    oe3.RLEN2 = 0.484125
    oe3.RWIDX1 = 0.04
    oe3.RWIDX2 = 0.04
    oe3.T_IMAGE = 0.0
    oe3.T_INCIDENCE = 89.5989295434
    oe3.T_REFLECTION = 89.5989295434
    oe3.T_SOURCE = 0.825
    #Empty element
    oe4.ALPHA = 90.0
    oe4.DUMMY = 100.0
    oe4.FWRITE = 3
    oe4.F_REFRAC = 2
    oe4.T_IMAGE = 0.0
    oe4.T_INCIDENCE = 0.0
    oe4.T_REFLECTION = 180.0
    oe4.T_SOURCE = 0.0
    #Secondary Slits
    oe5.DUMMY = 100.0
    oe5.FWRITE = 3
    oe5.F_REFRAC = 2
    oe5.F_SCREEN = 1
    oe5.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe5.N_SCREEN = 1
    oe5.RX_SLIT = numpy.array([0.0015, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe5.RZ_SLIT = numpy.array([0.0015, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe5.T_IMAGE = 0.0
    oe5.T_INCIDENCE = 0.0
    oe5.T_REFLECTION = 180.0
    oe5.T_SOURCE = 5.475
    #DCM-1
    oe6.DUMMY = 100.0
    oe6.FHIT_C = 1
    oe6.FILE_REFL = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/Si2_15.111'
    oe6.FWRITE = 1
    oe6.F_CENTRAL = 1
    oe6.F_CRYSTAL = 1
    oe6.PHOT_CENT = photon_energy
    oe6.RLEN1 = 0.008
    oe6.RLEN2 = 0.072
    oe6.RWIDX1 = 0.005
    oe6.RWIDX2 = 0.005
    oe6.R_LAMBDA = 5000.0
    oe6.T_IMAGE = 0.0
    oe6.T_INCIDENCE = 90 - dcmadeg(photon_energy)
    oe6.T_REFLECTION = 90 - dcmadeg(photon_energy)
    oe6.T_SOURCE = 1.9
    #DCM-2
    oe7.ALPHA = 180.0
    oe7.DUMMY = 100.0
    oe7.FHIT_C = 1
    oe7.FILE_REFL = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/Si2_15.111'
    oe7.FWRITE = 1
    oe7.F_CENTRAL = 1
    oe7.F_CRYSTAL = 1
    oe7.PHOT_CENT = photon_energy
    oe7.RLEN1 = 0.008
    oe7.RLEN2 = 0.072
    oe7.RWIDX1 = 0.005
    oe7.RWIDX2 = 0.005
    oe7.R_LAMBDA = 5000.0
    oe7.T_IMAGE = 0.0
    oe7.T_INCIDENCE = 90 - dcmadeg(photon_energy)
    oe7.T_REFLECTION = 90 - dcmadeg(photon_energy)
    oe7.T_SOURCE = 0.012
    #KB Slits    
    oe8.DUMMY = 100.0
    oe8.FWRITE = 3
    oe8.F_REFRAC = 2
    oe8.F_SCREEN = 1
    oe8.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe8.N_SCREEN = 1
    oe8.RX_SLIT = numpy.array([0.0003, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe8.RZ_SLIT = numpy.array([0.00054, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe8.T_IMAGE = 0.0
    oe8.T_INCIDENCE = 0.0
    oe8.T_REFLECTION = 180.0
    oe8.T_SOURCE = 11.488
    #VFM
    oe9.DUMMY = 100.0
    oe9.FCYL = 1
    oe9.FHIT_C = 1
    oe9.FILE_REFL = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/Ni.dat'
    oe9.FILE_RIP = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/VFM_sim_0.6urad.dat'
    oe9.FMIRR = 2
    oe9.FWRITE = 1
    oe9.F_DEFAULT = 0
    oe9.F_G_S = 2
    oe9.F_REFLEC = 1
    oe9.F_RIPPLE = 1
    oe9.RLEN1 = 0.045
    oe9.RLEN2 = 0.045
    oe9.RWIDX1 = 0.015
    oe9.RWIDX2 = 0.015
    oe9.SIMAG = foc_dist(working_distance, kbm='vfm')
    oe9.SSOUR = 51.5
    oe9.THETA = 89.6562253229
    oe9.T_IMAGE = 0.0575
    oe9.T_INCIDENCE = 89.6562253229
    oe9.T_REFLECTION = 89.6562253229
    oe9.T_SOURCE = 0.5
    #HFM
    oe10.ALPHA = 90.0
    oe10.DUMMY = 100.0
    oe10.FCYL = 1
    oe10.FHIT_C = 1
    oe10.FILE_REFL = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/Ni.dat'
    oe10.FILE_RIP = b'/home/esrf/reyesher/OASYS/OasysID21/ID21_SXM-II/HFM_sim_0.6urad.dat'
    oe10.FMIRR = 2
    oe10.FWRITE = 1
    oe10.F_DEFAULT = 0
    oe10.F_G_S = 2
    oe10.F_REFLEC = 1
    oe10.F_RIPPLE = 1
    oe10.RLEN1 = 0.025
    oe10.RLEN2 = 0.025
    oe10.RWIDX1 = 0.015
    oe10.RWIDX2 = 0.015
    oe10.SIMAG = foc_dist(working_distance, kbm='hfm')
    oe10.SSOUR = 51.59
    oe10.THETA = 89.6562253229
    oe10.T_IMAGE = foc_dist(working_distance, kbm='hfm')
    oe10.T_INCIDENCE = 89.6562253229
    oe10.T_REFLECTION = 89.6562253229
    oe10.T_SOURCE = 0.0325
    #Sample
    oe11.ALPHA = 90.0
    oe11.DUMMY = 100.0
    oe11.FWRITE = 3
    oe11.F_REFRAC = 2
    oe11.T_IMAGE = 0.0
    oe11.T_INCIDENCE = 0.0
    oe11.T_REFLECTION = 180.0
    oe11.T_SOURCE = 0.0  
    
    #Run SHADOW to create the source
    
    if iwrite:
        oe0.write("start.00")
    
    beam.genSource(oe0)
    
    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")    
    
    #
    #run optical element 1
    #
    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")
    
    beam.traceOE(oe1,1)
    
    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")    
    
    #
    #run optical element 2
    #
    print("    Running optical element: %d"%(2))
    if iwrite:
        oe2.write("start.02")
    
    beam.traceOE(oe2,2)
    
    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")    
    
    #
    #run optical element 3
    #
    print("    Running optical element: %d"%(3))
    if iwrite:
        oe3.write("start.03")
    
    beam.traceOE(oe3,3)
    
    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")    
    
    #
    #run optical element 4
    #
    print("    Running optical element: %d"%(4))
    if iwrite:
        oe4.write("start.04")
    
    beam.traceOE(oe4,4)
    
    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")    
    
    #
    #run optical element 5
    #
    print("    Running optical element: %d"%(5))
    if iwrite:
        oe5.write("start.05")
    
    beam.traceOE(oe5,5)
    
    if iwrite:
        oe5.write("end.05")
        beam.write("star.05")    
    
    #
    #run optical element 6
    #
    print("    Running optical element: %d"%(6))
    if iwrite:
        oe6.write("start.06")
    
    beam.traceOE(oe6,6)
    
    if iwrite:
        oe6.write("end.06")
        beam.write("star.06")    
    
    #
    #run optical element 7
    #
    print("    Running optical element: %d"%(7))
    if iwrite:
        oe7.write("start.07")
    
    beam.traceOE(oe7,7)
    
    if iwrite:
        oe7.write("end.07")
        beam.write("star.07")    
    
    #
    #run optical element 8
    #
    print("    Running optical element: %d"%(8))
    if iwrite:
        oe8.write("start.08")
    
    beam.traceOE(oe8,8)
    
    if iwrite:
        oe8.write("end.08")
        beam.write("star.08")    
    
    #
    #run optical element 9
    #
    print("    Running optical element: %d"%(9))
    if iwrite:
        oe9.write("start.09")
    
    beam.traceOE(oe9,9)
    
    if iwrite:
        oe9.write("end.09")
        beam.write("star.09")    
    
    #
    #run optical element 10
    #
    print("    Running optical element: %d"%(10))
    if iwrite:
        oe10.write("start.10")
    
    beam.traceOE(oe10,10)
    
    if iwrite:
        oe10.write("end.10")
        beam.write("star.10")    
    
    #
    #run optical element 11
    #
    print("    Running optical element: %d"%(11))
    if iwrite:
        oe11.write("start.11")
    
    beam.traceOE(oe11,11)
    
    if iwrite:
        oe11.write("end.11")
        beam.write("star.11")

    return beam, oe0

def write_output(working_distance, photon_energies, n_rays):
    """ This function writes an h5 output file, from a given photon
    energies array, and a working distance, please notice diffraction
    contribution is taken into account to get the final FWHMs"""
    #output tempo arrays
    out = numpy.zeros((8, photon_energies.size))
    
    for i, photon_energy in enumerate(photon_energies):
    
        beam, oe0 = run_shadow(photon_energy, working_distance, n_rays=n_rays)

        histo2 = beam.histo2(1, 3, nolost=1, nbins = 201)

        out[0,i] = photon_energy
        out[1,i] = histo2['fwhm_h'] * 1e6 # horizontal fwhm microns
        out[2,i] = histo2['fwhm_v'] * 1e6 # vertical fwhm microns
        out[3,i] = beam.intensity(nolost=1)/oe0.NTOTALPOINT #Beamline transmitivity
        out[4,i] = numpy.sqrt((histo2['fwhm_h']*1e6)**2 + (diff_fwhm(photon_energy, working_distance, kbm='hfm')*1e6)**2) #considering diff factor
        out[5,i] = numpy.sqrt((histo2['fwhm_v']*1e6)**2 + (diff_fwhm(photon_energy, working_distance, kbm='vfm')*1e6)**2) #considering diff factor
        out[6,i] = beam.intensity(nolost=1) #Beamline intensity
        out[7,i] = oe0.NTOTALPOINT #Total initial rays
    
    h5w = H5SimpleWriter.initialize_file(f'shadow_sxm2_work_dist_{working_distance*1e3}_mm.h5',
                                        creator="h5_basic_writer.py")
    h5w.create_entry(f'work_dist_{working_distance*1e3}_mm', nx_default="Beam Size")

    h5w.add_dataset(out[0,:], out[1,:], entry_name=f'work_dist_{working_distance*1e3}_mm',
                    dataset_name="Beam Size Horizontal",
                    title_x='Photon Energy (eV)', title_y='FWHM_H (um)')
    h5w.add_dataset(out[0,:], out[2,:], entry_name=f'work_dist_{working_distance*1e3}_mm',
                    dataset_name="Beam Size Vertical",
                    title_x='Photon Energy (eV)', title_y='FWHM_V (um)')
    h5w.add_dataset(out[0,:], out[3,:], entry_name=f'work_dist_{working_distance*1e3}_mm',
                    dataset_name="Transmitivity",
                    title_x='Photon Energy (eV)', title_y='Beamline transmitivity')
    h5w.add_dataset(out[0, :], out[4, :], entry_name=f'work_dist_{working_distance*1e3}_mm',
                    dataset_name = "Beam Size Horizontal (DL)",
                    title_x = 'Photon Energy (eV)', title_y = 'FWHM_H (um)')
    h5w.add_dataset(out[0, :], out[5, :], entry_name=f'work_dist_{working_distance*1e3}_mm',
                    dataset_name = "Beam Size Vertical (DL)",
                    title_x = 'Photon Energy (eV)', title_y = 'FWHM_V (um)')
    h5w.add_dataset(out[0,:], out[6,:], entry_name=f'work_dist_{working_distance*1e3}_mm',
                    dataset_name="Intensity",
                    title_x='Photon Energy (eV)', title_y='Intensity')
    h5w.add_dataset(out[0,:], out[7,:], entry_name=f'work_dist_{working_distance*1e3}_mm',
                    dataset_name="Total initial rays",
                    title_x='Photon Energy (eV)', title_y='Total initial rays')    


    print(f'shadow_sxm2_work_dist_{working_distance*1e3}_mm.h5 has been save on disk')  

if __name__ == "__main__":
    #test
    #photon_energies = numpy.array([5000, 9000])
    #working_distances = [0.03]

    #simulations

    #working_distances = [0.03, 0.04, 0.05, 0.06] # 4 working distances   
    #photon_energies = numpy.arange(2500, 10000, 800) # 10 steps    
#
    #for working_distance in working_distances:
    #    write_output(working_distance, photon_energies, n_rays=10e6)
    pass