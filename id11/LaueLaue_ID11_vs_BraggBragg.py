import numpy
# from oasys.util.oasys_util import get_fwhm
from srxraylib.plot.gol import plot, plot_show


def get_fwhm(histogram, bins, ret0=None):
    fwhm = ret0
    quote = ret0
    coordinates = None

    if histogram.size > 1:
        quote = numpy.max(histogram)*0.5
        cursor = numpy.where(histogram >= quote)

        if histogram[cursor].size > 1:
            bin_size    = bins[1]-bins[0]
            fwhm        = bin_size*(cursor[0][-1]-cursor[0][0])
            coordinates = (bins[cursor[0][0]], bins[cursor[0][-1]])

    return fwhm, quote, coordinates

def run_BraggBragg_energy_scan(ASYMMETRY_ANGLE=00.0, E0=100000.0, RMER=1000.0, THICKNESS=0.5, H=1, K=1, L=1):
    #
    # script to calculate crystal diffraction profiles (created by XOPPY:crystal)
    #

    import numpy
    from xoppylib.crystals.tools import bragg_calc2, run_diff_pat
    import xraylib
    from dabax.dabax_xraylib import DabaxXraylib
    dx = DabaxXraylib(dabax_repository="./DABAX",
                 file_f0="f0_InterTables.dat",
                 file_f1f2="f1f2_WindtWithCompton.dat",
                 file_CrossSec = "CrossSec_EPDL97.dat",
                 file_Crystals="Crystals.dat",)

    #
    # run bragg_calc (preprocessor) and create file xcrystal.bra
    #
    bragg_dictionary = bragg_calc2(
        descriptor="Si",
        hh=H,
        kk=K,
        ll=L,
        temper=1.0,
        emin=E0-50,
        emax=E0+50,
        estep=0.204,
        ANISO_SEL=0,
        fileout="xcrystal.bra",
        do_not_prototype=0,  # 0=use site groups (recommended), 1=use all individual sites
        verbose=False,
        material_constants_library=dx,
    )

    #
    # run external (fortran) diff_pat (note that some parameters may not be used)
    #
    run_diff_pat(
        bragg_dictionary,
        preprocessor_file="xcrystal.bra",
        descriptor="Si",
        MOSAIC=0,
        GEOMETRY=0,
        SCAN=3,
        UNIT=1,
        SCANFROM=E0-50,
        SCANTO  =E0+50,
        SCANPOINTS=2000,
        ENERGY=-E0,
        ASYMMETRY_ANGLE=ASYMMETRY_ANGLE,
        THICKNESS=THICKNESS,
        MOSAIC_FWHM=0.1,
        RSAG=15000.0,
        RMER=RMER,
        ANISOTROPY=0,
        POISSON=0.22,
        CUT="2 -1 -1 ; 1 1 1 ; 0 0 0",
        FILECOMPLIANCE="mycompliance.dat",
    )

    #
    # example plot
    #
    if 0:
        from srxraylib.plot.gol import plot
        data = numpy.loadtxt("diff_pat.dat", skiprows=5)
        plot(data[:, 0], data[:, -1], data[:, 0], data[:, -2], ytitle='Crystal reflectivity',
             legend=['s-polarization', 'p-polarization'])

    #
    # end script
    #
    data = numpy.loadtxt("diff_pat.dat", skiprows=5)
    return data[:, 0], data[:, -1]




def run_PP_energy_scan(ASYMMETRY_ANGLE=70.0, E0=100000.0, RMER=1000.0, THICKNESS=0.5, H=1, K=1, L=1):
    #
    # script to calculate crystal diffraction profiles (created by XOPPY:crystal)
    #

    import numpy
    from xoppylib.crystals.tools import bragg_calc2, run_diff_pat
    import xraylib
    from dabax.dabax_xraylib import DabaxXraylib
    dx = DabaxXraylib(dabax_repository="./DABAX",
                 file_f0="f0_InterTables.dat",
                 file_f1f2="f1f2_WindtWithCompton.dat",
                 file_CrossSec = "CrossSec_EPDL97.dat",
                 file_Crystals="Crystals.dat",)

    #
    # run bragg_calc (preprocessor) and create file xcrystal.bra
    #
    bragg_dictionary = bragg_calc2(
        descriptor="Si",
        hh=H,
        kk=K,
        ll=L,
        temper=1.0,
        emin=E0-1000,
        emax=E0+1000,
        estep=4.004,
        ANISO_SEL=0,
        fileout="xcrystal.bra",
        do_not_prototype=0,  # 0=use site groups (recommended), 1=use all individual sites
        verbose=False,
        material_constants_library=dx,
    )

    #
    # run external (fortran) diff_pat (note that some parameters may not be used)
    #
    run_diff_pat(
        bragg_dictionary,
        preprocessor_file="xcrystal.bra",
        descriptor="Si",
        MOSAIC=3,
        GEOMETRY=1,
        SCAN=3,
        UNIT=1,
        SCANFROM=E0-1000,
        SCANTO  =E0+1000,
        SCANPOINTS=20000,
        ENERGY=-E0,
        ASYMMETRY_ANGLE=ASYMMETRY_ANGLE,
        THICKNESS=THICKNESS,
        MOSAIC_FWHM=0.1,
        RSAG=15000.0,
        RMER=RMER,
        ANISOTROPY=0,
        POISSON=0.22,
        CUT="2 -1 -1 ; 1 1 1 ; 0 0 0",
        FILECOMPLIANCE="mycompliance.dat",
    )

    #
    # example plot
    #
    if 0:
        from srxraylib.plot.gol import plot
        data = numpy.loadtxt("diff_pat.dat", skiprows=5)
        plot(data[:, 0], data[:, -1], data[:, 0], data[:, -2], ytitle='Crystal reflectivity',
             legend=['s-polarization', 'p-polarization'])

    #
    # end script
    #
    data = numpy.loadtxt("diff_pat.dat", skiprows=5)
    return data[:, 0], data[:, -1]

def get_Rw(chi_deg=0.0, Energy_keV=50, F1=100.0, H=1, K=1, L=1):
    lattice_constant_A = 5.431  # Silicon lattice constant in Angstroms
    d_111_A = lattice_constant_A / numpy.sqrt(H*H+K*K+L*L)  # Si(111) interplanar spacing in Angstroms
    theta_B_rad = numpy.arcsin(1e10 * wavelength_m(Energy_keV) / (2 * d_111_A))  # in radians
    chi = numpy.radians(chi_deg)
    # Rw = F1  / (2 * numpy.cos( chi + theta_B_rad) )
    # Rw2 = F1 / (2 * numpy.cos(-chi + theta_B_rad) )

    Rw = F1  / (2 * numpy.cos( chi + theta_B_rad) )
    Rw2 = F1 / (2 * numpy.cos( chi - theta_B_rad) )
    return Rw, Rw2, theta_B_rad

import scipy.constants as codata


def wavelength_m(energy_keV):
    return codata.c * codata.h / codata.e / (1e3 * energy_keV)  # Wavelength in m

if __name__ == "__main__":

    H = 1
    K = 1
    L = 1
    beamline = 'ID11'
    npoints = 21

    if beamline == 'ID15':
        # ID 15
        F1 = 51.0
        THICKNESS = 0.3 # in cm
        CHI_ANGLE = 39.0
    elif beamline == 'ID31':
        # ID 31
        F1 = 105.2
        THICKNESS = 0.5 # in cm
        CHI_ANGLE = - 36.0
    elif beamline == 'ID11':
        # ID 31
        F1 = 31.6
        THICKNESS = 0.25 # in cm
        CHI_ANGLE = -12.5
    else:
        raise NotImplementedError()


    #
    #
    #

    ASYMMETRY_ANGLE = 90 - CHI_ANGLE

    R_m = numpy.linspace(30.0, 100, npoints)

    energies_kev = numpy.linspace(20, 80, 30)

    #
    # energy scan
    #

    ENERGY = []
    PEAK = []
    FWHM = []
    INTEG = []
    BANDWIDTH = []

    BB_PEAK = []
    BB_FWHM = []
    BB_INTEG = []
    BB_BANDWIDTH = []



    for ie, E0_keV in enumerate(energies_kev):

        Rw1, Rw2, theta_B_rad = get_Rw(chi_deg=CHI_ANGLE, Energy_keV=E0_keV, F1=F1, H=H, K=1, L=1)
        RMER = 2 * Rw1 * 100.0


        energy, ref_s = run_PP_energy_scan(ASYMMETRY_ANGLE=ASYMMETRY_ANGLE, E0=1e3 * E0_keV, RMER=RMER,
                                           THICKNESS=THICKNESS, H=H, K=K, L=L)
        ref_s = ref_s ** 2  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        peak = ref_s.max()
        fwhm, _, _ = get_fwhm(ref_s, energy)
        integ = numpy.trapz(ref_s, energy)
        bandwidth = 100 * fwhm / (E0_keV * 1e3)

        print("peak, fwhm, int, bandwidth: ", peak, fwhm, integ, bandwidth)

        ENERGY.append(E0_keV)
        PEAK.append(peak)
        FWHM.append(fwhm)
        INTEG.append(integ)
        BANDWIDTH.append(bandwidth)

        #
        # Bragg Bragg
        #
        energy, ref_s = run_BraggBragg_energy_scan(ASYMMETRY_ANGLE=0, E0=1e3 * E0_keV, RMER=RMER,
                                           THICKNESS=1, H=H, K=K, L=L)

        ref_s = ref_s ** 2  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        peak = ref_s.max()
        fwhm, _, _ = get_fwhm(ref_s, energy)
        integ = numpy.trapz(ref_s, energy)
        bandwidth = 100 * fwhm / (E0_keV * 1e3)

        print("peak, fwhm, int, bandwidth: ", peak, fwhm, integ, bandwidth)

        BB_PEAK.append(peak)
        BB_FWHM.append(fwhm)
        BB_INTEG.append(integ)
        BB_BANDWIDTH.append(bandwidth)

    title_laue = " chi=%.1f deg, t=%.1f mm" % (CHI_ANGLE, THICKNESS * 10)
    title = "%s Si%d%d%d" % (beamline, H, K, L)

    plot(energies_kev, numpy.array(INTEG),
         energies_kev, numpy.array(BB_INTEG),
         title=title, ytitle='INTEGRATED REFLECTIVITY [eV]',
         legend=["LaueLaue %s" % title_laue, "BraggBragg"], color=['r','g'],
         xtitle="Photon Energy [keV]", grid=1, show=0)

    # plot(energies_kev, numpy.array(BANDWIDTH),
    #      energies_kev, numpy.array(BB_BANDWIDTH),
    #      title=title, ytitle='BANDWIDTH [%]',
    #      legend=["LaueLaue %s" % title_laue, "BraggBragg"], color=['r','g'],
    #      xtitle="Photon Energy [keV]", grid=1, show=0)

    plot(energies_kev, numpy.array(BANDWIDTH),
         energies_kev, numpy.array(BB_BANDWIDTH),
         title=title, ytitle='BANDWIDTH [%]',
         legend=["LaueLaue %s" % title_laue, "BraggBragg"], color=['r','g'],
         xtitle="Photon Energy [keV]", grid=1, ylog=1, show=0)
    plot_show()