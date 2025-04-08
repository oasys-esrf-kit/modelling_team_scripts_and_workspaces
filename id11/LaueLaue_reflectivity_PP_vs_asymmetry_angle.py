import numpy
from oasys.util.oasys_util import get_fwhm
from srxraylib.plot.gol import plot, plot_show


def run_PP_energy_scan(ASYMMETRY_ANGLE=70.0, E0=100000.0, RMER=1000.0, THICKNESS=0.5):
    #
    # script to calculate crystal diffraction profiles (created by XOPPY:crystal)
    #

    import numpy
    from xoppylib.crystals.tools import bragg_calc2, run_diff_pat
    import xraylib
    from dabax.dabax_xraylib import DabaxXraylib

    #
    # run bragg_calc (preprocessor) and create file xcrystal.bra
    #
    bragg_dictionary = bragg_calc2(
        descriptor="Si",
        hh=1,
        kk=1,
        ll=1,
        temper=1.0,
        emin=E0-1000,
        emax=E0+1000,
        estep=4.004,
        ANISO_SEL=0,
        fileout="xcrystal.bra",
        do_not_prototype=0,  # 0=use site groups (recommended), 1=use all individual sites
        verbose=False,
        material_constants_library=xraylib,
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




def run_PP(ASYMMETRY_ANGLE=70.0, E0=100000.0, RMER=1000.0, THICKNESS=0.5):
    #
    # script to calculate crystal diffraction profiles (created by XOPPY:crystal)
    #

    import numpy
    from xoppylib.crystals.tools import bragg_calc2, run_diff_pat
    import xraylib
    from dabax.dabax_xraylib import DabaxXraylib

    #
    # run bragg_calc (preprocessor) and create file xcrystal.bra
    #
    bragg_dictionary = bragg_calc2(
        descriptor="Si",
        hh=1,
        kk=1,
        ll=1,
        temper=1.0,
        emin=E0-100,
        emax=E0+100,
        estep=0.4,
        ANISO_SEL=0,
        fileout="xcrystal.bra",
        do_not_prototype=0,  # 0=use site groups (recommended), 1=use all individual sites
        verbose=False,
        material_constants_library=xraylib,
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
        SCAN=2,
        UNIT=1,
        SCANFROM=-1000.0,
        SCANTO=1000.0,
        SCANPOINTS=20000,
        ENERGY=E0,
        ASYMMETRY_ANGLE=ASYMMETRY_ANGLE,
        THICKNESS=THICKNESS, # !!!!!!!
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

def get_Ru_Rl(chi_deg=0.0, Energy_keV=50, F1=100.0, F2=100):
    lattice_constant_A = 5.431  # Silicon lattice constant in Angstroms
    d_111_A = lattice_constant_A / numpy.sqrt(3)  # Si(111) interplanar spacing in Angstroms
    theta_B_rad = numpy.arcsin(1e10 * wavelength_m(Energy_keV) / (2 * d_111_A))  # in radians
    chi = numpy.radians(chi_deg)
    Ru = 2 * (numpy.cos(chi + theta_B_rad) / F1 + numpy.cos(chi - theta_B_rad) / F2) ** (-1)
    Rl = 2 * (numpy.cos(chi - theta_B_rad) / F1 + numpy.cos(chi + theta_B_rad) / F2) ** (-1)
    return Ru, Rl

def get_Rw(chi_deg=0.0, Energy_keV=50, F1=100.0):
    lattice_constant_A = 5.431  # Silicon lattice constant in Angstroms
    d_111_A = lattice_constant_A / numpy.sqrt(3)  # Si(111) interplanar spacing in Angstroms
    theta_B_rad = numpy.arcsin(1e10 * wavelength_m(Energy_keV) / (2 * d_111_A))  # in radians
    chi = numpy.radians(chi_deg)
    Rw = F1  / (2 * numpy.cos( chi + theta_B_rad) )
    Rw2 = F1 / (2 * numpy.cos(-chi + theta_B_rad) )
    return Rw, Rw2, theta_B_rad

import scipy.constants as codata


def wavelength_m(energy_keV):
    return codata.c * codata.h / codata.e / (1e3 * energy_keV)  # Wavelength in m

if __name__ == "__main__":

    OUT_ENERGY = []
    OUT_PEAK = []
    OUT_FWHM = []
    OUT_INTEG = []
    OUT_LEGEND = []
    OUT_BRAGG_ANGLE_RAD = []
    OUT_BANDWIDTH = []

    npoints = 21

    # OUT_ASYMMETRY_ANGLE = numpy.ndarray.tolist(numpy.linspace(90-60, 90+60, npoints))
    OUT_ASYMMETRY_ANGLE = numpy.linspace(90 - 80, 90 + 80, npoints)

    # ID31: F1=105.5, 5 mm thick.
    # ID11: F1 = 32.157
    # ID15: F1 = 51.0
    # F1 = 32.157 # 105.2 # 51.0
    # F2 = F1 # F1 / 3.0


    beamline = 'ID11'

    if beamline == 'ID15':
        # ID 15
        F1 = 51.0
        THICKNESS = 0.3 # in cm
        F2 = 65.0 - F1
        i_scan = 1 # 0: theta, 1=energy
    elif beamline == 'ID31':
        # ID 31
        F1 = 105.2
        THICKNESS = 0.5 # in cm
        F2 = 120.0 - F1
        i_scan = 1 # 0: theta, 1=energy
    elif beamline == 'ID11':
        # ID 31
        F1 = 32.157
        THICKNESS = 0.3 # in cm
        F2 = 92.0 - F1
        i_scan = 1 # 0: theta, 1=energy
    else:
        raise NotImplementedError()


    #
    # angle scan
    #
    if i_scan == 0:
        for ie, E0_keV in enumerate([90.0, 70.0, 50.0]):
            ENERGY = []
            PEAK = []
            FWHM = []
            INTEG = []
            R = []
            BANDWIDTH = []
            for ia, ASYMMETRY_ANGLE in enumerate(OUT_ASYMMETRY_ANGLE):

                RMER = 3000.0
                if 1:
                    # # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", ASYMMETRY_ANGLE-90, E0_keV, F1, F2)
                    # Ru, Rl = get_Ru_Rl(chi_deg=ASYMMETRY_ANGLE-90, Energy_keV=E0_keV, F1=F1, F2=F2)
                    Rw1, Rw2, theta_B_rad = get_Rw(chi_deg=ASYMMETRY_ANGLE - 90, Energy_keV=E0_keV, F1=F1)
                    # # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", Ru, Rl)
                    RMER = 2 * Rw1 * 100.0
                angle_urad, ref_s = run_PP(ASYMMETRY_ANGLE=ASYMMETRY_ANGLE, E0=1e3 * E0_keV, RMER=RMER)
                # plot(energy, ref_s)
                ref_s = ref_s ** 2 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                peak = ref_s.max()
                fwhm, _, _ = get_fwhm(ref_s, angle_urad)
                integ = numpy.trapz(ref_s, angle_urad * 1e-6)
                bandwidth = 100 * (fwhm * 1e-6) / numpy.tan(numpy.array(theta_B_rad))

                print("peak, fwhm, int, bandwidth: ", peak, fwhm, integ, bandwidth)

                ENERGY.append (E0_keV)
                PEAK.append (peak)
                FWHM.append (fwhm)
                INTEG.append (integ)
                R.append(RMER/100) # in m
                BANDWIDTH.append(bandwidth)

            OUT_LEGEND.append("E = %d keV" % E0_keV)

            # plot(ENERGY, FWHM,  title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='FWHM',  xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(FWHM).max()], show=0)
            # plot(ENERGY, PEAK,  title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='PEAK',  xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(PEAK).max()], show=0)
            # plot(ENERGY, INTEG, title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='INTEG', xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(INTEG).max()], show=0)
            # plot_show()

            OUT_ENERGY.append(ENERGY)
            OUT_FWHM.append(FWHM)
            OUT_PEAK.append(PEAK)
            OUT_INTEG.append(INTEG)
            OUT_BRAGG_ANGLE_RAD.append(theta_B_rad)
            OUT_BANDWIDTH.append(BANDWIDTH)

    #
    # energy scan
    #
    elif i_scan == 1:
        for ie, E0_keV in enumerate([90.0, 70.0, 50.0]):
            ENERGY = []
            PEAK = []
            FWHM = []
            INTEG = []
            R = []
            BANDWIDTH = []
            for ia, ASYMMETRY_ANGLE in enumerate(OUT_ASYMMETRY_ANGLE):

                RMER = 3000.0
                if 1:
                    # # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", ASYMMETRY_ANGLE-90, E0_keV, F1, F2)
                    # Ru, Rl = get_Ru_Rl(chi_deg=ASYMMETRY_ANGLE-90, Energy_keV=E0_keV, F1=F1, F2=F2)
                    Rw1, Rw2, theta_B_rad = get_Rw(chi_deg=ASYMMETRY_ANGLE - 90, Energy_keV=E0_keV, F1=F1)
                    # # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", Ru, Rl)
                    RMER = 2 * Rw1 * 100.0
                energy, ref_s = run_PP_energy_scan(ASYMMETRY_ANGLE=ASYMMETRY_ANGLE, E0=1e3 * E0_keV, RMER=RMER)
                # plot(energy, ref_s)
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
                R.append(RMER / 100)  # in m
                BANDWIDTH.append(bandwidth)

            OUT_LEGEND.append("E = %d keV" % E0_keV)

            # plot(ENERGY, FWHM,  title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='FWHM',  xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(FWHM).max()], show=0)
            # plot(ENERGY, PEAK,  title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='PEAK',  xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(PEAK).max()], show=0)
            # plot(ENERGY, INTEG, title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='INTEG', xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(INTEG).max()], show=0)
            # plot_show()

            OUT_ENERGY.append(ENERGY)
            OUT_FWHM.append(FWHM)
            OUT_PEAK.append(PEAK)
            OUT_INTEG.append(INTEG)
            OUT_BRAGG_ANGLE_RAD.append(theta_B_rad)
            OUT_BANDWIDTH.append(BANDWIDTH)


    if i_scan == 0:
        ytitle = 'INTEGRATED REFLECTIVITY [urad]'
    elif i_scan == 1:
        ytitle = 'INTEGRATED REFLECTIVITY [eV]'
    plot(OUT_ASYMMETRY_ANGLE - 90, numpy.array(OUT_INTEG[0]),
         OUT_ASYMMETRY_ANGLE - 90, numpy.array(OUT_INTEG[1]),
         OUT_ASYMMETRY_ANGLE - 90, numpy.array(OUT_INTEG[2]),
         title=beamline, legend=OUT_LEGEND, ytitle=ytitle,
         xtitle="Asymmetry angle [deg]", grid=1, yrange=[0, 1.2 * (numpy.array(OUT_INTEG)).max()], show=0)

    plot(OUT_ASYMMETRY_ANGLE - 90, R,
         title=beamline, ytitle=' R [m]',
         xtitle="Asymmetry angle [deg]", grid=1, show=0)

    plot(OUT_ASYMMETRY_ANGLE - 90, numpy.array(OUT_BANDWIDTH[0]),
         OUT_ASYMMETRY_ANGLE - 90, numpy.array(OUT_BANDWIDTH[1]),
         OUT_ASYMMETRY_ANGLE - 90, numpy.array(OUT_BANDWIDTH[2]),
         title=beamline, legend=OUT_LEGEND, ytitle='BANDWIDTH [%]',
         xtitle="Asymmetry angle [deg]", grid=1, yrange=[0, 1.2 * (numpy.array(OUT_BANDWIDTH)).max()], show=0)



    for ie, E0_keV in enumerate([90.0, 70.0, 50.0]):
        for ia, ASYMMETRY_ANGLE in enumerate(OUT_ASYMMETRY_ANGLE):
            # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", ASYMMETRY_ANGLE-90, E0_keV, F1, F2)
            Ru, Rl = get_Ru_Rl(chi_deg=ASYMMETRY_ANGLE-90, Energy_keV=E0_keV, F1=F1, F2=F2)
            Rw, Rw2, theta_B_rad = get_Rw(chi_deg=ASYMMETRY_ANGLE - 90, Energy_keV=E0_keV, F1=F1)
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", E0_keV, ASYMMETRY_ANGLE-90, Ru, Rl, Rw, Rw2)

    plot_show()