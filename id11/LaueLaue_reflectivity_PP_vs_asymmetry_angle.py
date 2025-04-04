import numpy
from oasys.util.oasys_util import get_fwhm
from srxraylib.plot.gol import plot, plot_show


def run_PP(ASYMMETRY_ANGLE=70.0, E0=100000.0, RMER=1000.0):
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
        SCANFROM=-600.0,
        SCANTO=600.0,
        SCANPOINTS=30000,
        ENERGY=E0,
        ASYMMETRY_ANGLE=ASYMMETRY_ANGLE,
        THICKNESS=0.3, # !!!!!!!
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

import scipy.constants as codata


def wavelength_m(energy_keV):
    return codata.c * codata.h / codata.e / (1e3 * energy_keV)  # Wavelength in m

if __name__ == "__main__":

    OUT_ENERGY = []
    OUT_PEAK = []
    OUT_FWHM = []
    OUT_INTEG = []
    OUT_LEGEND = []

    npoints = 21

    # OUT_ASYMMETRY_ANGLE = numpy.ndarray.tolist(numpy.linspace(90-60, 90+60, npoints))
    OUT_ASYMMETRY_ANGLE = numpy.linspace(90 - 80, 90 + 80, npoints)
    F1 = 32.0
    F2 = 92.0 - F1








    if 1:

        for ie, E0_keV in enumerate([90.0, 70.0, 50.0]):
            ENERGY = []
            PEAK = []
            FWHM = []
            INTEG = []
            for ia, ASYMMETRY_ANGLE in enumerate(OUT_ASYMMETRY_ANGLE):

                RMER = 500.0
                # # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", ASYMMETRY_ANGLE-90, E0_keV, F1, F2)
                # Ru, Rl = get_Ru_Rl(chi_deg=ASYMMETRY_ANGLE-90, Energy_keV=E0_keV, F1=F1, F2=F2)
                # # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", Ru, Rl)
                # RMER = Ru * 100.0
                energy, ref_s = run_PP(ASYMMETRY_ANGLE=ASYMMETRY_ANGLE, E0=1e3 * E0_keV, RMER=RMER)
                peak = ref_s.max()
                fwhm, _, _ = get_fwhm(ref_s, energy)
                integ = numpy.trapz(ref_s, energy)

                print("peak, fwhm, int: ", peak, fwhm, integ)

                ENERGY.append (E0_keV)
                PEAK.append (peak)
                FWHM.append (fwhm)
                INTEG.append (integ)
            OUT_LEGEND.append("E = %d keV" % E0_keV)

            # plot(ENERGY, FWHM,  title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='FWHM',  xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(FWHM).max()], show=0)
            # plot(ENERGY, PEAK,  title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='PEAK',  xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(PEAK).max()], show=0)
            # plot(ENERGY, INTEG, title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='INTEG', xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(INTEG).max()], show=0)
            # plot_show()

            OUT_ENERGY.append(ENERGY)
            OUT_FWHM.append(FWHM)
            OUT_PEAK.append(PEAK)
            OUT_INTEG.append(INTEG)


        plot(OUT_ASYMMETRY_ANGLE - 90, OUT_FWHM[0],
             OUT_ASYMMETRY_ANGLE - 90, OUT_FWHM[1],
             OUT_ASYMMETRY_ANGLE - 90, OUT_FWHM[2],
             # OUT_ASYMMETRY_ANGLE, OUT_FWHM[3],
             # OUT_ASYMMETRY_ANGLE, OUT_FWHM[4],
             title="", legend=OUT_LEGEND, ytitle='FWHM [urad]',
             xtitle="Asymmetry angle [deg]", grid=1, yrange=[0, 1.2 * numpy.array(OUT_FWHM).max()], show=0)

        plot(OUT_ASYMMETRY_ANGLE - 90, OUT_PEAK[0],
             OUT_ASYMMETRY_ANGLE - 90, OUT_PEAK[1],
             OUT_ASYMMETRY_ANGLE - 90, OUT_PEAK[2],
             # OUT_ASYMMETRY_ANGLE, OUT_PEAK[3],
             # OUT_ASYMMETRY_ANGLE, OUT_PEAK[4],
             title="", legend=OUT_LEGEND, ytitle='PEAK',
             xtitle="Asymmetry angle [deg]", grid=1, yrange=[0, 1.2 * numpy.array(OUT_PEAK).max()], show=0)

        plot(OUT_ASYMMETRY_ANGLE - 90, OUT_INTEG[0],
             OUT_ASYMMETRY_ANGLE - 90, OUT_INTEG[1],
             OUT_ASYMMETRY_ANGLE - 90, OUT_INTEG[2],
             # OUT_ASYMMETRY_ANGLE, OUT_INTEG[3],
             # OUT_ASYMMETRY_ANGLE, OUT_INTEG[4],
             title="", legend=OUT_LEGEND, ytitle='INTEGRATED REFLECTIVITY',
             xtitle="Asymmetry angle [deg]", grid=1, yrange=[0, 1.2 * numpy.array(OUT_INTEG).max()], show=0)



    for ie, E0_keV in enumerate([90.0, 70.0, 50.0]):
        for ia, ASYMMETRY_ANGLE in enumerate(OUT_ASYMMETRY_ANGLE):
            RMER = 10000.0
            # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", ASYMMETRY_ANGLE-90, E0_keV, F1, F2)
            Ru, Rl = get_Ru_Rl(chi_deg=ASYMMETRY_ANGLE-90, Energy_keV=E0_keV, F1=F1, F2=F2)
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", E0_keV, ASYMMETRY_ANGLE-90, Ru, Rl)

    plot_show()