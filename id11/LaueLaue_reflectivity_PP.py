import numpy
from oasys.util.oasys_util import get_fwhm
from srxraylib.plot.gol import plot, plot_show


def run_PP(ASYMMETRY_ANGLE=70.0, E0=100000.0):
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
        SCANFROM=-400.0,
        SCANTO=400.0,
        SCANPOINTS=20000,
        ENERGY=E0,
        ASYMMETRY_ANGLE=ASYMMETRY_ANGLE,
        THICKNESS=0.2,
        MOSAIC_FWHM=0.1,
        RSAG=15000.0,
        RMER=4200.0,
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


if __name__ == "__main__":

    OUT_ENERGY = []
    OUT_PEAK = []
    OUT_FWHM = []
    OUT_INTEG = []
    OUT_LEGEND = []


    OUT_ASYMMETRY_ANGLE = [70.0, 80.0, 90, 100, 110.]
    for ASYMMETRY_ANGLE in OUT_ASYMMETRY_ANGLE:
        ENERGY = []
        PEAK = []
        FWHM = []
        INTEG = []
        for i, E0_keV in enumerate([40, 60, 80, 100, 120, 140]):
            energy, ref_s = run_PP(ASYMMETRY_ANGLE=ASYMMETRY_ANGLE, E0=1e3 * E0_keV)
            peak = ref_s.max()
            fwhm, _, _ = get_fwhm(ref_s, energy)
            integ = numpy.trapz(ref_s, energy)

            print("peak, fwhm, int: ", peak, fwhm, integ)

            ENERGY.append (E0_keV)
            PEAK.append (peak)
            FWHM.append (fwhm)
            INTEG.append (integ)
        OUT_LEGEND.append("Asym angle: %d deg" % ASYMMETRY_ANGLE)

        # plot(ENERGY, FWHM,  title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='FWHM',  xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(FWHM).max()], show=0)
        # plot(ENERGY, PEAK,  title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='PEAK',  xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(PEAK).max()], show=0)
        # plot(ENERGY, INTEG, title="asymmetry: %d deg" % (ASYMMETRY_ANGLE), ytitle='INTEG', xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(INTEG).max()], show=0)
        # plot_show()

        OUT_ENERGY.append(ENERGY)
        OUT_FWHM.append(FWHM)
        OUT_PEAK.append(PEAK)
        OUT_INTEG.append(INTEG)


    plot(OUT_ENERGY[0], OUT_FWHM[0],
         OUT_ENERGY[1], OUT_FWHM[1],
         OUT_ENERGY[2], OUT_FWHM[2],
         OUT_ENERGY[3], OUT_FWHM[3],
         OUT_ENERGY[4], OUT_FWHM[4],
         title="", legend=OUT_LEGEND, ytitle='FWHM [urad]',
         xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(OUT_FWHM).max()], show=0)

    plot(OUT_ENERGY[0], OUT_PEAK[0],
         OUT_ENERGY[1], OUT_PEAK[1],
         OUT_ENERGY[2], OUT_PEAK[2],
         OUT_ENERGY[3], OUT_PEAK[3],
         OUT_ENERGY[4], OUT_PEAK[4],
         title="", legend=OUT_LEGEND, ytitle='PEAK',
         xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(OUT_PEAK).max()], show=0)

    plot(OUT_ENERGY[0], OUT_INTEG[0],
         OUT_ENERGY[1], OUT_INTEG[1],
         OUT_ENERGY[2], OUT_INTEG[2],
         OUT_ENERGY[3], OUT_INTEG[3],
         OUT_ENERGY[4], OUT_INTEG[4],
         title="", legend=OUT_LEGEND, ytitle='INTEGRATED REFLECTIVITY',
         xtitle="Energy [keV]", grid=1, yrange=[0, 1.2 * numpy.array(OUT_INTEG).max()], show=0)

    plot_show()
