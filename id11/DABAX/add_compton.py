def get_f2():
    #
    # script to make the calculations (created by XOPPY:xf1f2)
    #
    from xoppylib.scattering_functions.xoppy_calc_f1f2 import xoppy_calc_f1f2
    import xraylib
    from dabax.dabax_xraylib import DabaxXraylib

    out_dict = xoppy_calc_f1f2(
        descriptor="Si",
        density=2.33,
        MAT_FLAG=0,
        CALCULATE=1,
        GRID=0,
        GRIDSTART=5000.0,
        GRIDEND=25000.0,
        GRIDN=100,
        THETAGRID=0,
        ROUGH=0.0,
        THETA1=2.0,
        THETA2=5.0,
        THETAN=50,
        DUMP_TO_FILE=0,
        FILE_NAME="f1f2.dat",
        material_constants_library=DabaxXraylib(file_f1f2="f1f2_Windt.dat"),
    )

    #
    # example plot
    #
    if 0:
        from srxraylib.plot.gol import plot, plot_image
        try:
            plot(out_dict["data"][0, :], out_dict["data"][-1, :],
                 xtitle=out_dict["labels"][0], ytitle=out_dict["labels"][1], title="xf1f2",
                 xlog=True, ylog=True, show=True)
        except:
            plot_image(out_dict["data2D"], out_dict["dataX"], out_dict["dataY"],
                       xtitle='Energy [eV]', ytitle='Theta [mrad]', title='Reflectivity',
                       aspect='auto', show=True)
    #
    # end script
    #
    return out_dict["data"][0, :], out_dict["data"][-1, :]

if __name__ == "__main__":
    import xraylib
    import numpy
    from srxraylib.plot.gol import plot

    a = numpy.loadtxt("Si_windt.txt")


    energies_eV, y = a[:, 0], a[:, 2] # get_f2()

    energies_keV = energies_eV * 1e-3
    # print(energies_keV, y)

    Z_Si = xraylib.SymbolToAtomicNumber("Si")

    # Compute cross section ratios
    ratios = []
    ratio = 1.0
    for E in energies_keV:
        try:
            compton = xraylib.CS_Compt(Z_Si, E)
            photoelectric = xraylib.CS_Photo(Z_Si, E)
            ratio = (photoelectric +  compton) / photoelectric if photoelectric != 0 else float('inf')
            ratios.append(ratio)
        except:
            ratios.append(ratio)

        # print(">>", E, ratio)

    ratios = numpy.array(ratios)

    plot(energies_keV, y,
         energies_keV, ratios,
         energies_keV, y * ratios,
         # xrange=[10,100],
         xlog=1, ylog=1, show=True, legend=['original f2', 'weight=(compton+photoel)/photoel', 'weighted f2'],
         xtitle="Photon energy [eV]")

    for i in range(energies_keV.size):
        print(energies_eV[i], + a[i, 1] , a[i, 2] * ratios[i])
