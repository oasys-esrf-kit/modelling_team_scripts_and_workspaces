import numpy
import scipy.constants as codata

def sinc(x):
    return numpy.sin(x) / x

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

def calculate_1d_scan(energy=17000.0, D=400e-6, f=0.050, npoints=2000, window=200e-9):
    wavelength = codata.h * codata.c / codata.e / energy
    x = numpy.linspace(-window / 2, window / 2, npoints)
    k = 2 * numpy.pi / wavelength
    a = k * D * x / 2 / f
    y = sinc(a)**2
    return x, a, y, k



if __name__ == "__main__":
    from srxraylib.plot.gol import plot
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 20})


    D = 400e-6
    f = 0.050
    npoints = 2000
    window = 200e-9
    energies = numpy.array([7000.0, 17000, 35000, 40000])

    SCANS = []
    for energy in energies:
        x, a, y, k = calculate_1d_scan(energy=energy, D=D, f=f, npoints=npoints, window=window)

        SCANS.append(y)
        fwhm_sinc2, _, _ = get_fwhm(y, a)
        fwhm, _, _ = get_fwhm(y, 1e9 * x)

        fwhm_theor1 = 1e9 * fwhm_sinc2 / (k * D / 2 / f)
        fwhm_theor2 = 1e9 * 2.78 / (k * D / 2 / f)

        print("energy: %d keV, fwhm: numerical: %.3f nm, theor1: %.3f nm, theor2: %.3f nm"  %(1e-3 * energy, fwhm, fwhm_theor1, fwhm_theor2))



    plot(1e9 * x, SCANS[0],
         1e9 * x, SCANS[1],
         1e9 * x, SCANS[2],
         1e9 * x, SCANS[3],
         xtitle="x [nm]", ytitle="Intensity", legend=["7 keV","17 keV", "35 keV", "40 keV"],
         title="D: %.1f $\mu$m, f=%.1f mm" % (D*1e6, f*1e3), figsize=(14,8), show=0)
    plt.grid()


    energies = numpy.linspace(7000, 40000, 100)
    FWHM = []
    for energy in energies:
        wavelength = codata.h * codata.c / codata.e / energy
        k = 2 * numpy.pi / wavelength
        fwhm_theor2 = 1e9 * 2.78 / (k * D / 2 / f)
        FWHM.append(fwhm_theor2)

    plot(1e-3 * energies, FWHM, xtitle="Photon energy [keV]", ytitle="FWHM [nm]",
         title="Diffraction limited size; D: %.1f $\mu$m, f=%.1f mm" % (D*1e6, f*1e3), figsize=(14,8), show=0)
    plt.grid()
    plt.show()