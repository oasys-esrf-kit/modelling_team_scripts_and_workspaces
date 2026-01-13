import numpy


def get_fwhm(histogram, bins, ret0=None, height=0.5):
    fwhm = ret0
    quote = ret0
    coordinates = None

    if histogram.size > 1:
        quote = numpy.max(histogram)*height
        cursor = numpy.where(histogram >= quote)

        if histogram[cursor].size > 1:
            bin_size    = bins[1]-bins[0]
            fwhm        = bin_size*(cursor[0][-1]-cursor[0][0])
            coordinates = (bins[cursor[0][0]], bins[cursor[0][-1]])

    return fwhm, quote, coordinates

if __name__ in ["__main__"]:
    from srxraylib.plot.gol import plot, plot_show

    from scipy.optimize import curve_fit

    def gaussian(x, A, mu, sigma):
        return A * numpy.exp(-(x - mu)**2 / (2 * sigma**2))


    do_gaussian_fit = 0
    rotation_axis = 'y'
    rotation_axis = 'x'

    #
    #
    #

    if rotation_axis == 'x':
        OFFSET5, INTENSITY5, FWHM5 = numpy.loadtxt("rotation_x_17keV_5.dat", unpack=True)
        OFFSET6, INTENSITY6, FWHM6 = numpy.loadtxt("rotation_x_17keV_6.dat", unpack=True)
    elif rotation_axis == 'y':
        OFFSET5, INTENSITY5, FWHM5 = numpy.loadtxt("rotation_y_17keV_5.dat", unpack=True)
        OFFSET6, INTENSITY6, FWHM6 = numpy.loadtxt("rotation_y_17keV_6.dat", unpack=True)


    width_at_p9_5, _, _ = get_fwhm(INTENSITY5, OFFSET5, ret0=None, height=0.9)
    width_at_p9_6, _, _ = get_fwhm(INTENSITY6, OFFSET6, ret0=None, height=0.9)

    if do_gaussian_fit:
        p5 = [INTENSITY5.max(), OFFSET5[numpy.argmax(INTENSITY5)], numpy.std(OFFSET5)]
        p6 = [INTENSITY6.max(), OFFSET6[numpy.argmax(INTENSITY6)], numpy.std(OFFSET6)]



        popt5, pcov5 = curve_fit(gaussian, OFFSET5, INTENSITY5, p0=p5)
        popt6, pcov6 = curve_fit(gaussian, OFFSET6, INTENSITY6, p0=p6)


        A5, mu5, sigma5 = popt5
        A6, mu6, sigma6 = popt6

        for height in [0.5, 0.9]:

            width_at_height5 = 2 * sigma5 * numpy.sqrt(-2 * numpy.log(height))
            width_at_height6 = 2 * sigma6 * numpy.sqrt(-2 * numpy.log(height))
            print("Width at %f of height5 = %f urad = %f mdeg" % (height, 1e6 * width_at_height5, 1e3 * numpy.degrees(width_at_height5)))
            print("Width at %f of height6 = %f urad = %f mdeg" % (height, 1e6 * width_at_height6, 1e3 * numpy.degrees(width_at_height6)))

        plot(1e6 * OFFSET5 * 1  , INTENSITY5,
             1e6 * OFFSET5 * 1  , gaussian(OFFSET5, A5, mu5, sigma5),
             1e6 * OFFSET6 * 1 , INTENSITY6,
             1e6 * OFFSET6 * 1 , gaussian(OFFSET6, A6, mu6, sigma6),
             xtitle="rotation_x [urad]", ytitle="intensity [a.u.]",
             legend=["shadow4 rot channelcut 1 (width at 0.9 height = %f deg)" % numpy.degrees(width_at_p9_5), "Gaussian fit xtal 1 width at 0.9 max=%.1f urad" % (1e6 * width_at_height5),
                     "shadow4 rot channelcut 2 (width at 0.9 height = %f deg)" % numpy.degrees(width_at_p9_6), "Gaussian fit xtal 2 width at 0.9 max=%.1f urad" % (1e6 * width_at_height6),
                     ],
                    color=['b','b','r','r','g','g','k','k'],
                    linestyle=[None,":",None,":",None,":",None,":",],
                    grid=1, show=0)

        plot(numpy.degrees(OFFSET5), INTENSITY5,
             numpy.degrees(OFFSET5), gaussian(OFFSET5, A5, mu5, sigma5),
             numpy.degrees(OFFSET6), INTENSITY6,
             numpy.degrees(OFFSET6), gaussian(OFFSET6, A6, mu6, sigma6),
             xtitle="rotation_x [deg]", ytitle="intensity [a.u.]",
             legend=["shadow4 rot channelcut 1 (width at 0.9 height = %f deg)" % numpy.degrees(width_at_p9_5), "Gaussian fit xtal 1 width at 0.9 max=%.1f urad" % (1e6 * width_at_height5),
                     "shadow4 rot channelcut 2 (width at 0.9 height = %f deg)" % numpy.degrees(width_at_p9_6), "Gaussian fit xtal 2 width at 0.9 max=%.1f urad" % (1e6 * width_at_height6),
                     ],
                    color=['b','b','r','r','g','g','k','k'],
                    linestyle=[None,":",None,":",None,":",None,":",],
                    grid=1, show=1)
    else:

        plot(1e6 * OFFSET5, INTENSITY5,
             1e6 * OFFSET6, INTENSITY6,
             xtitle="rotation_%s [urad]" % rotation_axis, ytitle="intensity [a.u.]",
             legend=["shadow4 rot channelcut 1 (width at 0.9 height = %f urad)" % (1e6 * width_at_p9_5),
                     "shadow4 rot channelcut 2 (width at 0.9 height = %f urad)" % (1e6 * width_at_p9_6),
                     ],
                    color=['b','r',],
                    linestyle=[None,None,], figsize=(10, 5),
                    grid=1, show=0)

        plot(numpy.degrees(OFFSET5), INTENSITY5,
             numpy.degrees(OFFSET6), INTENSITY6,
             xtitle="rotation_%s [deg]" % rotation_axis, ytitle="intensity [a.u.]",
             legend=["shadow4 rot channelcut 1 (width at 0.9 height = %f deg)" % numpy.degrees(width_at_p9_5),
                     "shadow4 rot channelcut 2 (width at 0.9 height = %f deg)" % numpy.degrees(width_at_p9_6),
                     ],
                    color=['b','r',],
                    linestyle=[None,None,], figsize=(10, 5),
                    grid=1, show=1)