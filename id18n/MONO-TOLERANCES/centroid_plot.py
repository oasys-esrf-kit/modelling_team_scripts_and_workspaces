import numpy

if __name__ in ["__main__"]:
    from srxraylib.plot.gol import plot, plot_show

    #
    # main
    #


    rotation_axis = 'y'

    file_name5 = "centroid_rotation_%s_17keV_5.dat" % rotation_axis
    file_name6 = "centroid_rotation_%s_17keV_6.dat" % rotation_axis

    # file_name5 = "centroid_und_rotation_%s_17keV_5.dat" % rotation_axis
    # file_name6 = "centroid_und_rotation_%s_17keV_6.dat" % rotation_axis


    OFFSET5, CEN_x5, CEN_z5, CEN_xp5, CEN_zp5, SD_x5, SD_z5, SD_xp5, SD_zp5 = numpy.loadtxt(file_name5, unpack=True)
    OFFSET6, CEN_x6, CEN_z6, CEN_xp6, CEN_zp6, SD_x6, SD_z6, SD_xp6, SD_zp6 = numpy.loadtxt(file_name6, unpack=True)

    print(">>>>>>>>>>>>>>>>>>>5x: ", CEN_x5, CEN_xp5)
    print(">>>>>>>>>>>>>>>>>>>5z: ", CEN_z5, CEN_zp5)
    print(">>>>>>>>>>>>>>>>>>>6x: ", CEN_x6, CEN_xp6)
    print(">>>>>>>>>>>>>>>>>>>6z: ", CEN_z6, CEN_zp6)

    if rotation_axis == 'y': # use deg
        plot(numpy.degrees(OFFSET5), 1e6 * CEN_x5,
             numpy.degrees(OFFSET5), 1e6 * CEN_z5,
             numpy.degrees(OFFSET6), 1e6 * CEN_x6,
             numpy.degrees(OFFSET6), 1e6 * CEN_z6,
             # 1e6 * (OFFSET5), 1e6 * CEN_xp5,
             # 1e6 * (OFFSET5), 1e6 * CEN_zp5,
             # numpy.degrees(OFFSET5), -1e6 * 2 * 1 * OFFSET5 * numpy.cos(numpy.radians(6.67)),
             # numpy.degrees(OFFSET5), 1e6 * 2 * 1 * OFFSET6 * numpy.cos(numpy.radians(6.67)),
             title="rot axis '%s' (roll)" % (rotation_axis), yrange=[-2500, 2500],
             xtitle="rotation [deg]", ytitle="centroid [um]", legend=[
                "misalign channelcut 1 x",
                "misalign channelcut 1 z",
                "misalign channelcut 2 x",
                "misalign channelcut 2 z",
                "TEST 1",
                "TEST 2",
            ], grid=1, show=1)
    else: # use urad
        plot(1e6 * (OFFSET5), 1e6 * CEN_x5,
             1e6 * (OFFSET5), 1e6 * CEN_z5,
             1e6 * (OFFSET6), 1e6 * CEN_x6,
             1e6 * (OFFSET6), 1e6 * CEN_z6,
             # 1e6 * (OFFSET5), 1e6 * CEN_xp5,
             # 1e6 * (OFFSET5), 1e6 * CEN_zp5,
             # numpy.degrees(OFFSET5), -1e6 * 2 * 1 * OFFSET5 * numpy.cos(numpy.radians(6.67)),
             # numpy.degrees(OFFSET5), 1e6 * 2 * 1 * OFFSET6 * numpy.cos(numpy.radians(6.67)),
             title="rot axis '%s' (pitch)" % (rotation_axis), yrange=[-2500, 2500],
             xtitle="rotation [urad]", ytitle="centroid [um]", legend=[
                "misalign channelcut 1 x",
                "misalign channelcut 1 z",
                "misalign channelcut 2 x",
                "misalign channelcut 2 z",
                "TEST 1",
                "TEST 2",
            ], grid=1, show=1)

