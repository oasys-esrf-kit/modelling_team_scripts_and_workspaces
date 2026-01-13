# import shadow4
# print(dir(shadow4))
# from shadow4.beam.s4_beam import S4Beam
# from shadow4.beamline.s4_beamline import S4Beamline
# print(dir(shadow4))
import numpy

def run_beamline_17keV(crystal_number=-1, rotation_axis=None, rotation=1.74533e-05):
    rotation_x_1 = 0.0
    rotation_x_2 = 0.0
    rotation_x_3 = 0.0
    rotation_x_4 = 0.0

    rotation_y_1 = 0.0
    rotation_y_2 = 0.0
    rotation_y_3 = 0.0
    rotation_y_4 = 0.0

    if rotation_axis == 'x':
        if crystal_number == 1:
            rotation_x_1 = rotation
        elif crystal_number == 2:
            rotation_x_2 = rotation
        elif crystal_number == 3:
            rotation_x_3 = rotation
        elif crystal_number == 4:
            rotation_x_4 = rotation
        elif crystal_number == 5:
            rotation_x_1 = rotation
            rotation_x_2 = -rotation
        elif crystal_number == 6:
            rotation_x_3 = rotation
            rotation_x_4 = -rotation
    elif rotation_axis == 'y':
        if crystal_number == 1:
            rotation_y_1 = rotation
        elif crystal_number == 2:
            rotation_y_2 = rotation
        elif crystal_number == 3:
            rotation_y_3 = rotation
        elif crystal_number == 4:
            rotation_y_4 = rotation
        elif crystal_number == 5:
            rotation_y_1 = rotation
            rotation_y_2 = rotation
        elif crystal_number == 6:
            rotation_y_3 = rotation
            rotation_y_4 = rotation

    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    #
    #
    #
    from shadow4.sources.s4_light_source_from_file import S4LightSourceFromFile
    light_source = S4LightSourceFromFile(name='Shadow4 File Reader', file_name='C:/Users/srio/Oasys/tmp.h5',
                                         simulation_name='run001', beam_name='begin')
    beam = light_source.get_beam()

    beamline.set_light_source(light_source)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Generic Crystal',
                                     boundary_shape=boundary_shape, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=17000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.45421465, angle_azimuthal=0, angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_z=0,
                                           rotation_y=rotation_y_1,
                                           rotation_x=rotation_x_1)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Generic Crystal',
                                     boundary_shape=boundary_shape, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=17000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.45421465, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_z=0,
                                           rotation_y=rotation_y_2,
                                           rotation_x=rotation_x_2)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Generic Crystal',
                                     boundary_shape=boundary_shape, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=17000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.45421465, angle_azimuthal=0, angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_z=0,
                                           rotation_y=rotation_y_3,
                                           rotation_x=rotation_x_3)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Generic Crystal',
                                     boundary_shape=boundary_shape, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=17000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.454214616, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.454214616)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_z=0,
                                           rotation_y=rotation_y_4,
                                           rotation_x=rotation_x_4)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # test plot
    if 0:
        from srxraylib.plot.gol import plot_scatter
        plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
                     title='(Intensity,Photon Energy)', plot_histograms=0)
        plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')

    print(">>>>>>", beam.intensity(nolost=1))

    return beam

if __name__ in ["__main__"]:
    from srxraylib.plot.gol import plot, plot_show

    #
    # main
    #
    do_calculate = 1
    rotation_axis = 'y'
    crystal_number = 6 # 5 = move first channel cut, 6=move second one
    xrange = [17000-10, 17000+10]

    file_name = "rotation_%s_17keV_%d.dat" % (rotation_axis, crystal_number)
    if do_calculate:
        # WARNING: NO incremental result allowed!!"
        # OFFSET = numpy.linspace(-100e-6, 100e-6, 51)
        if rotation_axis == 'x':
            OFFSET = numpy.linspace(-50e-6, 50e-6, 51)
        elif rotation_axis == 'y':
            OFFSET = numpy.linspace(numpy.radians(-3), numpy.radians(3), 51)
        INTENSITY = numpy.zeros_like(OFFSET)
        FWHM = numpy.zeros_like(OFFSET)

        for i, rotation in enumerate(OFFSET):
            print ("iteration %d of %d" % (i, OFFSET.size))
            beam = run_beamline_17keV(crystal_number=crystal_number, rotation_axis=rotation_axis, rotation=rotation)

            ticket = beam.histo1(26, nbins=100, xrange=xrange, nolost=1, ref=23)
            INTENSITY[i] = ticket['intensity']
            if ticket['fwhm'] is not None: FWHM[i] = ticket['fwhm']

            if 0:
                title = "I: %.1f " % ticket['intensity']
                if ticket['fwhm'] is not None: title += "FWHM: %f " % ticket['fwhm']
                plot(ticket['bin_path'], ticket['histogram_path'],
                     title=title, xtitle="column 26", show=1)

        numpy.savetxt(file_name, numpy.column_stack((OFFSET, INTENSITY, FWHM)))
        print("File %s written to disk." % file_name)



    else:
        OFFSET, INTENSITY, FWHM = numpy.loadtxt(file_name, unpack=True)


    # plot(1e6 * OFFSET, 1e6 * FWHM,
    #      xtitle="rotation_x [urad]", ytitle="fwhm [urad]", show=0)
    #
    # plot(1e6 * OFFSET, INTENSITY,
    #      xtitle="rotation_x [urad]", ytitle="intensity [a.u.]", show=0)


    from scipy.optimize import curve_fit

    def gaussian(x, A, mu, sigma):
        return A * numpy.exp(-(x - mu)**2 / (2 * sigma**2))

    p0 = [INTENSITY.max(), OFFSET[numpy.argmax(INTENSITY)], numpy.std(OFFSET)]
    popt, pcov = curve_fit(gaussian, OFFSET, INTENSITY, p0=p0)

    A, mu, sigma = popt

    height = 0.5
    width_at_height = 2 * sigma * numpy.sqrt(-2 * numpy.log(height))
    print("Width at %f of height: %f urad" % (height, 1e6 * width_at_height))

    height = 0.9
    width_at_height = 2 * sigma * numpy.sqrt(-2 * numpy.log(height))
    print("Width at %f of height: %f urad = %f mdeg" % (height, 1e6 * width_at_height, 1e3 * numpy.degrees(width_at_height)))


    # plot(1e6 * OFFSET, 1e6 * FWHM,
    #      xtitle="rotation_x [urad]", ytitle="FWHM [urad]", grid=1, show=0)


    plot(1e6 * OFFSET, INTENSITY,
         1e6 * OFFSET, gaussian(OFFSET, A, mu, sigma),
         xtitle="rotation_x [urad]", ytitle="intensity [a.u.]", legend=["ray tracing", "Gaussian fit FWHM=%.1f urad" % (2.355e6 * sigma)], grid=1, show=0)

    plot(numpy.degrees(OFFSET), INTENSITY,
         numpy.degrees(OFFSET), gaussian(OFFSET, A, mu, sigma),
         xtitle="rotation_x [deg", ytitle="intensity [a.u.]", legend=["ray tracing", "Gaussian fit FWHM=%.1f urad" % (2.355e6 * sigma)], grid=1, show=1)

