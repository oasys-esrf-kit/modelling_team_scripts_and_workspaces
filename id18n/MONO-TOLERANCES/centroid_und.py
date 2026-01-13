# import shadow4
# print(dir(shadow4))
# from shadow4.beam.s4_beam import S4Beam
# from shadow4.beamline.s4_beamline import S4Beamline
# print(dir(shadow4))
import numpy

def run_beamline_17keV(q=0.0, crystal_number=-1, rotation_axis=None, rotation=1.74533e-05):
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
    coordinates = ElementCoordinates(p=0, q=0.01, angle_radial=1.45421465, angle_azimuthal=0,
                                     angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0,
                                           rotation_x=rotation_x_1,
                                           rotation_y=rotation_y_1,
                                           rotation_z=0)
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
    coordinates = ElementCoordinates(p=0, q=0.01, angle_radial=1.45421465, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0,
                                           rotation_x=rotation_x_2,
                                           rotation_y=rotation_y_2,
                                           rotation_z=0)
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
    coordinates = ElementCoordinates(p=0, q=0.01, angle_radial=1.45421465, angle_azimuthal=0,
                                     angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0,
                                           rotation_x=rotation_x_3,
                                           rotation_y=rotation_y_3,
                                           rotation_z=0)
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
    coordinates = ElementCoordinates(p=0, q=0.01, angle_radial=1.45421465, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0,
                                           rotation_x=rotation_x_4,
                                           rotation_y=rotation_y_4,
                                           rotation_z=0)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=165.697, q=0, angle_radial=0, angle_azimuthal=4.71238898,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # test plot
    if 0:
        from srxraylib.plot.gol import plot_scatter
        plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
                     title='(Intensity,Photon Energy)', plot_histograms=0)
        plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')

    return beam

if __name__ in ["__main__"]:
    from srxraylib.plot.gol import plot, plot_show

    #
    # main
    #
    do_calculate = 1

    rotation_axis = 'x'
    crystal_number = 6
    q=0.01

    file_name = "centroid_und_rotation_%s_17keV_%d.dat" % (rotation_axis, crystal_number)
    if do_calculate:
        # WARNING: NO incremental result allowed!!"
        # OFFSET = numpy.linspace(-100e-6, 100e-6, 51)
        if rotation_axis == 'x':
            OFFSET = numpy.linspace(-50e-6, 50e-6, 11)
        elif rotation_axis == 'y':
            OFFSET = numpy.linspace(numpy.radians(-1.5), numpy.radians(1.5), 11)
        CEN_x = numpy.zeros_like(OFFSET)
        CEN_z = numpy.zeros_like(OFFSET)
        CEN_xp = numpy.zeros_like(OFFSET)
        CEN_zp = numpy.zeros_like(OFFSET)
        SD_x = numpy.zeros_like(OFFSET)
        SD_z = numpy.zeros_like(OFFSET)
        SD_xp = numpy.zeros_like(OFFSET)
        SD_zp = numpy.zeros_like(OFFSET)


        for i, rotation in enumerate(OFFSET):
            print ("iteration %d of %d" % (i, OFFSET.size))
            beam1 = run_beamline_17keV(q=q,
                                       crystal_number=crystal_number,
                                       rotation_axis=rotation_axis,
                                       rotation=rotation)

            # beam1.retrace(1.0)

            x = beam1.get_column(1, nolost=1)
            z = beam1.get_column(3, nolost=1)
            xp = beam1.get_column(4, nolost=1)
            zp = beam1.get_column(6, nolost=1)
            w = beam1.get_column(23, nolost=1)
            titles = ['x', 'z', 'xp', 'zp']

            for j, array in enumerate([x, z, xp, zp]):
                average = numpy.average(array, weights=w)
                variance = numpy.average((array - average) ** 2, weights=w)
                print(OFFSET[i], titles[j], average, "+/-", numpy.sqrt(variance))
                if j==0:
                    CEN_x[i] = average
                    SD_x[i] = numpy.sqrt(variance)
                elif j==1:
                    CEN_z[i] = average
                    SD_z[i] = numpy.sqrt(variance)
                elif j==2:
                    CEN_xp[i] = average
                    SD_xp[i] = numpy.sqrt(variance)
                elif j==3:
                    CEN_zp[i] = average
                    SD_zp[i] = numpy.sqrt(variance)



        numpy.savetxt(file_name, numpy.column_stack((OFFSET, CEN_x, CEN_z, CEN_xp, CEN_zp, SD_x, SD_z, SD_xp, SD_zp)))
        print("File %s written to disk." % file_name)



    else:
        OFFSET, CEN_x, CEN_z, CEN_xp, CEN_zp, SD_x, SD_z, SD_xp, SD_zp = numpy.loadtxt(file_name, unpack=True)


    plot(1e6 * OFFSET, 1e6 * CEN_x,
         1e6 * OFFSET, 1e6 * CEN_z,
         title="rot axis '%s'" % (rotation_axis),
         xtitle="rotation [urad]", ytitle="centroid [um]", legend=["x", "z"], grid=1, show=0)

    plot(numpy.degrees(OFFSET), 1e6 * CEN_x,
         numpy.degrees(OFFSET), 1e6 * CEN_z,
         title="rot axis '%s'" % (rotation_axis),
         xtitle="rotation [deg]", ytitle="centroid [um]", legend=["x", "z"], grid=1, show=1)

