# import shadow4
# print(dir(shadow4))
# from shadow4.beam.s4_beam import S4Beam
# from shadow4.beamline.s4_beamline import S4Beamline
# print(dir(shadow4))
import numpy

def run_beamline_17keV(
    q1=0.0,
    q2=1.0,
    rotation_x_1 = 0.0,
    rotation_x_2 = 0.0,
    rotation_y_1 = 0.0,
    rotation_y_2 = 0.0,
    rotation_z_1 = 0.0,
    rotation_z_2 = 0.0, ):


    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    #
    #
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    light_source = SourceGeometrical(name='Geometrical Source', nrays=50000, seed=5676561)
    light_source.set_spatial_type_point()
    light_source.set_depth_distribution_off()
    light_source.set_angular_distribution_flat(hdiv1=0.000000, hdiv2=0.000000, vdiv1=0.000000, vdiv2=0.000000)
    light_source.set_energy_distribution_uniform(value_min=16990.000000, value_max=17010.000000, unit='eV')
    light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
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
    coordinates = ElementCoordinates(p=0, q=q1, angle_radial=1.45421465, angle_azimuthal=0,
                                     angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0,
                                           rotation_x=rotation_x_1,
                                           rotation_y=rotation_y_1,
                                           rotation_z=rotation_z_1)
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
    coordinates = ElementCoordinates(p=0, q=q2, angle_radial=1.454214616, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.454214616)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0,
                                           rotation_x=rotation_x_2,
                                           rotation_y=rotation_y_2,
                                           rotation_z=rotation_z_2)
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

    return beam

if __name__ in ["__main__"]:
    from srxraylib.plot.gol import plot, plot_show

    #
    # main
    #
    do_calculate = 0
    rotation_axis = 'z'
    q1 = 0.1
    q2 = 1.0 - q1
    rot2_over_rot1 = -1
    file_name = "check_s4rotations_rotation_%s_17keV.dat" % (rotation_axis)
    if do_calculate:
        # WARNING: NO incremental result allowed!!"

        if rotation_axis == 'x':
            OFFSET = numpy.linspace(numpy.radians(-0.001), numpy.radians(0.001), 5)
        elif rotation_axis == 'y':
            OFFSET = numpy.linspace(numpy.radians(-0.1), numpy.radians(0.1), 5)
        elif rotation_axis == 'z':
            OFFSET = numpy.linspace(numpy.radians(-1), numpy.radians(1), 5)
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
            if rotation_axis == 'x':
                beam1 = run_beamline_17keV(q1=q1, q2=q2, rotation_x_1=rotation, rotation_x_2=rot2_over_rot1 * rotation)
            elif rotation_axis == 'y':
                beam1 = run_beamline_17keV(q1=q1, q2=q2, rotation_y_1=rotation, rotation_y_2=rot2_over_rot1 * rotation)
            elif rotation_axis == 'z':
                beam1 = run_beamline_17keV(q1=q1, q2=q2, rotation_z_1=rotation, rotation_z_2=rot2_over_rot1 * rotation)

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


    if rotation_axis == 'x':
        theory = -(2 * q1 * OFFSET + 2 * q2 * (OFFSET + rot2_over_rot1 * OFFSET))

        plot(numpy.degrees(OFFSET), 1e6 * CEN_x,
             numpy.degrees(OFFSET), 1e6 * CEN_z,
             numpy.degrees(OFFSET), 1e6 * theory,
             title="rot axis '%s', q1=%.3f m, q2=%.3f m, rotation2 = %.3f * rotation1" % (rotation_axis, q1, q2, rot2_over_rot1),
             xtitle="rotation1 [deg]", ytitle="centroid [um]", legend=["x", "z", "pitch equation"],
             linestyle=[None, None, ':'],
             yrange=[-50,50],
             figsize=(10,4), grid=1, show=1)
    elif rotation_axis == 'y':
        theory = (2 * q1 * OFFSET + 2 * q2 * (OFFSET - rot2_over_rot1 * OFFSET)) * numpy.sin(numpy.radians(6.7))

        plot(numpy.degrees(OFFSET), 1e6 * CEN_x,
             numpy.degrees(OFFSET), 1e6 * CEN_z,
             numpy.degrees(OFFSET), 1e6 * theory,
             title="rot axis '%s', q1=%.3f m, q2=%.3f m, rotation2 = %.3f * rotation1" % (rotation_axis, q1, q2, rot2_over_rot1),
             xtitle="rotation1 [deg]", ytitle="centroid [um]", legend=["x", "z", "roll equation"],
             linestyle=[None, None, ':'],
             # yrange=[-50,50],
             figsize=(10,4), grid=1, show=1)

    elif rotation_axis == 'z':
        theory = OFFSET * 0
        plot(numpy.degrees(OFFSET), 1e6 * CEN_x,
             numpy.degrees(OFFSET), 1e6 * CEN_z,
             numpy.degrees(OFFSET), 1e6 * theory,
             title="rot axis '%s', q1=%.3f m, q2=%.3f m, rotation2 = %.3f * rotation1" % (rotation_axis, q1, q2, rot2_over_rot1),
             xtitle="rotation1 [deg]", ytitle="centroid [um]", legend=["x", "z", "yaw equation"],
             linestyle=[None, None, ':'],
             # yrange=[-50,50],
             figsize=(10,4), grid=1, show=1)

