# import shadow4
# print(dir(shadow4))
# from shadow4.beam.s4_beam import S4Beam
# from shadow4.beamline.s4_beamline import S4Beamline
# print(dir(shadow4))
import numpy
from shadow4.beamline.optical_elements.compound.s4_compound import S4Compound, S4CompoundElement


def get_optical_element_instance_channel_cut(
        crystal_separation=0,
        roll=0, pitch=0, yaw=0,  # applied in this order
        T=[0, 0, 0],
        use_mirrors=False,
):
    from shadow4.optical_surfaces.s4_conic import S4Conic
    from shadow4.beamline.optical_elements.crystals.s4_conic_crystal import S4ConicCrystal, S4ConicCrystalElement
    from shadow4.beamline.optical_elements.mirrors.s4_conic_mirror import S4ConicMirror, S4ConicMirrorElement

    try:
        name = self.getNode().title
    except:
        name = "Channel Cut Crystal Monochromator"

    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_conic_crystal import S4ConicCrystal

    ccc1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -0.5 * crystal_separation]
    ccc2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.5 * crystal_separation]

    # Rx, beta
    R_pitch = [[1, 0, 0],
               [0, numpy.cos(pitch), -numpy.sin(pitch)],
               [0, numpy.sin(pitch), numpy.cos(pitch)]]

    # Ry, gamma
    R_roll = [[numpy.cos(roll), 0, numpy.sin(roll)],
              [0, 1, 0],
              [-numpy.sin(roll), 0, numpy.cos(roll)]]

    # Rz, alpha
    R_yaw = [[numpy.cos(yaw), -numpy.sin(yaw), 0],
             [numpy.sin(yaw), numpy.cos(yaw), 0],
             [0, 0, 1]]

    # R = Rz Rx Ry
    R = numpy.array(R_yaw) @ numpy.array(R_pitch) @ numpy.array(R_roll)

    print(">>> R: ", R)

    conic_coefficients1 = S4Conic.rotate_and_translate_coefficients(ccc1, R, T)

    conic_coefficients2 = S4Conic.rotate_and_translate_coefficients(ccc2, R, T)

    if use_mirrors:

        optical_element1 = S4ConicMirror(name='Generic Mirror', boundary_shape=boundary_shape,
                                         conic_coefficients=conic_coefficients1,
                                         f_reflec=0, f_refl=6, file_refl='<none>',
                                         refraction_index=0.99999 + 0.001j,
                                         coating_material='Ni', coating_density=8.902, coating_roughness=0)
        optical_element2 = S4ConicMirror(name='Generic Mirror', boundary_shape=boundary_shape,
                                         conic_coefficients=conic_coefficients2,
                                         f_reflec=0, f_refl=6, file_refl='<none>',
                                         refraction_index=0.99999 + 0.001j,
                                         coating_material='Ni', coating_density=8.902, coating_roughness=0)
    else:
        optical_element1 = S4ConicCrystal(name='Generic Crystal',
                                          boundary_shape=boundary_shape,
                                          conic_coefficients=conic_coefficients1,
                                          material='Si', miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                          f_bragg_a=False, asymmetry_angle=0.0,
                                          is_thick=1, thickness=0.001,
                                          f_central=1, f_phot_cent=0, phot_cent=5000.0,
                                          file_refl='bragg.dat',
                                          f_ext=0,
                                          material_constants_library_flag=1,
                                          # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                          )

        optical_element2 = S4ConicCrystal(name='Generic Crystal',
                                          boundary_shape=boundary_shape,
                                          conic_coefficients=conic_coefficients2,
                                          material='Si', miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                          f_bragg_a=False, asymmetry_angle=0.0,
                                          is_thick=1, thickness=0.001,
                                          f_central=0, f_phot_cent=0, phot_cent=5000.0,
                                          file_refl='bragg.dat',
                                          f_ext=0,
                                          material_constants_library_flag=1,
                                          # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                          )
    #
    #
    #
    return S4Compound(name=name, oe_list=[optical_element1, optical_element2])

def run_beamline_17keV(
    crystal_separation=0,
    q2=0.0,
    theta_bragg_deg=10.0,
    use_mirrors=1,
    use_undulator=0,
    pitch=0,
    roll=0,
    yaw=0, ):

    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    if use_undulator:  # undulator
        # electron beam
        from shadow4.sources.s4_electron_beam import S4ElectronBeam

        electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
        electron_beam.set_sigmas_all(sigma_x=3.01836e-05, sigma_y=3.63641e-06, sigma_xp=4.36821e-06,
                                     sigma_yp=1.37498e-06)
        electron_beam.set_dispersion_all(0, 0, 0, 0)

        # magnetic structure
        from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian

        source = S4UndulatorGaussian(
            period_length=0.042,  # syned Undulator parameter (length in m)
            number_of_periods=38.571,  # syned Undulator parameter
            photon_energy=17000.0,  # Photon energy (in eV)
            delta_e=20.0,  # Photon energy width (in eV)
            ng_e=100,  # Photon energy scan number of points
            flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_energy_spread=0,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
            harmonic_number=1,  # harmonic number
            flag_autoset_flux_central_cone=0,  # value to set the flux peak
            flux_central_cone=10000000000.0,  # value to set the flux peak
        )

        # light source
        from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource

        light_source = S4UndulatorGaussianLightSource(name='Undulator Gaussian', electron_beam=electron_beam,
                                                      magnetic_structure=source, nrays=50000, seed=5676561)
    else:  # pencil

        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

        light_source = SourceGeometrical(name='Geometrical Source', nrays=500, seed=5676561)
        light_source.set_spatial_type_point()
        light_source.set_depth_distribution_off()
        light_source.set_angular_distribution_flat(hdiv1=0.000000, hdiv2=0.000000, vdiv1=0.000000, vdiv2=0.000000)
        light_source.set_energy_distribution_uniform(value_min=16990.0, value_max=17010.0, unit='eV')
        light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)

    beam = light_source.get_beam()
    beam0 = beam.duplicate()

    beamline.set_light_source(light_source)

    optical_element = get_optical_element_instance_channel_cut(
        crystal_separation=crystal_separation,
        roll=roll, pitch=pitch, yaw=yaw,  # applied in this order
        T=[0, 0, 0],
        use_mirrors=use_mirrors,
    )

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=0, q=q2,
                                     angle_radial=numpy.radians(90 - theta_bragg_deg),
                                     angle_azimuthal=0,
                                     angle_radial_out=numpy.radians(90 + theta_bragg_deg))

    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0,
                                           rotation_x=0, rotation_y=0, rotation_z=0)

    beamline_element = S4CompoundElement(
        optical_element=optical_element,
        coordinates=coordinates,
        movements=movements,
        input_beam=beam)

    beam, footprints = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # print(beamline.info())

    return beam0, beam, footprints

def plot_2d(beam, footprints, irange=[-0.0001, 0.0001], frange=[-0.009,0.009]):
    from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show

    # image
    ticket = beam.histo2(1, 3, nbins_h=100, nbins_v=100,
                         xrange=irange, yrange=irange, nolost=1, ref=23)

    title = "BEAM I: %.1f " % ticket['intensity']
    if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
    if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

    plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                               title=title, xtitle="column 1", ytitle="column 3",
                               cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=0)
    # footprint 1
    ticket = footprints[0].histo2(2, 1, nbins_h=100, nbins_v=100,
                                  xrange=frange, yrange=frange, nolost=1, ref=23)

    title = "FOOTPRINT FIRST CRYSTAL: %.1f " % ticket['intensity']
    if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
    if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

    plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                               title=title, xtitle="column 2", ytitle="column 1",
                               cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=0)

    # footprint 2
    ticket = footprints[1].histo2(2, 1, nbins_h=100, nbins_v=100,
                                  xrange=frange, yrange=frange, nolost=1, ref=23)

    title = "FOOTPRINT SECOND CRYSTAL: %.1f " % ticket['intensity']
    if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
    if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

    plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                               title=title, xtitle="column 2", ytitle="column 1",
                               cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=1)

def print_centroids(beam1, title='', factor=1.0):
    x = beam1.get_column(1, nolost=1)
    y = beam1.get_column(2, nolost=1)
    z = beam1.get_column(3, nolost=1)
    xp = beam1.get_column(4, nolost=1)
    yp = beam1.get_column(5, nolost=1)
    zp = beam1.get_column(6, nolost=1)
    w = beam1.get_column(23, nolost=1)

    arrays = [x, y, z, xp, yp, zp]
    titles = ['x', 'y', 'z', 'xp', 'yp', 'zp']

    # arrays = [x, y, z]
    # titles = ['x', 'y', 'z']

    print("\n-------------------------", title, "factor= %f" % factor)
    for j, array in enumerate(arrays):
        average = numpy.average(array, weights=w)
        variance = numpy.average((array - average) ** 2, weights=w)
        print(titles[j], factor * average, "+/-", factor * numpy.sqrt(variance))
    print("intensity: ", beam1.intensity(nolost=1))

if __name__ in ["__main__"]:
    from srxraylib.plot.gol import plot, plot_show

    #
    # main
    #
    do_calculate = 1
    rotation_axis = 'y'
    use_mirrors = 1
    use_undulator = 1
    theta_bragg_deg = 6.679637445440589
    theta_bragg = numpy.radians(theta_bragg_deg)


    # q1 = 0.1
    # crystal_separation = q1 * numpy.sin(numpy.radians(theta_bragg_deg))
    # q2 = 0.9

    q1 = 0.1
    crystal_separation = 0.005 # * numpy.sin(numpy.radians(theta_bragg_deg))
    q2 = 166.0
    q1 = crystal_separation / numpy.sin(theta_bragg)


    print(">>> crystal_separation, q1: ", crystal_separation, q1)


    theory_beam_coordinate = crystal_separation / numpy.sin(theta_bragg) * numpy.cos(numpy.pi/2 - 2 * theta_bragg)
    print(">>> theory beam separation: ", theory_beam_coordinate )
    print(">>> theory beam separation 2: ", crystal_separation / numpy.sin(theta_bragg) * numpy.sin(2 * theta_bragg),
          2 * crystal_separation * numpy.cos(theta_bragg))


    file_name = "check_s4rotations_rotation_channelcut_%s_17keV.dat" % (rotation_axis)
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
                beam0, beam1, footprints = run_beamline_17keV(use_undulator=use_undulator, use_mirrors=use_mirrors,
                                                              crystal_separation=crystal_separation, pitch=OFFSET[i],
                                                              theta_bragg_deg=theta_bragg_deg, q2=q2)
            elif rotation_axis == 'y':
                beam0, beam1, footprints = run_beamline_17keV(use_undulator=use_undulator, use_mirrors=use_mirrors,
                                                              crystal_separation=crystal_separation, roll=OFFSET[i],
                                                              theta_bragg_deg=theta_bragg_deg, q2=q2)
            elif rotation_axis == 'z':
                beam0, beam1, footprints = run_beamline_17keV(use_undulator=use_undulator, use_mirrors=use_mirrors,
                                                              crystal_separation=crystal_separation, yaw=OFFSET[i],
                                                              theta_bragg_deg=theta_bragg_deg, q2=q2)

            print_centroids(beam0, title='SOURCE', factor=1e6)
            print_centroids(footprints[0], title='CRYSTAL 1', factor=1e3)
            print_centroids(footprints[1], title='CRYSTAL 2', factor=1e3)
            print_centroids(beam1, title='IMAGE', factor=1e6)


            x = beam1.get_column(1, nolost=1)
            z = beam1.get_column(3, nolost=1)
            xp = beam1.get_column(4, nolost=1)
            zp = beam1.get_column(6, nolost=1)
            w = beam1.get_column(23, nolost=1)
            titles = ['x', 'z', 'xp', 'zp']

            # plot_2d(beam1, footprints, irange=[-1e-3, 1e-3])

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

        # numpy.savetxt(file_name, numpy.column_stack((OFFSET, CEN_x, CEN_z, CEN_xp, CEN_zp, SD_x, SD_z, SD_xp, SD_zp)))
        # print("File %s written to disk." % file_name)



    else:
        OFFSET, CEN_x, CEN_z, CEN_xp, CEN_zp, SD_x, SD_z, SD_xp, SD_zp = numpy.loadtxt(file_name, unpack=True)





    print(">>>> CEN_z: ", CEN_z)
    if rotation_axis == 'x':
        theory_old = -(2 * q1 * OFFSET + 2 * q2 * (OFFSET - OFFSET))

        theory = -(2 * crystal_separation * numpy.sin(theta_bragg) * OFFSET + 2 * q2 * (0.0 * OFFSET))

        print(">>>> crystal_separation, q1: ", crystal_separation, q1)

        plot(1e6 * OFFSET, 1e6 * (CEN_x - CEN_x[CEN_x.size // 2]),
             1e6 * OFFSET, 1e6 * (CEN_z - CEN_z[CEN_z.size // 2]), # theory_beam_coordinate),
             1e6 * OFFSET, 1e6 * theory_old,
             1e6 * OFFSET, 1e6 * theory,
             title="rot axis '%s', d=%.3f m, q2=%.3f m" % (rotation_axis, crystal_separation, q2),
             xtitle="rotation1 [urad]", ytitle="centroid [um]", legend=["x", "z", "pitch equation old", "pitch equation NEW"],
             linestyle=[None, '--', ':',':'],
             yrange=[-50,50],
             figsize=(10,4), grid=1, show=1)

    elif rotation_axis == 'y':

        theory = (2 * crystal_separation / numpy.sin(theta_bragg) * OFFSET + 2 * q2 * (0.0 * OFFSET)) * numpy.sin(theta_bragg)

        plot(numpy.degrees(OFFSET), 1e6 * (CEN_x - CEN_x[CEN_x.size // 2]),
             numpy.degrees(OFFSET), 1e6 * (CEN_z - CEN_z[CEN_z.size // 2]), # theory_beam_coordinate),
             numpy.degrees(OFFSET), 1e6 * theory,
             title="rot axis '%s', d=%.3f m, q2=%.3f m" % (rotation_axis, crystal_separation, q2),
             xtitle="rotation1 [deg]", ytitle="centroid [um]", legend=["x", "z", "roll equation"],
             linestyle=[None, None, ':'],
             yrange=[-50,50],
             figsize=(10,4), grid=1, show=1)

    elif rotation_axis == 'z':
        theory = OFFSET * 0
        plot(numpy.degrees(OFFSET), 1e6 * CEN_x,
             numpy.degrees(OFFSET), 1e6 * (CEN_z - theory_beam_coordinate),
             numpy.degrees(OFFSET), 1e6 * theory,
             title="rot axis '%s', q1=%.3f m, q2=%.3f m" % (rotation_axis, q1, q2),
             xtitle="rotation1 [deg]", ytitle="centroid [um]", legend=["x", "z", "yaw equation"],
             linestyle=[None, None, ':'],
             yrange=[-50,50],
             figsize=(10,4), grid=1, show=1)

