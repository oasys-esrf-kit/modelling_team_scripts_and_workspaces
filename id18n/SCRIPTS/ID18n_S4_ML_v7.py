import numpy
def run_beamline(angle_radial=numpy.pi/2-4.66667090125186e-3, pencil=False, external_M2_shape=False):
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    if not pencil:
        # electron beam
        from shadow4.sources.s4_electron_beam import S4ElectronBeam

        electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
        electron_beam.set_sigmas_all(sigma_x=2.87429e-05, sigma_y=5.15531e-06, sigma_xp=4.18025e-06, sigma_yp=1.93975e-06)

        # magnetic structure
        from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian

        source = S4UndulatorGaussian(
            period_length=0.0186,  # syned Undulator parameter (length in m)
            number_of_periods=134.40860215053763,  # syned Undulator parameter
            photon_energy=17000.0,  # Photon energy (in eV)
            delta_e=8.0,  # Photon energy width (in eV)
            ng_e=100,  # Photon energy scan number of points
            flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_energy_spread=1,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
            harmonic_number=3,  # harmonic number
            flag_autoset_flux_central_cone=1,  # value to set the flux peak
            flux_central_cone=1729338018179448.2,  # value to set the flux peak
        )

        # light source
        from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource

        light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam,
                                                      magnetic_structure=source, nrays=150000, seed=0)
        beam = light_source.get_beam()

    else:

        #
        #
        #
        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

        light_source = SourceGeometrical(name='SourceGeometrical', nrays=5000, seed=5676561)
        light_source.set_spatial_type_point()
        light_source.set_depth_distribution_off()
        light_source.set_angular_distribution_flat(hdiv1=0.000000, hdiv2=0.000000, vdiv1=0.000000, vdiv2=0.000000)
        light_source.set_energy_distribution_singleline(1000.000000, unit='eV')
        light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
        beam = light_source.get_beam()

    beamline.set_light_source(light_source)

    # optical element number XX

    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen

    optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator', boundary_shape=boundary_shape,
                               i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                               i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=27.3, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement

    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror

    optical_element = S4PlaneMirror(name='Plane Mirror', boundary_shape=boundary_shape,
                                    f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                    coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=2.2, q=0.3214294464, angle_radial=angle_radial, angle_azimuthal=4.71238898,
                                     angle_radial_out=1.566129656)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement

    beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements,
                                            input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beam1 = beam.duplicate()


    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    if external_M2_shape:
        # optical element number XX
        from syned.beamline.shape import Rectangle
        boundary_shape = Rectangle(x_left=-0.005, x_right=0.005, y_bottom=-0.365, y_top=0.365)

        from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror
        optical_element = S4PlaneMirror(name='Plane Mirror', boundary_shape=boundary_shape,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)
        ideal_mirror = optical_element
        boundary_shape = None

        from shadow4.beamline.optical_elements.mirrors.s4_numerical_mesh_mirror import S4NumericalMeshMirror
        optical_element = S4NumericalMeshMirror(name='Numerical Mesh Mirror', boundary_shape=boundary_shape,
                                                xx=None, yy=None, zz=None,
                                                surface_data_file='M2_shape.h5',
                                                f_reflec=0, f_refl=0, file_refl='', refraction_index=1,
                                                coating_material='', coating_density=1, coating_roughness=0)
        numerical_mesh_mirror = optical_element
        from syned.beamline.shape import Rectangle
        boundary_shape = Rectangle(x_left=-0.005, x_right=0.005, y_bottom=-0.365, y_top=0.365)

        from shadow4.beamline.optical_elements.mirrors.s4_additional_numerical_mesh_mirror import \
            S4AdditionalNumericalMeshMirror
        optical_element = S4AdditionalNumericalMeshMirror(name='ideal + error Mirror', ideal_mirror=ideal_mirror,
                                                          numerical_mesh_mirror=numerical_mesh_mirror)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=0.3214294464, q=24.25714111, angle_radial=1.566129656,
                                         angle_azimuthal=3.141592654, angle_radial_out=angle_radial)
        movements = None
        from shadow4.beamline.optical_elements.mirrors.s4_additional_numerical_mesh_mirror import \
            S4AdditionalNumericalMeshMirrorElement
        beamline_element = S4AdditionalNumericalMeshMirrorElement(optical_element=optical_element,
                                                                  coordinates=coordinates, movements=movements,
                                                                  input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)
    else:
        boundary_shape = None

        from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirror

        optical_element = S4SphereMirror(name='Sphere Mirror', boundary_shape=boundary_shape,
                                         surface_calculation=1, is_cylinder=1, cylinder_direction=0,
                                         convexity=1, radius=5786.911762, p_focus=30.142859, q_focus=24.257141,
                                         grazing_angle=0.004667,
                                         f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                         coating_material='Si', coating_density=2.33, coating_roughness=0)

        from syned.beamline.element_coordinates import ElementCoordinates

        coordinates = ElementCoordinates(p=0.3214294464, q=24.25714111, angle_radial=1.566129656, angle_azimuthal=3.141592654,
                                         angle_radial_out=angle_radial)
        movements = None
        from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirrorElement

        beamline_element = S4SphereMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements,
                                                 input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

    footprint1 = footprint.duplicate()

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty

    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=1.570796327, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement

    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen

    optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator HSS', boundary_shape=boundary_shape,
                               i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                               i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement

    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    return beam, beam1, footprint1, beamline

#
# main
#
if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show

    # theta_grazing = 0.0035
    theta_grazing = 4.66667090125186e-3
    # theta_grazing = 0.007
    angle_radial=numpy.pi/2 - theta_grazing

    beam, beam1, footprint1, beamline = run_beamline(angle_radial=angle_radial, pencil=False, external_M2_shape=1)
    print(beamline.info())

    ##########
    ticket = footprint1.histo2(2, 1, nbins_h=100, nbins_v=100, xrange=[-0.4,0.4], yrange=[-0.0025,0.0025], nolost=1, ref=23)

    title = "FOOTPRINT - I: %.1f " % ticket['intensity']
    if ticket['fwhm_h'] is not None: title += "FWHM H: %.4f " % ticket['fwhm_h']
    if ticket['fwhm_v'] is not None: title += "FWHM V: %.4f " % ticket['fwhm_v']

    plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
        title=title, xtitle="M2 (length) [m]", ytitle="M2 (width) [m]",
        cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=0)


    ##############
    # ticket = beam1.histo2(3, 1, nbins_h=100, nbins_v=100, xrange=[-0.0015,0.0015], yrange=[-0.0015,0.0015], nolost=1, ref=23)
    #
    # title = "INTERMEDIATE SCREEN I: %.1f " % ticket['intensity']
    # if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
    # if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']
    #
    # plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
    #     title=title, xtitle="column 3", ytitle="column 1",
    #     cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=0)


    ##############
    centroid_x = numpy.average(beam.get_column(1, nolost=1), weights=beam.get_column(23, nolost=1))
    ticket = beam.histo2(1, 3, nbins_h=100, nbins_v=100,
                         xrange = [-1000e-6 + centroid_x, 1000e-6 + centroid_x], yrange = [-1000e-6, 1000e-6],
                         nolost=1, ref=23)

    title = "IMAGE I: %.1f " % ticket['intensity']
    if ticket['fwhm_h'] is not None: title += "FWHM H: %.1f " % (1e6 * ticket['fwhm_h'])
    if ticket['fwhm_v'] is not None: title += "FWHM V: %.1f " % (1e6 * ticket['fwhm_v'])

    plot_image_with_histograms(ticket['histogram'], 1e6 * ticket['bin_h_center'], 1e6 * ticket['bin_v_center'],
        title=title, xtitle="H [um]", ytitle="V [um]",
        cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=1)







