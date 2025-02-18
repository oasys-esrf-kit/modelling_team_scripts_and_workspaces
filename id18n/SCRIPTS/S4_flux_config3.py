

def run_config_a():
    global SEED
    print(">>>>  saved seed: ", SEED)
    SEED += 555
    print(">>>>  using seed: ", SEED)

    # DO NOT FORGET TO CHANGE seed=SEED!!!!!!!!!!!!!!!!
    # DO NOT FORGET TO return beam
    # DO NOT FORGET TO OPEN THE 50um SLIT
    ##################
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

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
        flag_autoset_flux_central_cone=0,  # value to set the flux peak
        flux_central_cone=1350000000000000.0,  # value to set the flux peak
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
    light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam,
                                                  magnetic_structure=source, nrays=150000, seed=SEED)
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
    coordinates = ElementCoordinates(p=2.2, q=0, angle_radial=1.53937691, angle_azimuthal=4.71238898,
                                     angle_radial_out=1.53937691)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
    beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                            movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.0065, x_right=0.0065, y_bottom=-0.06, y_top=0.06)

    from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirror
    optical_element = S4SphereMirror(name='Sphere Mirror', boundary_shape=boundary_shape,
                                     surface_calculation=0, is_cylinder=1, cylinder_direction=0,
                                     convexity=1, radius=1.000000, p_focus=28.300000, q_focus=11.700000,
                                     grazing_angle=0.031419,
                                     f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                     coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=2, q=22.5, angle_radial=1.53937691, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.53937691)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirrorElement
    beamline_element = S4SphereMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=1.570796327,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-1.4595e-05, x_right=1.4595e-05, y_bottom=-0.5, y_top=0.5)

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

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.03, y_top=0.03)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.014990,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=145.9, q=0.0335, angle_radial=1.555806392, angle_azimuthal=0,
                                     angle_radial_out=1.555806392)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.012537, y_top=0.012537)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=145.950000, q_focus=0.050000, grazing_angle=0.014990,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.0165, q=0, angle_radial=1.555806392, angle_azimuthal=1.570796327,
                                     angle_radial_out=1.555806392)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.011, x_right=0.011, y_bottom=-0.000375, y_top=0.000375)

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='CORRECTION !!!!  Generic Beam Screen/Slit/Stopper/Attenuator HSS',
                               boundary_shape=boundary_shape,
                               i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                               i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0.05, angle_radial=0, angle_azimuthal=4.71238898,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)


    return beam

# MONO AFTER
def run_config_b():
    global SEED
    print(">>>>  saved seed: ", SEED)
    SEED += 555
    print(">>>>  using seed: ", SEED)
    ##################

    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

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
        flag_autoset_flux_central_cone=0,  # value to set the flux peak
        flux_central_cone=1350000000000000.0,  # value to set the flux peak
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
    light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam,
                                                  magnetic_structure=source, nrays=150000, seed=SEED)
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
    coordinates = ElementCoordinates(p=2.2, q=0, angle_radial=1.53937691, angle_azimuthal=4.71238898,
                                     angle_radial_out=1.53937691)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
    beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                            movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.0065, x_right=0.0065, y_bottom=-0.06, y_top=0.06)

    from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirror
    optical_element = S4SphereMirror(name='Sphere Mirror', boundary_shape=boundary_shape,
                                     surface_calculation=0, is_cylinder=1, cylinder_direction=0,
                                     convexity=1, radius=1.000000, p_focus=28.300000, q_focus=11.700000,
                                     grazing_angle=0.031419,
                                     f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                     coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=2, q=22.5, angle_radial=1.53937691, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.53937691)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirrorElement
    beamline_element = S4SphereMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
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
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
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
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
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
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
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
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=3.141592654,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=1.570796327,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-1.4595e-05, x_right=1.4595e-05, y_bottom=-0.5, y_top=0.5)

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

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.03, y_top=0.03)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.014990,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=145.9, q=0.0335, angle_radial=1.555806392, angle_azimuthal=0,
                                     angle_radial_out=1.555806392)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.012537, y_top=0.012537)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=145.950000, q_focus=0.050000, grazing_angle=0.014990,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.0165, q=0, angle_radial=1.555806392, angle_azimuthal=1.570796327,
                                     angle_radial_out=1.555806392)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.011, x_right=0.011, y_bottom=-0.000375, y_top=0.000375)

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='CORRECTION !!!!  Generic Beam Screen/Slit/Stopper/Attenuator HSS',
                               boundary_shape=boundary_shape,
                               i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                               i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0.05, angle_radial=0, angle_azimuthal=4.71238898,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)


    return beam

# MONO BEFORE
def run_config_c():
    global SEED
    print(">>>>  saved seed: ", SEED)
    SEED += 555
    print(">>>>  using seed: ", SEED)
    ##################
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

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
        flag_autoset_flux_central_cone=0,  # value to set the flux peak
        flux_central_cone=1350000000000000.0,  # value to set the flux peak
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
    light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam,
                                                  magnetic_structure=source, nrays=150000, seed=SEED)
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
    coordinates = ElementCoordinates(p=2.2, q=0, angle_radial=1.53937691, angle_azimuthal=4.71238898,
                                     angle_radial_out=1.53937691)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
    beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                            movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
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
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
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
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
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
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
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
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=3.141592654,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.0065, x_right=0.0065, y_bottom=-0.06, y_top=0.06)

    from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirror
    optical_element = S4SphereMirror(name='Sphere Mirror', boundary_shape=boundary_shape,
                                     surface_calculation=0, is_cylinder=1, cylinder_direction=0,
                                     convexity=1, radius=1.000000, p_focus=28.300000, q_focus=11.700000,
                                     grazing_angle=0.031419,
                                     f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                     coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=2, q=22.5, angle_radial=1.53937691, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.53937691)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirrorElement
    beamline_element = S4SphereMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=1.570796327,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-1.4595e-05, x_right=1.4595e-05, y_bottom=-0.5, y_top=0.5)

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

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.03, y_top=0.03)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.014990,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=145.9, q=0.0335, angle_radial=1.555806392, angle_azimuthal=0,
                                     angle_radial_out=1.555806392)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.012537, y_top=0.012537)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=145.950000, q_focus=0.050000, grazing_angle=0.014990,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.0165, q=0, angle_radial=1.555806392, angle_azimuthal=1.570796327,
                                     angle_radial_out=1.555806392)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.011, x_right=0.011, y_bottom=-0.000375, y_top=0.000375)

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='CORRECTION !!!!  Generic Beam Screen/Slit/Stopper/Attenuator HSS',
                               boundary_shape=boundary_shape,
                               i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                               i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0.05, angle_radial=0, angle_azimuthal=4.71238898,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)



    return beam

# if __name__ == "__main__":
if True:
    import numpy
    from srxraylib.plot.gol import plot


    N = 50

    #
    # a
    #
    SEED = 1234567
    INTENSITY_a = []
    FWHM_a = []
    FWHMH_a = []
    FWHMV_a = []
    for i in range(N):
        beam = run_config_a()
        tkt = beam.histo1(26, xrange=[16997,17003], nbins=100, nolost=1, ref=23,
                   write=None, factor=1.0, calculate_widths=1, calculate_hew=0)
        tkt_h = beam.histo1(1, xrange=[-15e-9,15e-9], nbins=100, nolost=1, ref=23,
                   write=None, factor=1.0, calculate_widths=1, calculate_hew=0)
        tkt_v = beam.histo1(3, xrange=[-15e-9,15e-9], nbins=100, nolost=1, ref=23,
                   write=None, factor=1.0, calculate_widths=1, calculate_hew=0)
        # plot(tkt["bin_path"], tkt["histogram_path"])
        print("Intensity: ", beam.get_intensity(nolost=1))
        print("FWHM: eV, H, V:", tkt["fwhm"], tkt_h["fwhm"], tkt_v["fwhm"],)
        INTENSITY_a.append(beam.get_intensity(nolost=1))
        FWHM_a.append(tkt["fwhm"])
        FWHMH_a.append(tkt_h["fwhm"])
        FWHMV_a.append(tkt_v["fwhm"])

    INTENSITY_a = numpy.array(INTENSITY_a)
    FWHM_a = numpy.array(FWHM_a)
    FWHMH_a = 1e9 * numpy.array(FWHMH_a)
    FWHMV_a = 1e9 * numpy.array(FWHMV_a)


    print(">>>>>>>>>>>>>>>> CONFIG A <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_a, FWHM_a)
    print("Intensity: ", INTENSITY_a.mean(), " +/- ", numpy.std(INTENSITY_a))
    print("Energy FWHM: ", FWHM_a.mean(), " +/- ", numpy.std(FWHM_a))
    print("Energy FWHM H: ", FWHMH_a.mean(), " +/- ", numpy.std(FWHMH_a))
    print("Energy FWHM V: ", FWHMV_a.mean(), " +/- ", numpy.std(FWHMV_a))

    #
    # b
    #
    SEED = 1234567
    INTENSITY_b = []
    FWHM_b = []
    FWHMH_b = []
    FWHMV_b = []
    for i in range(N):
        beam = run_config_b()
        tkt = beam.histo1(26, xrange=[16997,17003], nbins=100, nolost=1, ref=23,
                   write=None, factor=1.0, calculate_widths=1, calculate_hew=0)
        tkt_h = beam.histo1(1, xrange=[-15e-9,15e-9], nbins=100, nolost=1, ref=23,
                   write=None, factor=1.0, calculate_widths=1, calculate_hew=0)
        tkt_v = beam.histo1(3, xrange=[-15e-9,15e-9], nbins=100, nolost=1, ref=23,
                   write=None, factor=1.0, calculate_widths=1, calculate_hew=0)
        # plot(tkt["bin_path"], tkt["histogram_path"])
        print("Intensity: ", beam.get_intensity(nolost=1))
        print("FWHM: eV, H, V:", tkt["fwhm"], tkt_h["fwhm"], tkt_v["fwhm"],)
        INTENSITY_b.append(beam.get_intensity(nolost=1))
        FWHM_b.append(tkt["fwhm"])
        FWHMH_b.append(tkt_h["fwhm"])
        FWHMV_b.append(tkt_v["fwhm"])

    INTENSITY_b = numpy.array(INTENSITY_b)
    FWHM_b = numpy.array(FWHM_b)
    FWHMH_b = 1e9 * numpy.array(FWHMH_b)
    FWHMV_b = 1e9 * numpy.array(FWHMV_b)


    print(">>>>>>>>>>>>>>>> CONFIG B <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_b, FWHM_b)
    print("Intensity: ", INTENSITY_b.mean(), " +/- ", numpy.std(INTENSITY_b))
    print("Energy FWHM: ", FWHM_b.mean(), " +/- ", numpy.std(FWHM_b))
    print("Energy FWHM H: ", FWHMH_b.mean(), " +/- ", numpy.std(FWHMH_b))
    print("Energy FWHM V: ", FWHMV_b.mean(), " +/- ", numpy.std(FWHMV_b))


    #
    # c
    #
    SEED = 1234567
    INTENSITY_c = []
    FWHM_c = []
    FWHMH_c = []
    FWHMV_c = []
    for i in range(N):
        beam = run_config_c()
        tkt = beam.histo1(26, xrange=[16997,17003], nbins=100, nolost=1, ref=23,
                   write=None, factor=1.0, calculate_widths=1, calculate_hew=0)
        tkt_h = beam.histo1(1, xrange=[-15e-9,15e-9], nbins=100, nolost=1, ref=23,
                   write=None, factor=1.0, calculate_widths=1, calculate_hew=0)
        tkt_v = beam.histo1(3, xrange=[-15e-9,15e-9], nbins=100, nolost=1, ref=23,
                   write=None, factor=1.0, calculate_widths=1, calculate_hew=0)
        # plot(tkt["bin_path"], tkt["histogram_path"])
        print("Intensity: ", beam.get_intensity(nolost=1))
        print("FWHM: eV, H, V:", tkt["fwhm"], tkt_h["fwhm"], tkt_v["fwhm"],)
        INTENSITY_c.append(beam.get_intensity(nolost=1))
        FWHM_c.append(tkt["fwhm"])
        FWHMH_c.append(tkt_h["fwhm"])
        FWHMV_c.append(tkt_v["fwhm"])

    INTENSITY_c = numpy.array(INTENSITY_c)
    FWHM_c = numpy.array(FWHM_c)
    FWHMH_c = 1e9 * numpy.array(FWHMH_c)
    FWHMV_c = 1e9 * numpy.array(FWHMV_c)


    print(">>>>>>>>>>>>>>>> CONFIG C <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_c, FWHM_c)
    print("Intensity: ", INTENSITY_c.mean(), " +/- ", numpy.std(INTENSITY_c))
    print("Energy FWHM: ", FWHM_c.mean(), " +/- ", numpy.std(FWHM_c))
    print("Energy FWHM H: ", FWHMH_c.mean(), " +/- ", numpy.std(FWHMH_c))
    print("Energy FWHM V: ", FWHMV_c.mean(), " +/- ", numpy.std(FWHMV_c))

#################

    print(">>>>>>>>>>>>>>>> CONFIG A <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_a, FWHM_a)
    print("Intensity: ", INTENSITY_a.mean(), " +/- ", numpy.std(INTENSITY_a))
    print("Energy FWHM: ", FWHM_a.mean(), " +/- ", numpy.std(FWHM_a))
    print("Energy FWHM H: ", FWHMH_a.mean(), " +/- ", numpy.std(FWHMH_a))
    print("Energy FWHM V: ", FWHMV_a.mean(), " +/- ", numpy.std(FWHMV_a))


    print(">>>>>>>>>>>>>>>> CONFIG B <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_b, FWHM_b)
    print("Intensity: ", INTENSITY_b.mean(), " +/- ", numpy.std(INTENSITY_b))
    print("Energy FWHM: ", FWHM_b.mean(), " +/- ", numpy.std(FWHM_b))
    print("Energy FWHM H: ", FWHMH_b.mean(), " +/- ", numpy.std(FWHMH_b))
    print("Energy FWHM V: ", FWHMV_b.mean(), " +/- ", numpy.std(FWHMV_b))

    print(">>>>>>>>>>>>>>>> CONFIG C <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_c, FWHM_c)
    print("Intensity: ", INTENSITY_c.mean(), " +/- ", numpy.std(INTENSITY_c))
    print("Energy FWHM: ", FWHM_c.mean(), " +/- ", numpy.std(FWHM_c))
    print("Energy FWHM H: ", FWHMH_c.mean(), " +/- ", numpy.std(FWHMH_c))
    print("Energy FWHM V: ", FWHMV_c.mean(), " +/- ", numpy.std(FWHMV_c))

"""

>>>>>>>>>>>>>>>> CONFIG A <<<<<<<<<<<<<<<<<<<<<<<
[1538. 1590. 1617. 1568. 1528. 1549. 1569. 1570. 1565. 1545. 1641. 1571.
 1543. 1607. 1563. 1604. 1650. 1524. 1551. 1596. 1547. 1550. 1583. 1634.
 1639. 1583. 1560. 1568. 1623. 1626. 1634. 1551. 1520. 1564. 1559. 1569.
 1534. 1583. 1557. 1564. 1556. 1534. 1550. 1508. 1509. 1576. 1565. 1556.
 1633. 1555.] [5.94 5.88 5.88 5.82 5.76 5.88 5.88 5.94 5.94 5.82 5.94 5.82 5.88 5.94
 5.76 5.82 5.94 5.82 5.82 5.7  5.94 5.94 5.88 5.94 5.88 5.76 5.94 5.88
 5.94 5.76 5.58 5.82 5.94 5.94 5.58 5.82 5.94 5.94 5.88 5.1  5.88 5.64
 5.82 5.64 5.94 5.88 5.88 5.88 5.82 5.76]
Intensity:  1571.579999999999  +/-  35.53369668357066
Energy FWHM:  5.834400000127353  +/-  0.14125381411088003
Energy FWHM H:  8.249999999999984  +/-  0.5977457653551369
Energy FWHM V:  6.305999999999988  +/-  0.6074240693288326
>>>>>>>>>>>>>>>> CONFIG B <<<<<<<<<<<<<<<<<<<<<<<
[383.29296191 357.16811853 393.31564604 351.23447984 341.61478137
 383.35968297 356.56376689 379.21380906 374.62338754 385.16389975
 357.32963466 370.10552701 344.28932848 340.77947982 380.87804416
 399.29067104 399.13198584 371.88249158 363.80489816 360.2881192
 366.08230387 377.66197988 371.19376014 369.63783551 369.18902849
 384.5279964  349.541421   373.45321029 376.66018708 382.77807145
 375.71984634 386.42554926 348.08080928 348.64553008 366.75237625
 368.69634268 348.13045322 381.16219339 370.83909504 396.88743935
 354.57518989 341.13793915 359.45117658 374.99201034 357.3199088
 384.71487307 380.22870847 362.77501166 391.42610873 386.5529552 ] [1.62 1.74 1.8  1.92 1.86 1.8  1.92 1.86 1.92 1.74 1.86 1.68 1.74 1.8
 1.86 1.62 1.92 1.92 1.92 1.86 1.86 1.8  1.98 1.92 1.74 1.8  1.86 1.8
 1.74 1.92 1.44 1.86 1.86 1.56 1.38 1.86 1.98 1.86 1.98 1.74 1.8  2.04
 1.8  1.62 1.68 1.5  1.98 1.62 1.8  1.86]
Intensity:  369.9714004954837  +/-  15.717627499392416
Energy FWHM:  1.8000000000392902  +/-  0.13994284547934185
Energy FWHM H:  7.865999999999985  +/-  0.9549052308999025
Energy FWHM V:  6.107999999999987  +/-  0.989108689679752
>>>>>>>>>>>>>>>> CONFIG C <<<<<<<<<<<<<<<<<<<<<<<
[374.36323658 345.15566932 374.96300978 341.32493122 326.63158642
 368.79478817 342.82290209 362.58467617 356.86511699 370.68104863
 348.00985895 353.19101452 329.13943356 330.49049317 360.98430701
 383.95399508 378.95343262 354.40325578 342.84392254 346.91622179
 350.34123749 365.69231713 355.14013565 353.81584602 354.38822289
 363.55830144 333.4426179  364.10766629 360.45671335 366.18464114
 359.3779064  367.9684847  332.57347616 336.87432416 352.16577506
 352.77015869 326.63594113 368.75103264 355.72469679 377.7180158
 339.02179692 330.36491129 339.97118609 357.53861537 347.26496895
 369.56731129 367.04746237 346.8168378  372.6657591  375.23337202] [1.62 1.68 1.68 1.92 1.68 1.8  1.74 1.68 1.92 1.62 1.86 1.5  1.68 1.68
 1.68 1.44 1.8  1.86 1.68 1.86 1.86 1.8  1.86 1.68 1.74 1.56 1.8  1.8
 1.68 1.74 1.26 1.62 1.8  1.56 1.38 1.68 1.62 1.8  1.8  1.74 1.8  2.04
 1.68 1.62 1.68 1.5  1.98 1.62 1.74 1.74]
Intensity:  354.7250526480217  +/-  14.995103973043207
Energy FWHM:  1.711200000037352  +/-  0.14411995004476916
Energy FWHM H:  8.051999999999984  +/-  1.0772631990372619
Energy FWHM V:  6.119999999999988  +/-  1.034021276376843
Thu 13 Feb 2025 06:19:24 PM CET


"""