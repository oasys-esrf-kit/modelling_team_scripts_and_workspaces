

def run_config_2xtals():
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
        photon_energy=7000.0,  # Photon energy (in eV)
        delta_e=2.0,  # Photon energy width (in eV)
        ng_e=100,  # Photon energy scan number of points
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread=1,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
        harmonic_number=1,  # harmonic number
        flag_autoset_flux_central_cone=1,  # value to set the flux peak
        flux_central_cone=3228083726728806.0,  # value to set the flux peak
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
                                     f_central=1, f_phot_cent=0, phot_cent=7000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.284411084, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.284411084)
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
                                     f_central=1, f_phot_cent=0, phot_cent=7000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.284411084, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.284411084)
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
                                     convexity=1, radius=1.000000, p_focus=33.500000, q_focus=11.100000,
                                     grazing_angle=0.031419,
                                     f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                     coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=4, q=11.1, angle_radial=1.53937691, angle_azimuthal=3.141592654,
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

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.03, y_top=0.03)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.036400,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=155.3, q=0.0335, angle_radial=1.534396327, angle_azimuthal=0,
                                     angle_radial_out=1.534396327)
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
                                        p_focus=155.350000, q_focus=0.050000, grazing_angle=0.036400,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.0165, q=0, angle_radial=1.534396327, angle_azimuthal=1.570796327,
                                     angle_radial_out=1.534396327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.011, x_right=0.011, y_bottom=-0.0009125, y_top=0.0009125)

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

# 4 crystals
def run_config_4xtals():
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
        photon_energy=7000.0,  # Photon energy (in eV)
        delta_e=2.0,  # Photon energy width (in eV)
        ng_e=100,  # Photon energy scan number of points
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread=1,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
        harmonic_number=1,  # harmonic number
        flag_autoset_flux_central_cone=1,  # value to set the flux peak
        flux_central_cone=3228083726728806.0,  # value to set the flux peak
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
                                     f_central=1, f_phot_cent=0, phot_cent=7000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.284411084, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.284411084)
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
                                     f_central=1, f_phot_cent=0, phot_cent=7000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.284411084, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.284411084)
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
                                     f_central=1, f_phot_cent=0, phot_cent=7000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.284411084, angle_azimuthal=0,
                                     angle_radial_out=1.284411084)
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
                                     f_central=1, f_phot_cent=0, phot_cent=7000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.284411084, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.284411084)
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
                                     convexity=1, radius=1.000000, p_focus=33.500000, q_focus=11.100000,
                                     grazing_angle=0.031419,
                                     f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                     coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=4, q=11.1, angle_radial=1.53937691, angle_azimuthal=3.141592654,
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

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.03, y_top=0.03)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.036400,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=155.3, q=0.0335, angle_radial=1.534396327, angle_azimuthal=0,
                                     angle_radial_out=1.534396327)
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
                                        p_focus=155.350000, q_focus=0.050000, grazing_angle=0.036400,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.0165, q=0, angle_radial=1.534396327, angle_azimuthal=1.570796327,
                                     angle_radial_out=1.534396327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.011, x_right=0.011, y_bottom=-0.0009125, y_top=0.0009125)

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
    # 2 crystals
    #
    SEED = 1234567
    INTENSITY_a = []
    FWHM_a = []
    FWHMH_a = []
    FWHMV_a = []
    for i in range(N):
        beam = run_config_2xtals()
        tkt = beam.histo1(26, xrange=[6999,7001], nbins=100, nolost=1, ref=23,
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


    print(">>>>>>>>>>>>>>>> CONFIG 2 crystals <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_a, FWHM_a)
    print("Intensity: ", INTENSITY_a.mean(), " +/- ", numpy.std(INTENSITY_a))
    print("Energy FWHM: ", FWHM_a.mean(), " +/- ", numpy.std(FWHM_a))
    print("Energy FWHM H: ", FWHMH_a.mean(), " +/- ", numpy.std(FWHMH_a))
    print("Energy FWHM V: ", FWHMV_a.mean(), " +/- ", numpy.std(FWHMV_a))

    #
    # 4 crystals
    #
    SEED = 1234567
    INTENSITY_b = []
    FWHM_b = []
    FWHMH_b = []
    FWHMV_b = []
    for i in range(N):
        beam = run_config_4xtals()
        tkt = beam.histo1(26, xrange=[6999,7001], nbins=100, nolost=1, ref=23,
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


    print(">>>>>>>>>>>>>>>> CONFIG 4 crystals <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_b, FWHM_b)
    print("Intensity: ", INTENSITY_b.mean(), " +/- ", numpy.std(INTENSITY_b))
    print("Energy FWHM: ", FWHM_b.mean(), " +/- ", numpy.std(FWHM_b))
    print("Energy FWHM H: ", FWHMH_b.mean(), " +/- ", numpy.std(FWHMH_b))
    print("Energy FWHM V: ", FWHMV_b.mean(), " +/- ", numpy.std(FWHMV_b))



#################

    print(">>>>>>>>>>>>>>>> CONFIG 2 crystals <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_a, FWHM_a)
    print("Intensity: ", INTENSITY_a.mean(), " +/- ", numpy.std(INTENSITY_a))
    print("Energy FWHM: ", FWHM_a.mean(), " +/- ", numpy.std(FWHM_a))
    print("Energy FWHM H: ", FWHMH_a.mean(), " +/- ", numpy.std(FWHMH_a))
    print("Energy FWHM V: ", FWHMV_a.mean(), " +/- ", numpy.std(FWHMV_a))


    print(">>>>>>>>>>>>>>>> CONFIG 4 crystals <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_b, FWHM_b)
    print("Intensity: ", INTENSITY_b.mean(), " +/- ", numpy.std(INTENSITY_b))
    print("Energy FWHM: ", FWHM_b.mean(), " +/- ", numpy.std(FWHM_b))
    print("Energy FWHM H: ", FWHMH_b.mean(), " +/- ", numpy.std(FWHMH_b))
    print("Energy FWHM V: ", FWHMV_b.mean(), " +/- ", numpy.std(FWHMV_b))


"""
>>>>>>>>>>>>>>>> CONFIG 2 crystals <<<<<<<<<<<<<<<<<<<<<<<
[2746.92707746 2602.84572306 2725.25350102 2623.66281943 2628.52793708
 2615.05634617 2665.84173483 2654.8996056  2714.98405037 2610.71918422
 2632.70929054 2586.2232694  2700.23620078 2674.07567911 2661.60342285
 2707.54276827 2708.37696733 2623.55461887 2651.80290051 2709.19613194
 2657.71329283 2672.96046281 2654.70410313 2719.24414589 2683.7579475
 2636.15971753 2586.18013848 2691.82958529 2646.03385488 2694.19088677
 2732.12753594 2661.30655403 2714.32079861 2624.38549352 2671.90218445
 2634.74239697 2679.89119753 2658.71626101 2666.25638138 2679.98152638
 2625.21225571 2625.86556766 2654.3652813  2674.84120257 2718.91537435
 2626.72752993 2693.77114508 2653.65747181 2696.94414366 2615.5521792 ] [0.68 0.72 0.72 0.7  0.72 0.72 0.74 0.74 0.72 0.72 0.7  0.72 0.7  0.68
 0.72 0.72 0.7  0.7  0.72 0.7  0.72 0.72 0.68 0.72 0.72 0.66 0.72 0.7
 0.74 0.72 0.72 0.7  0.74 0.7  0.68 0.74 0.68 0.7  0.72 0.7  0.72 0.72
 0.7  0.68 0.68 0.72 0.72 0.68 0.72 0.7 ]
Intensity:  2663.9259169003153  +/-  38.99040120922881
Energy FWHM:  0.7092000000154803  +/-  0.01885099466912856
Energy FWHM H:  6.0899999999999865  +/-  0.4695742752749548
Energy FWHM V:  6.953999999999985  +/-  0.5150572783681432
>>>>>>>>>>>>>>>> CONFIG 4 crystals <<<<<<<<<<<<<<<<<<<<<<<
[1831.87986166 1733.65618228 1821.33610157 1746.76595349 1741.13908632
 1746.46465047 1765.56871579 1762.4628669  1797.28700369 1742.91111087
 1744.85232357 1723.25392232 1799.62394438 1780.19113672 1767.16669288
 1803.22428334 1800.82465626 1740.86940055 1751.63673183 1805.59351878
 1760.95744956 1776.56086179 1770.70587116 1809.39247118 1776.30549429
 1758.71857036 1715.19812738 1791.7491608  1753.584262   1793.49514887
 1817.98880046 1772.15449832 1800.16961607 1745.01570952 1777.69380628
 1741.4532883  1780.78994152 1766.43591238 1767.60815982 1786.25575777
 1755.56662886 1746.14095203 1770.18964098 1782.46738072 1805.33765719
 1757.59450845 1781.42277688 1760.89383501 1788.23171261 1727.91763738] [0.6  0.6  0.62 0.64 0.62 0.62 0.62 0.64 0.64 0.6  0.64 0.6  0.6  0.6
 0.58 0.62 0.6  0.58 0.6  0.6  0.64 0.66 0.62 0.66 0.62 0.58 0.6  0.56
 0.6  0.6  0.62 0.6  0.64 0.64 0.58 0.66 0.6  0.56 0.64 0.58 0.58 0.62
 0.58 0.58 0.58 0.62 0.64 0.64 0.58 0.6 ]
Intensity:  1770.8940756317697  +/-  26.627129778157858
Energy FWHM:  0.610000000013315  +/-  0.025690465157891026
Energy FWHM H:  5.999999999999987  +/-  0.5055689863905805
Energy FWHM V:  6.863999999999987  +/-  0.48979995916700425


"""

