

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
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.0065, x_right=0.0065, y_bottom=-0.06, y_top=0.06)

    from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirror
    optical_element = S4SphereMirror(name='Sphere Mirror', boundary_shape=boundary_shape,
                                     surface_calculation=0, is_cylinder=1, cylinder_direction=0,
                                     convexity=1, radius=1.000000, p_focus=31.500000, q_focus=22.500000,
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
                                        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.036400,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=145.9, q=0.0335, angle_radial=1.534396327, angle_azimuthal=0,
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
                                        p_focus=145.950000, q_focus=0.050000, grazing_angle=0.036400,
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
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.0065, x_right=0.0065, y_bottom=-0.06, y_top=0.06)

    from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirror
    optical_element = S4SphereMirror(name='Sphere Mirror', boundary_shape=boundary_shape,
                                     surface_calculation=0, is_cylinder=1, cylinder_direction=0,
                                     convexity=1, radius=1.000000, p_focus=31.500000, q_focus=22.500000,
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
                                        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.036400,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=145.9, q=0.0335, angle_radial=1.534396327, angle_azimuthal=0,
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
                                        p_focus=145.950000, q_focus=0.050000, grazing_angle=0.036400,
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
[6116.77262378 5810.16411107 6049.0668608  6000.34697576 5964.43592701
 5866.3994785  6034.50892935 6025.11993126 6015.59455828 6019.20434033
 6021.14138199 5977.85404603 5996.6295317  5950.25665747 5919.35535765
 5917.38959415 6041.94972295 6035.93571876 5960.98161901 6079.35524471
 6027.20805808 5978.62986515 6019.89170132 6023.64674408 5966.68405518
 6009.90563768 5879.91890233 6025.17524774 6014.57857773 6000.54324171
 6032.44814944 6026.15599354 6135.75154178 5975.47339591 6046.6077067
 6033.85736565 6028.51355444 5931.2149792  5956.23656591 5960.54539671
 5926.14405196 5879.06121133 5944.06720558 5976.10718127 6062.60240338
 5960.64337159 6062.17874067 5935.53899728 5939.56951344 5936.93526689] [0.7  0.72 0.72 0.7  0.72 0.72 0.74 0.72 0.78 0.72 0.74 0.72 0.72 0.68
 0.7  0.74 0.7  0.7  0.72 0.7  0.72 0.74 0.68 0.74 0.7  0.7  0.72 0.72
 0.76 0.68 0.76 0.72 0.74 0.72 0.74 0.74 0.72 0.72 0.72 0.72 0.7  0.7
 0.7  0.68 0.76 0.7  0.72 0.7  0.7  0.7 ]
Intensity:  5989.965944685132  +/-  62.103702227220296
Energy FWHM:  0.7176000000156637  +/-  0.021777052142581264
Energy FWHM H:  14.61599999999997  +/-  0.8400857099129809
Energy FWHM V:  7.1159999999999854  +/-  0.3601999444752865
>>>>>>>>>>>>>>>> CONFIG 4 crystals <<<<<<<<<<<<<<<<<<<<<<<
[3973.47872279 3772.04392144 3931.33318901 3905.38858004 3874.48824008
 3822.06692111 3929.81242347 3921.07162727 3886.8574084  3916.15722583
 3906.43039554 3893.78109937 3898.80181062 3887.11939826 3865.30498004
 3848.42700199 3928.28228449 3924.64091031 3859.87721311 3963.31854949
 3916.68375618 3886.89909094 3933.00610998 3920.72501784 3882.10216855
 3919.74191613 3822.27442108 3912.8777416  3886.39647375 3900.36060471
 3916.97850718 3911.86492864 3978.79740324 3885.20650044 3934.5437424
 3906.28501138 3922.69894604 3863.51175609 3873.3731375  3889.62390013
 3872.64078229 3814.28147018 3853.09517504 3892.78159309 3945.96555277
 3900.20996431 3962.47010486 3871.3304448  3879.54559099 3867.46652188] [0.6  0.58 0.62 0.62 0.56 0.58 0.62 0.62 0.62 0.62 0.62 0.62 0.6  0.58
 0.58 0.62 0.6  0.58 0.6  0.62 0.62 0.64 0.62 0.62 0.6  0.62 0.6  0.58
 0.64 0.58 0.64 0.6  0.62 0.6  0.6  0.66 0.6  0.58 0.6  0.58 0.58 0.58
 0.58 0.62 0.6  0.6  0.62 0.6  0.56 0.6 ]
Intensity:  3896.6484047333524  +/-  40.0838644826502
Energy FWHM:  0.604000000013184  +/-  0.02116601048897873
Energy FWHM H:  14.567999999999973  +/-  0.9102615008886165
Energy FWHM V:  7.037999999999986  +/-  0.40343029137634095


"""