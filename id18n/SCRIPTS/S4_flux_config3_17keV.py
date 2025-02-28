

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


    print(">>>>>>>>>>>>>>>> CONFIG 4 crystals <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_b, FWHM_b)
    print("Intensity: ", INTENSITY_b.mean(), " +/- ", numpy.std(INTENSITY_b))
    print("Energy FWHM: ", FWHM_b.mean(), " +/- ", numpy.std(FWHM_b))
    print("Energy FWHM H: ", FWHMH_b.mean(), " +/- ", numpy.std(FWHMH_b))
    print("Energy FWHM V: ", FWHMV_b.mean(), " +/- ", numpy.std(FWHMV_b))



#################

    print(">>>>>>>>>>>>>>>> CONFIG 2 crystals <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_a, FWHM_a)
    print("Intensity: %d +/- %d" %(INTENSITY_a.mean(), numpy.std(INTENSITY_a)))
    print("Energy FWHM: %.2f +/- %.2f" %  (FWHM_a.mean(), numpy.std(FWHM_a)))
    print("FWHM HxV: %.1f x %.1f " % (FWHMH_a.mean(), FWHMV_a.mean()))
    print("FWHM H: ", FWHMH_a.mean(), " +/- ", numpy.std(FWHMH_a))
    print("FWHM V: ", FWHMV_a.mean(), " +/- ", numpy.std(FWHMV_a))


    print(">>>>>>>>>>>>>>>> CONFIG 4 crystals <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_b, FWHM_b)
    print("Intensity: %d +/- %d" %(INTENSITY_b.mean(), numpy.std(INTENSITY_b)))
    print("Energy FWHM: %.2f +/- %.2f" %  (FWHM_b.mean(), numpy.std(FWHM_b)))
    print("FWHM HxV: %.1f x %.1f " % (FWHMH_b.mean(), FWHMV_b.mean()))
    print("Energy FWHM H: ", FWHMH_b.mean(), " +/- ", numpy.std(FWHMH_b))
    print("Energy FWHM V: ", FWHMV_b.mean(), " +/- ", numpy.std(FWHMV_b))

"""

>>>>>>>>>>>>>>>> CONFIG 2 crystals <<<<<<<<<<<<<<<<<<<<<<<
[639.19054682 650.31827302 672.16277108 643.57643366 621.3340832
 646.82453892 637.68872186 668.34169716 683.47328464 677.71492135
 666.61881963 630.33236102 643.02583085 618.80540987 638.76353692
 674.0043619  653.40709777 666.09460224 651.80010943 680.59446671
 635.88359534 645.22422575 657.98816662 675.8434942  680.16485473
 652.5784671  647.12845098 681.25609544 705.68774701 714.51656204
 611.12654541 662.07800252 627.74790043 619.72849073 643.00819576
 654.42361255 664.78803764 668.73800397 668.18319591 692.89259166
 639.67739187 635.00070667 622.35522809 646.26654301 647.90518339
 682.15702988 640.19560539 632.46985155 653.44767479 627.85460036] [1.98 1.98 1.98 2.16 2.16 1.92 2.04 1.86 2.04 1.8  2.16 1.98 2.1  2.1
 1.98 2.04 2.1  2.16 1.92 1.86 2.04 1.92 1.8  2.04 2.1  1.8  2.04 2.04
 2.22 2.04 1.92 2.04 1.86 2.04 1.86 2.1  2.04 1.8  1.98 1.92 2.1  2.16
 1.86 2.16 1.86 1.92 2.1  2.1  1.74 1.92]
Intensity: 654 +/- 22
Energy FWHM: 2.00 +/- 0.12
FWHM HxV: 8.1 x 6.3 
FWHM H:  8.123999999999985  +/-  0.6230762393158629
FWHM V:  6.281999999999987  +/-  0.7646410922779378
>>>>>>>>>>>>>>>> CONFIG 4 crystals <<<<<<<<<<<<<<<<<<<<<<<
[542.41088411 552.72549498 578.34505653 544.65737324 516.74504354
 553.66039151 546.12066495 575.1615148  581.49999318 582.12058164
 558.78436589 545.77778435 546.37456852 522.13328987 535.40900602
 577.33093625 556.98084097 571.53668222 550.03489255 583.05763686
 540.13611266 556.16957457 556.10757028 564.76676788 583.4016778
 556.11709069 551.57380133 582.30587677 588.0173374  613.4355391
 522.20398349 567.75697464 530.11558689 526.48122659 553.15508073
 549.06082622 563.68886446 571.71137916 576.22904005 593.73987219
 546.67148922 530.77241349 531.5814362  555.67565103 554.95499772
 583.58528923 541.0483651  531.96475574 557.25787552 532.17649813] [1.86 1.74 1.86 1.86 1.92 1.86 1.86 1.86 1.92 1.68 1.92 1.74 1.98 1.68
 1.86 1.8  1.8  1.92 1.8  1.74 1.86 1.86 1.68 1.92 1.86 1.8  1.86 1.8
 1.74 1.8  1.86 1.86 1.86 1.74 1.86 1.86 1.98 1.8  1.8  1.86 1.92 2.04
 1.86 1.92 1.8  1.8  1.86 1.86 1.74 1.86]
Intensity: 556 +/- 20
Energy FWHM: 1.84 +/- 0.08
FWHM HxV: 8.1 x 6.2 
Energy FWHM H:  8.075999999999985  +/-  0.5934846249061541
Energy FWHM V:  6.2099999999999875  +/-  0.9434511116109813


"""