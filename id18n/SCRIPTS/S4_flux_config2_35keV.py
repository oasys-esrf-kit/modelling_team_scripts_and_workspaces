

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
        photon_energy=35000.0,  # Photon energy (in eV)
        delta_e=16.0,  # Photon energy width (in eV)
        ng_e=100,  # Photon energy scan number of points
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread=1,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
        harmonic_number=5,  # harmonic number
        flag_autoset_flux_central_cone=1,  # value to set the flux peak
        flux_central_cone=751468106877072.0,  # value to set the flux peak
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
                                     f_central=1, f_phot_cent=0, phot_cent=35000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.514269058, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.514269058)
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
                                     f_central=1, f_phot_cent=0, phot_cent=35000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.514269058, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.514269058)
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
                                        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.007290,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=155.3, q=0.0335, angle_radial=1.563506327, angle_azimuthal=0,
                                     angle_radial_out=1.563506327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.012, x_right=0.012, y_bottom=-0.0004362, y_top=0.0004362)

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
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.012537, y_top=0.012537)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=155.350000, q_focus=0.050000, grazing_angle=0.007290,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.0165, q=0, angle_radial=1.563506327, angle_azimuthal=1.570796327,
                                     angle_radial_out=1.563506327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.011, x_right=0.011, y_bottom=-0.000183, y_top=0.000183)

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
        photon_energy=35000.0,  # Photon energy (in eV)
        delta_e=16.0,  # Photon energy width (in eV)
        ng_e=100,  # Photon energy scan number of points
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread=1,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
        harmonic_number=5,  # harmonic number
        flag_autoset_flux_central_cone=1,  # value to set the flux peak
        flux_central_cone=751468106877072.0,  # value to set the flux peak
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
                                     f_central=1, f_phot_cent=0, phot_cent=35000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.514269058, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.514269058)
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
                                     f_central=1, f_phot_cent=0, phot_cent=35000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.514269058, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.514269058)
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
                                     f_central=1, f_phot_cent=0, phot_cent=35000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.514269058, angle_azimuthal=0,
                                     angle_radial_out=1.514269058)
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
                                     f_central=1, f_phot_cent=0, phot_cent=35000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.514269058, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.514269058)
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
                                        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.007290,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=155.3, q=0.0335, angle_radial=1.563506327, angle_azimuthal=0,
                                     angle_radial_out=1.563506327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.012, x_right=0.012, y_bottom=-0.0004362, y_top=0.0004362)

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
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.012537, y_top=0.012537)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                        surface_calculation=0,
                                        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                        p_focus=155.350000, q_focus=0.050000, grazing_angle=0.007290,
                                        is_cylinder=1, cylinder_direction=0, convexity=1,
                                        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.0165, q=0, angle_radial=1.563506327, angle_azimuthal=1.570796327,
                                     angle_radial_out=1.563506327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.011, x_right=0.011, y_bottom=-0.000183, y_top=0.000183)

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
        tkt = beam.histo1(26, xrange=[34995,35005], nbins=100, nolost=1, ref=23,
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
        tkt = beam.histo1(26, xrange=[34995,35005], nbins=100, nolost=1, ref=23,
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


    print(">>>>>>>>>>>>>>>> CONFIG 2 crystals <<<<<<<<<<<<<<<<<<<<<<<")
    print(INTENSITY_a, FWHM_a)
    print("Intensity: ", INTENSITY_a.mean(), " +/- ", numpy.std(INTENSITY_a))
    print("Energy FWHM: ", FWHM_a.mean(), " +/- ", numpy.std(FWHM_a))
    print("Energy FWHM H: ", FWHMH_a.mean(), " +/- ", numpy.std(FWHMH_a))
    print("Energy FWHM V: ", FWHMV_a.mean(), " +/- ", numpy.std(FWHMV_a))


#################

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
[186.90313691 170.83388445 150.86797896 170.20022383 165.80812472
 168.31850198 182.60908422 167.02152444 172.39993784 173.62605794
 169.34366107 181.97164693 171.54998302 157.84585909 173.13473232
 197.42464157 210.32851736 175.85518893 150.3177271  170.49551572
 175.25891532 185.33663086 165.52905694 164.41534635 170.44127493
 177.13999884 182.05834833 188.53953555 179.49521527 172.99911178
 184.66945413 183.55484614 163.77952525 162.30636042 175.26929914
 171.90857814 174.60907545 181.13858385 171.30960205 179.43126377
 162.90978732 165.59264993 166.60074138 155.75529054 162.10950138
 178.44908161 184.93678503 179.90937092 193.12471579 189.96483205] [2.8 2.6 3.7 3.7 2.8 3.1 4.2 4.  3.7 3.9 4.  1.8 4.  3.6 5.1 3.6 4.2 2.4
 1.3 2.2 4.3 4.7 4.6 4.  2.5 4.1 3.  4.2 3.9 4.1 2.1 3.5 2.8 3.  4.5 3.8
 4.3 4.1 4.3 3.2 1.5 5.3 2.7 2.9 4.4 3.8 3.6 3.5 4.2 3.6]
Intensity: 174 +/- 11
Energy FWHM: 3.54 +/- 0.88
FWHM HxV: 5.7 x 5.2 
FWHM H:  5.669999999999988  +/-  1.4930840565755135
FWHM V:  5.201999999999989  +/-  1.3139999999999972
>>>>>>>>>>>>>>>> CONFIG 4 crystals <<<<<<<<<<<<<<<<<<<<<<<
[149.39146108 122.79387653 117.19480471 125.1696552  122.32292081
 130.39530431 134.49125243 122.81737046 130.98339576 137.25630181
 127.81976999 139.43527203 116.52248867 116.28598313 129.96894774
 149.28574661 161.15802201 125.17919851 122.36999243 127.37206779
 131.10748524 142.29554447 119.19975344 117.62693401 125.68010423
 134.41143492 143.74338356 136.81269209 141.85693329 130.87916962
 146.03498406 147.10786717 126.72904164 118.94434762 120.11976785
 124.1594625  122.23408014 139.69324868 127.10352318 136.10609839
 125.39485185 113.17464616 126.28921121 114.33781826 117.36412507
 140.12458148 135.84259124 131.11915181 140.4499265  143.72790335] [2.8 2.3 3.4 3.7 0.2 2.1 3.8 2.4 2.9 2.4 3.2 1.8 3.1 3.1 3.7 2.1 3.2 1.1
 1.3 2.2 3.2 3.5 2.8 3.2 2.5 3.7 2.5 3.  3.3 3.2 0.7 3.1 2.4 2.6 1.7 2.9
 1.8 3.6 3.6 2.6 2.9 4.1 2.7 2.9 2.8 2.9 3.1 2.9 2.8 3. ]
Intensity: 130 +/- 10
Energy FWHM: 2.74 +/- 0.78
FWHM HxV: 4.1 x 5.0 
Energy FWHM H:  4.127999999999991  +/-  1.403999999999997
Energy FWHM V:  4.97399999999999  +/-  1.3853245107194168


"""