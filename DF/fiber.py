import numpy as np


def run_beamline(rotation_x=0):
    import numpy as np
    from dabax.dabax_xraylib import DabaxXraylib
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    #
    #
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    light_source = SourceGeometrical(name='Geometrical Source', nrays=5000, seed=5676561)
    light_source.set_spatial_type_point()
    light_source.set_depth_distribution_off()
    light_source.set_angular_distribution_cone(cone_max=0.250000, cone_min=0.000000)
    light_source.set_energy_distribution_singleline(15550.000000, unit='A')
    light_source.set_polarization(polarization_degree=0.000000, phase_diff=0.008727, coherent_beam=1)
    beam = light_source.get_beam()

    beamline.set_light_source(light_source)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Screen from source to lens (pass 1)', boundary_shape=boundary_shape,
                               i_abs=0,  # attenuation: 0=No, 1=prerefl file, 2=xraylib, 3=dabax
                               i_stop=0,  # 0=slit or aperture, 1=beam stop
                               thick=0,  # for i_abs>0
                               file_abs='<specify file name>',  # for i_abs=1
                               material='Au', density=19.3,  # for i_abs=2,3
                               dabax=None,
                               # if using dabax (i_abs=3), instance of DabaxXraylib() (use None for default)
                               )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.02, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None
    from shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror import S4ParaboloidMirror
    optical_element = S4ParaboloidMirror(name='Parabolic Mirror 1 (pass 1)', boundary_shape=boundary_shape,
                                         at_infinity=1,
                                         surface_calculation=0,
                                         p_focus=0.020000, q_focus=0.060000, grazing_angle=np.pi/2,
                                         # for surface_calculation=0
                                         parabola_parameter=1.000000, pole_to_focus=1.000000,  # for external
                                         is_cylinder=0, cylinder_direction=0, convexity=1,
                                         f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                                         f_refl=6,
                                         # for f_reflec=1: file: 0=prerefl, 2=(mrad, refl), 3=(eV, refl), 4=(eV, mrad, refl); 1=refr index, 5=xraylib, 6=dabax
                                         file_refl='<none>',  # for f_refl=0,2,3,4
                                         refraction_index=0.99999 + 0.001j,  # for f_refl=1
                                         coating_material='Ni', coating_density=8.902,  # for f_refl=5,6
                                         coating_roughness=0,  # for f_refl=0,1,5,6
                                         dabax=None,
                                         # if using dabax (f_reflec=1,f_refl=6), instance of DabaxXraylib() (use None for default)
                                         )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=0)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror import S4ParaboloidMirrorElement
    beamline_element = S4ParaboloidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                 movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror
    optical_element = S4PlaneMirror(name='Plane Mirror (backreflecting pass1)', boundary_shape=boundary_shape,
                                    f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                                    f_refl=1,
                                    # for f_reflec=1: file: 0=prerefl, 2=(mrad, refl), 3=(eV, refl), 4=(eV, mrad, refl); 1=refr index, 5=xraylib, 6=dabax
                                    file_refl='<none>',  # for f_refl=0,2,3,4
                                    refraction_index=0 + 1.51j,  # for f_refl=1
                                    coating_material='Ni', coating_density=8.902,  # for f_refl=5,6
                                    coating_roughness=0,  # for f_refl=0,1,5,6
                                    dabax=None,
                                    # if using dabax (f_reflec=1,f_refl=6), instance of DabaxXraylib() (use None for default)
                                    )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.1, q=0.1, angle_radial=0, angle_azimuthal=0, angle_radial_out=0)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_x=rotation_x,
                                           rotation_y=0, rotation_z=0)
    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
    beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                            movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None
    from shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror import S4ParaboloidMirror
    optical_element = S4ParaboloidMirror(name='Parabolic Mirror 2 (pass 1)', boundary_shape=boundary_shape,
                                         at_infinity=0,
                                         surface_calculation=0,
                                         p_focus=0.020000, q_focus=0.020000, grazing_angle=np.pi/2,
                                         # for surface_calculation=0
                                         parabola_parameter=1.000000, pole_to_focus=1.000000,  # for external
                                         is_cylinder=0, cylinder_direction=0, convexity=1,
                                         f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                                         f_refl=6,
                                         # for f_reflec=1: file: 0=prerefl, 2=(mrad, refl), 3=(eV, refl), 4=(eV, mrad, refl); 1=refr index, 5=xraylib, 6=dabax
                                         file_refl='<none>',  # for f_refl=0,2,3,4
                                         refraction_index=0.99999 + 0.001j,  # for f_refl=1
                                         coating_material='Ni', coating_density=8.902,  # for f_refl=5,6
                                         coating_roughness=0,  # for f_refl=0,1,5,6
                                         dabax=None,
                                         # if using dabax (f_reflec=1,f_refl=6), instance of DabaxXraylib() (use None for default)
                                         )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=0)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror import S4ParaboloidMirrorElement
    beamline_element = S4ParaboloidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                 movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name=' Screen at source position (pass 1)', boundary_shape=boundary_shape,
                               i_abs=0,  # attenuation: 0=No, 1=prerefl file, 2=xraylib, 3=dabax
                               i_stop=0,  # 0=slit or aperture, 1=beam stop
                               thick=0,  # for i_abs>0
                               file_abs='<specify file name>',  # for i_abs=1
                               material='Au', density=19.3,  # for i_abs=2,3
                               dabax=None,
                               # if using dabax (i_abs=3), instance of DabaxXraylib() (use None for default)
                               )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.02, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Ellipse
    boundary_shape = Ellipse(a_axis_min=-4.5e-06, a_axis_max=4.5e-06, b_axis_min=-4.5e-06, b_axis_max=4.5e-06)

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Stopper', boundary_shape=boundary_shape,
                               i_abs=0,  # attenuation: 0=No, 1=prerefl file, 2=xraylib, 3=dabax
                               i_stop=1,  # 0=slit or aperture, 1=beam stop
                               thick=0,  # for i_abs>0
                               file_abs='<specify file name>',  # for i_abs=1
                               material='Au', density=19.3,  # for i_abs=2,3
                               dabax=None,
                               # if using dabax (i_abs=3), instance of DabaxXraylib() (use None for default)
                               )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()
    beam_stopper = beam.duplicate()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror
    optical_element = S4PlaneMirror(name='Plane mirror st source position (from pass 1 to pass 2)',
                                    boundary_shape=boundary_shape,
                                    f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                                    f_refl=1,
                                    # for f_reflec=1: file: 0=prerefl, 2=(mrad, refl), 3=(eV, refl), 4=(eV, mrad, refl); 1=refr index, 5=xraylib, 6=dabax
                                    file_refl='<none>',  # for f_refl=0,2,3,4
                                    refraction_index=0 + 1.51j,  # for f_refl=1
                                    coating_material='Ni', coating_density=8.902,  # for f_refl=5,6
                                    coating_roughness=0,  # for f_refl=0,1,5,6
                                    dabax=None,
                                    # if using dabax (f_reflec=1,f_refl=6), instance of DabaxXraylib() (use None for default)
                                    )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=0)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
    beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                            movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Screen from source to lens (pass 2)', boundary_shape=boundary_shape,
                               i_abs=0,  # attenuation: 0=No, 1=prerefl file, 2=xraylib, 3=dabax
                               i_stop=0,  # 0=slit or aperture, 1=beam stop
                               thick=0,  # for i_abs>0
                               file_abs='<specify file name>',  # for i_abs=1
                               material='Au', density=19.3,  # for i_abs=2,3
                               dabax=None,
                               # if using dabax (i_abs=3), instance of DabaxXraylib() (use None for default)
                               )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.02, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None
    from shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror import S4ParaboloidMirror
    optical_element = S4ParaboloidMirror(name='Parabolic Mirror 1 (pass 2)', boundary_shape=boundary_shape,
                                         at_infinity=1,
                                         surface_calculation=0,
                                         p_focus=0.020000, q_focus=0.060000, grazing_angle=np.pi/2,
                                         # for surface_calculation=0
                                         parabola_parameter=1.000000, pole_to_focus=1.000000,  # for external
                                         is_cylinder=0, cylinder_direction=0, convexity=1,
                                         f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                                         f_refl=6,
                                         # for f_reflec=1: file: 0=prerefl, 2=(mrad, refl), 3=(eV, refl), 4=(eV, mrad, refl); 1=refr index, 5=xraylib, 6=dabax
                                         file_refl='<none>',  # for f_refl=0,2,3,4
                                         refraction_index=0.99999 + 0.001j,  # for f_refl=1
                                         coating_material='Ni', coating_density=8.902,  # for f_refl=5,6
                                         coating_roughness=0,  # for f_refl=0,1,5,6
                                         dabax=None,
                                         # if using dabax (f_reflec=1,f_refl=6), instance of DabaxXraylib() (use None for default)
                                         )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=0)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror import S4ParaboloidMirrorElement
    beamline_element = S4ParaboloidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                 movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror
    optical_element = S4PlaneMirror(name='Plane Mirror (backreflecting pass2)', boundary_shape=boundary_shape,
                                    f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                                    f_refl=1,
                                    # for f_reflec=1: file: 0=prerefl, 2=(mrad, refl), 3=(eV, refl), 4=(eV, mrad, refl); 1=refr index, 5=xraylib, 6=dabax
                                    file_refl='<none>',  # for f_refl=0,2,3,4
                                    refraction_index=0 + 1.51j,  # for f_refl=1
                                    coating_material='Ni', coating_density=8.902,  # for f_refl=5,6
                                    coating_roughness=0,  # for f_refl=0,1,5,6
                                    dabax=None,
                                    # if using dabax (f_reflec=1,f_refl=6), instance of DabaxXraylib() (use None for default)
                                    )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.1, q=0.1, angle_radial=0, angle_azimuthal=0, angle_radial_out=0)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_x=rotation_x,
                                           rotation_y=0, rotation_z=0)
    from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
    beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                            movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None
    from shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror import S4ParaboloidMirror
    optical_element = S4ParaboloidMirror(name='Parabolic Mirror 2 (pass 2)', boundary_shape=boundary_shape,
                                         at_infinity=0,
                                         surface_calculation=0,
                                         p_focus=0.020000, q_focus=0.020000, grazing_angle=np.pi/2,
                                         # for surface_calculation=0
                                         parabola_parameter=1.000000, pole_to_focus=1.000000,  # for external
                                         is_cylinder=0, cylinder_direction=0, convexity=1,
                                         f_reflec=0,  # reflectivity of surface: 0=no reflectivity, 1=full polarization
                                         f_refl=6,
                                         # for f_reflec=1: file: 0=prerefl, 2=(mrad, refl), 3=(eV, refl), 4=(eV, mrad, refl); 1=refr index, 5=xraylib, 6=dabax
                                         file_refl='<none>',  # for f_refl=0,2,3,4
                                         refraction_index=0.99999 + 0.001j,  # for f_refl=1
                                         coating_material='Ni', coating_density=8.902,  # for f_refl=5,6
                                         coating_roughness=0,  # for f_refl=0,1,5,6
                                         dabax=None,
                                         # if using dabax (f_reflec=1,f_refl=6), instance of DabaxXraylib() (use None for default)
                                         )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=0)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror import S4ParaboloidMirrorElement
    beamline_element = S4ParaboloidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                 movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name=' Screen at source position (pass 2)', boundary_shape=boundary_shape,
                               i_abs=0,  # attenuation: 0=No, 1=prerefl file, 2=xraylib, 3=dabax
                               i_stop=0,  # 0=slit or aperture, 1=beam stop
                               thick=0,  # for i_abs>0
                               file_abs='<specify file name>',  # for i_abs=1
                               material='Au', density=19.3,  # for i_abs=2,3
                               dabax=None,
                               # if using dabax (i_abs=3), instance of DabaxXraylib() (use None for default)
                               )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.02, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Ellipse
    boundary_shape = Ellipse(a_axis_min=-4.5e-06, a_axis_max=4.5e-06, b_axis_min=-4.5e-06, b_axis_max=4.5e-06)

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Pinhole', boundary_shape=boundary_shape,
                               i_abs=0,  # attenuation: 0=No, 1=prerefl file, 2=xraylib, 3=dabax
                               i_stop=0,  # 0=slit or aperture, 1=beam stop
                               thick=0,  # for i_abs>0
                               file_abs='<specify file name>',  # for i_abs=1
                               material='Au', density=19.3,  # for i_abs=2,3
                               dabax=None,
                               # if using dabax (i_abs=3), instance of DabaxXraylib() (use None for default)
                               )

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)


    return beam, beam_stopper


#
# main
#
from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show

do_calculate = 0

##

beam, beam_stopper = run_beamline(rotation_x=0)
I1 = beam_stopper.intensity(nolost=2)  # stopped
R = beam_stopper.get_column(20)
weight = np.ones_like(R) * 0.1263
weight[R < 125e-6] = 0.0353
weight[R > 1200e-6] = 0.0
c23 = beam.get_column(23)
cflag = beam.get_column(10)
I2 = beam.intensity(nolost=1)
I2w = (c23[cflag > 0] * weight[cflag > 0]).sum()
print("------------------angle=0, I1, I2, I2w: ", I1, I2, I2w, I1 + I2w)


if do_calculate:
    ##
    angles = np.linspace(-20e-3, 20e-3, 1001)

    I1_ARRAY = np.zeros_like(angles)
    I2_ARRAY = np.zeros_like(angles)
    I2w_ARRAY = np.zeros_like(angles)

    for i, angle in enumerate(angles):
        beam, beam_stopper = run_beamline(rotation_x=angle)
        I1 = beam_stopper.intensity(nolost=2) # stopped
        R = beam_stopper.get_column(20)
        weight = np.ones_like(R) * 0.1263
        weight[R < 125e-6] = 0.0353
        weight[R > 1200e-6] = 0.0
        c23 = beam.get_column(23)
        cflag = beam.get_column(10)
        I2 = beam.intensity(nolost=1)
        I2w =  (c23[cflag > 0] * weight[cflag > 0]).sum()
        print("------------------i, angle, I1, I2, I2w: ", i, angle, I1, I2, I2w )
        I1_ARRAY[i] = I1
        I2_ARRAY[i] = I2
        I2w_ARRAY[i] = I2w


    Itot = I1_ARRAY + I2w_ARRAY
    # Stack columns: each column will appear in the file
    data = np.column_stack((angles, I1_ARRAY, I2_ARRAY, I2w_ARRAY, Itot))

    np.savetxt(
        "fiber.dat",
        data,
        header="# angle_rad I1_single_ref I2_double_ref I2w_weighted I_total",
        fmt="%.6e"
    )
    print("File written to disk: fiber.dat")

else:
    angles, I1_ARRAY, I2_ARRAY, I2w_ARRAY, Itot = np.loadtxt("fiber.dat", unpack=True)


plot(
     1e3 * angles, Itot, xtitle="angle [mrad]", ytitle="intensity", show=0, grid=1,
     )


plot(1e3 * angles, I1_ARRAY,
     1e3 * angles, I2_ARRAY,
     1e3 * angles, I2w_ARRAY,
     1e3 * angles, Itot,
     legend=['I1 single ref','I2 double reflection in pinhole','I2w weighted','Itotal I1+I2w'], xtitle="angle [mrad]",
     grid=1, linestyle=[':',':',':',None]
     )