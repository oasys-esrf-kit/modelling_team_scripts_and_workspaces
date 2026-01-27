import numpy


def run_beamline(pitch1=0,pitch2=0,roll1=0,roll2=0,t2=0):
    from shadow4.beamline.s4_beamline import S4Beamline
    import numpy as np

    beamline = S4Beamline()

    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam
    electron_beam = S4ElectronBeam(energy_in_GeV=6,energy_spread=0.001,current=0.2)
    electron_beam.set_sigmas_all(sigma_x=2.87429e-05, sigma_y=5.15531e-06, sigma_xp=4.18025e-06, sigma_yp=1.93975e-06)
    electron_beam.set_dispersion_all(0,0,0, 0)

    # magnetic structure
    from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
    source = S4UndulatorGaussian(
        period_length     = 0.0186,     # syned Undulator parameter (length in m)
        number_of_periods = 134.40860215053763, # syned Undulator parameter
        photon_energy     = 17000.0, # Photon energy (in eV)
        delta_e           = 20.0, # Photon energy width (in eV)
        ng_e              = 100, # Photon energy scan number of points
        flag_emittance    = 1, # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread = 1, # when sampling rays: Use e- energy spread (0=No, 1=Yes)
        harmonic_number    = 1, # harmonic number
        flag_autoset_flux_central_cone  = 1, # value to set the flux peak
        flux_central_cone  = 555727766916010.9, # value to set the flux peak
        )


    # light source
    from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
    light_source = S4UndulatorGaussianLightSource(name='Undulator Gaussian 17 keV', electron_beam=electron_beam, magnetic_structure=source,nrays=50000,seed=12345)
    beam = light_source.get_beam()

    beamline.set_light_source(light_source)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.001, x_right=0.001, y_bottom=-0.0005, y_top=0.0005)

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='FE V and H slit', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=23.513, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Primary Slits', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=3.778, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Attenuator (3Axis)', boundary_shape=boundary_shape,
        i_abs=3, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.758, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number DMM1
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.15, y_top=0.15)

    from shadow4.beamline.optical_elements.multilayers.s4_plane_multilayer import S4PlaneMultilayer
    optical_element = S4PlaneMultilayer(name='Generic Multilayer',boundary_shape=boundary_shape,
        f_refl=0,file_refl='mlayerPdCx15.dat', structure='[C,Pt]x30+Si', period=50.000000, Gamma=0.400000)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=1.492, q=0, angle_radial=1.566083938, angle_azimuthal=4.71238898, angle_radial_out=1.566083938)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    # srio
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0.02, offset_z=0,
                                           rotation_x=pitch1, rotation_y=roll1, rotation_z=0)

    from shadow4.beamline.optical_elements.multilayers.s4_plane_multilayer import S4PlaneMultilayerElement
    beamline_element = S4PlaneMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number DMM2
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.3643, y_top=0.3643)

    from shadow4.beamline.optical_elements.multilayers.s4_sphere_multilayer import S4SphereMultilayer
    optical_element = S4SphereMultilayer(name='Generic Multilayer',boundary_shape=boundary_shape,
        surface_calculation=0,is_cylinder=1,cylinder_direction=0,
        convexity=1,radius=1.000000,p_focus=29.541000,q_focus=24.859000,
        grazing_angle=0.004712,
        f_refl=0,file_refl='mlayerPdCx15.dat', structure='[C,Pt]x30+Si', period=50.000000, Gamma=0.400000)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.566083938, angle_azimuthal=3.141592654, angle_radial_out=1.566083938)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    #srio2
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=t2, offset_z=0,
                                           rotation_x=pitch2, rotation_y=roll2, rotation_z=0)

    from shadow4.beamline.optical_elements.multilayers.s4_sphere_multilayer import S4SphereMultilayerElement
    beamline_element = S4SphereMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=4.71238898, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='White beam Stop', boundary_shape=boundary_shape,
        i_abs=3, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=1.231, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='TF Screen (before)', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.724, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Cleaning slit', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.827, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=1.48, q=0, angle_radial=0, angle_azimuthal=1.570796327, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=4.71238898, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='First Safety Aperture', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.957, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Attenuator (5Axis)', boundary_shape=boundary_shape,
        i_abs=3, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=17.832, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='HSS H Secondary Source', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=1.768, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Secondary Safety Aperture', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.502, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='KB entrance screen', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=144.998, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.03, y_top=0.03)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Generic Mirror', boundary_shape=boundary_shape,
        surface_calculation=0,
        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
        p_focus=199.900000, q_focus=0.100000, grazing_angle=0.036400,
        is_cylinder=1, cylinder_direction=0, convexity=1,
        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999+0.001j,
        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0.0335, angle_radial=1.534396327, angle_azimuthal=0, angle_radial_out=1.534396327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.01, x_right=0.01, y_bottom=-0.012537, y_top=0.012537)

    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Generic Mirror', boundary_shape=boundary_shape,
        surface_calculation=0,
        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
        p_focus=145.550000, q_focus=0.050000, grazing_angle=0.036400,
        is_cylinder=1, cylinder_direction=0, convexity=1,
        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999+0.001j,
        coating_material='Si', coating_density=2.33, coating_roughness=0)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.0165, q=0, angle_radial=1.534396327, angle_azimuthal=1.570796327, angle_radial_out=1.534396327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.011, x_right=0.011, y_bottom=-0.0009125, y_top=0.0009125)

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='CORRECTION !!!!  Generic Beam Screen/Slit/Stopper/Attenuator HSS', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
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
    coordinates = ElementCoordinates(p=0, q=0.05, angle_radial=0, angle_azimuthal=4.71238898, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)
    return beam, footprint

#
# main
#
from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show

# ANGLE = numpy.linspace(-4e-6, 4e-6, 50)
# title = "pitch M2"
#
ANGLE = numpy.linspace(-500e-6, 500e-6, 15)
title = "roll M2"

# ANGLE = numpy.linspace(-0.002, 0.002, 50) # it is not angle, it is translation
# title = "translation M2"

# reference
beam, footprint = run_beamline()
ticket = beam.histo2(1, 3, nbins_h=300, nbins_v=300, xrange=[-100e-9, 100e-9], yrange=[-100e-9, 100e-9], nolost=1, ref=23)

INTENSITY_REF = beam.intensity(nolost=1)
if ticket['fwhm_h'] is not None: FWHM_H_REF = ticket['fwhm_h']
if ticket['fwhm_v'] is not None: FWHM_V_REF = ticket['fwhm_v']

# print(ticket['histogram'].shape, ticket['bin_h_center'].shape, ticket['bin_v_center'].shape)
# plot_image_with_histograms(ticket['histogram'], 1e9 * ticket['bin_h_center'], 1e9 * ticket['bin_v_center'],
#     title=title, xtitle="column 1", ytitle="column 3",
#     xrange=[-100, 100], yrange=[-100, 100],
#     cmap='jet', add_colorbar=0, figsize=(8, 8), histo_path_flag=1, show=1)

#
# LOOP
#

INTENSITY = numpy.zeros_like(ANGLE)
FWHM_H = numpy.zeros_like(ANGLE)
FWHM_V = numpy.zeros_like(ANGLE)

for i, angle in enumerate(ANGLE):
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> i: ", i, " of", ANGLE.size)
    beam, footprint = run_beamline(t2=angle) # roll2=angle)
    ticket = beam.histo2(1, 3, nbins_h=300, nbins_v=300, xrange=[-100e-9, 100e-9], yrange=[-100e-9, 100e-9], nolost=1,
                         ref=23)
    INTENSITY[i] = beam.intensity(nolost=1)
    if ticket['fwhm_h'] is not None: FWHM_H[i] = ticket['fwhm_h']
    if ticket['fwhm_v'] is not None: FWHM_V[i] = ticket['fwhm_v']

from srxraylib.plot.gol import plot, plot_show

plot(1e6 * ANGLE, INTENSITY / INTENSITY_REF,
     title=title,
     xtitle="angle [urad]", ytitle="Intensity / Intensity_reference",
     figsize=(10, 6), grid=1, show=0)

plot(1e6 * ANGLE, (FWHM_H / FWHM_H_REF),
     1e6 * ANGLE, (FWHM_V / FWHM_V_REF),
     # 1e6 * ANGLE, 1e9 * INTENSITY / INTENSITY.max(),
     title=title,
     xtitle="angle [urad]", ytitle="FWHM / FWHM_reference", legend=["H", "V"],
     linestyle=[None, None,],
     figsize=(10, 6), grid=1, show=1)

# plot(1e3 * ANGLE, INTENSITY / INTENSITY_REF,
#      title=title,
#      xtitle="translation [mm]", ytitle="Intensity / Intensity_reference",
#      figsize=(10, 6), grid=1, show=0)
#
# plot(1e3 * ANGLE, (FWHM_H / FWHM_H_REF),
#      1e3 * ANGLE, (FWHM_V / FWHM_V_REF),
#      # 1e6 * ANGLE, 1e9 * INTENSITY / INTENSITY.max(),
#      title=title,
#      xtitle="translation [mm]", ytitle="FWHM / FWHM_reference", legend=["H", "V"],
#      linestyle=[None, None,],
#      figsize=(10, 6), grid=1, show=1)