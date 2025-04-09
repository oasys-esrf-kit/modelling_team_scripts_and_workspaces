import numpy
# import Shadow
import xraylib

def transfocator_calculate_focal_distance(deltas=[0.999998], nlenses=[1], radii=[500e-6]):
    inverse_focal_distance = 0.0
    for i,nlensesi in enumerate(nlenses):
        if nlensesi > 0:
            focal_distance_i = radii[i] / (2. * nlensesi * deltas[i])
            inverse_focal_distance += 1.0 / focal_distance_i
    return 1.0 / inverse_focal_distance


def run_beamline(
        n_lens=[1, 2, 4, 8, 1, 2, 1],
        piling_thickness=[0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003],  # syned stuff
        material=['Be', 'Be', 'Be', 'Be', 'Be', 'Be', 'Be'],
        density=[1.848, 1.848, 1.848, 1.848, 1.848, 1.848, 1.848],
        thickness=[4.9999999999999996e-05, 4.9999999999999996e-05, 4.9999999999999996e-05,
                   4.9999999999999996e-05, 4.9999999999999996e-05, 4.9999999999999996e-05,
                   4.9999999999999996e-05],
        radius=[0.0005, 0.0005, 0.0005, 0.0005, 0.001, 0.001, 0.0015],
        empty_space_after_last_interface=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ):
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam
    electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
    electron_beam.set_sigmas_all(sigma_x=2.80624e-05, sigma_y=5.14918e-06, sigma_xp=5.01922e-06, sigma_yp=1.94206e-06)

    # magnetic structure
    from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
    source = S4UndulatorGaussian(
        period_length=0.035,  # syned Undulator parameter (length in m)
        number_of_periods=65.714,  # syned Undulator parameter
        photon_energy=14000.0,  # Photon energy (in eV)
        delta_e=0.0,  # Photon energy width (in eV)
        ng_e=100,  # Photon energy scan number of points
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread=1,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
        harmonic_number=3,  # harmonic number
        flag_autoset_flux_central_cone=0,  # value to set the flux peak
        flux_central_cone=10000000000.0,  # value to set the flux peak
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
    light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam,
                                                  magnetic_structure=source, nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    beamline.set_light_source(light_source)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.00025, x_right=0.00025, y_bottom=-0.00025, y_top=0.00025)

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
    from shadow4.beamline.optical_elements.refractors.s4_transfocator import S4Transfocator

    optical_element = S4Transfocator(name='TF',
                                     n_lens=n_lens,
                                     piling_thickness=piling_thickness,  # syned stuff
                                     boundary_shape=boundary_shape,
                                     # syned stuff, replaces "diameter" in the shadow3 append_lens
                                     material=material,
                                     # the material for ri_calculation_mode > 1
                                     density=density,
                                     # the density for ri_calculation_mode > 1
                                     thickness=thickness,
                                     radius=radius,
                                     empty_space_after_last_interface=empty_space_after_last_interface,
                                     # syned stuff, lens thickness [m] (distance between the two interfaces at the center of the lenses)
                                     surface_shape=[2, 2, 2, 2, 2, 2, 2],
                                     # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                     # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
                                     convex_to_the_beam=[0, 0, 0, 0, 0, 0, 0],
                                     # for surface_shape: convexity of the first interface exposed to the beam 0=No, 1=Yes
                                     cylinder_angle=[0, 0, 0, 0, 0, 0, 0],
                                     # for surface_shape: 0=not cylindricaL, 1=meridional 2=sagittal
                                     ri_calculation_mode=[2, 2, 2, 2, 2, 2, 2],
                                     # source of refraction indices and absorption coefficients
                                     # 0=User, 1=prerefl file, 2=xraylib, 3=dabax
                                     prerefl_file=['NONE SPECIFIED', 'NONE SPECIFIED', 'NONE SPECIFIED',
                                                   'NONE SPECIFIED', 'NONE SPECIFIED', 'NONE SPECIFIED',
                                                   'NONE SPECIFIED'],
                                     # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                                     refraction_index=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                                     # for ri_calculation_mode=1: n (real)
                                     attenuation_coefficient=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                     # for ri_calculation_mode=1: mu in cm^-1 (real)
                                     dabax=None,  # the pointer to dabax library
                                     # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                                     conic_coefficients1=[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                                     # for surface_shape = 3: the conic coefficients of the single lens interface 1
                                     conic_coefficients2=[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                                     # for surface_shape = 3: the conic coefficients of the single lens interface 2
                                     )

    import numpy
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=34.752, q=39.378, angle_radial=0, angle_azimuthal=0,
                                     angle_radial_out=3.141592654)
    movements = None
    from shadow4.beamline.optical_elements.refractors.s4_transfocator import S4TransfocatorElement
    beamline_element = S4TransfocatorElement(optical_element=optical_element, coordinates=coordinates,
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

if __name__ == "__main__":

    print("setting Transfocator for ID30B")


    #
    # transfocator id30B
    #

    # geometry of the TF:  ALL VARIABLES ARE ARRAYS!!!!
    n_axes = 7
    n_lens_in_axis = [1,  2,  4,  8,   1,   2,   1]  # number of lenses per axis
    radius = [500e-6, 500e-6, 500e-6, 500e-6, 1000e-6, 1000e-6, 1500e-6]
    piling_thickness = [0.003] * n_axes            #total piling thickness of each single lens in m
    material = ['Be'] * n_axes
    thickness = [50e-06] * n_axes
    empty_space_after_last_interface = [0.0] * n_axes

    tf_on_off  = [1, 1, 1, 1, 1, 1, 1]                      # set (1) or unset (0)

    # position of the TF measured from the center of the transfocator
    slit_position = 27.3
    tf_p = 62.052 # m
    tf_q = 39.378 # m

    print("\n=============== TRANSFOCATOR INPUTS ================================")
    print("                               n_axes: ", n_axes)
    print("                       n_lens_in_axis: ", n_lens_in_axis)
    print("                            tf_on_off: ", tf_on_off)
    print("                            radii [m]: ", radius)
    print("            Lens piling thickness [m]: ", piling_thickness)
    print("                       Lens materials: ", material)
    print("tf_p: distance from previous oe to center of TF: %f cm"%(tf_p))
    print("    tf_q: distance from center of TF to next oe: %f cm"%(tf_q))
    print("====================================================================")

    n_lens = numpy.array(n_lens_in_axis) * numpy.array(tf_on_off)
    empty_slots = numpy.array(n_lens_in_axis) - n_lens
    empty_space_after_last_interface = empty_slots * numpy.array(piling_thickness)



    #
    # ideal calculations
    #

    photon_energy_ev = 14000.0

    density = []
    for i in range(n_axes):
        xrl_density = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(material[i]))
        density.append(xrl_density)

    deltas = []
    for i in range(n_axes):
        delta = 1.0 - xraylib.Refractive_Index_Re(material[i], photon_energy_ev * 1e-3, density[i])
        deltas.append(delta)

    # print("                       Lens densities: ", density)
    # print("                       Lens delta: ", deltas)


    #
    # analytical calculations
    #

    focal_distance = transfocator_calculate_focal_distance(nlenses=n_lens, deltas=deltas, radii=radius)
    focal_q = 1.0/( (1.0 / focal_distance) - (1.0 / tf_p) )
    demagnification_factor = tf_p / focal_q

    print("\n=============== ANALYTICAL CALCULATIONS ================================")
    print("Photon energy is: %f eV" % (photon_energy_ev))
    print("Focal distance is: %f m" % (focal_distance))
    print("for p: %f m we get focal q at: %f m (%f m from q position that was %f m)" % (tf_p, focal_q, focal_q - tf_q, tf_q))
    print("Demagnification factor: %f" % (demagnification_factor))

    if 1:
        #
        # shadow calculations
        #
        print("\n=============== SHADOW INPUTS ================================")
        print(">>>>>>>         tf_p tf_q: ", tf_p, tf_q)
        print(">>>>>>>         tf_p - slit_position", tf_p - slit_position)
        print(">>>>>>>         focal to focus p q: ", 1.0 / (1 / tf_p + 1 / tf_q))


        print(">>>>> n_lens: ", n_lens)
        print(">>>>> piling_thickness: ", piling_thickness)
        print(">>>>> material: ", material)
        print(">>>>> density: ", density)
        print(">>>>> thickness: ", thickness)
        print(">>>>> radius: ", radius)
        print(">>>>> empty_space_after_last_interface: ", empty_space_after_last_interface)
        # print(">>>>>>> deltas: ", deltas)


        from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show

        # WARNING: NO incremental result allowed!!"
        beam = run_beamline()

        ticket = beam.histo2(1, 3, nbins_h=100, nbins_v=100, xrange=[-0.0019283344812185753, 0.0019074550298156779],
                             yrange=[-0.0019283344812185753, 0.0019074550298156779], nolost=1, ref=23)

        title = "I: %.1f " % ticket['intensity']
        if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
        if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

        plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                                   title=title, xtitle="column 1", ytitle="column 3",
                                   cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=1)