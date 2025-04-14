import numpy
import xraylib

def transfocator_compute_parameters(
            photon_energy_ev,
            nlenses_target=[1,1,1],
            symbol=["Be","Be","Be"],
            density=[1.845,1.845,1.845],
            nlenses_radii = [500e-6,1000e-6,1500e-6],
            sigmaz=6.00e-6,
            sigmazp=10.0e-6,
            tf_p=59.60,
            tf_q=38.00 ):

    deltas = [(1.0 - xraylib.Refractive_Index_Re(symbol[i], photon_energy_ev * 1e-3, density[i])) \
              for i in range(len(symbol))]

    focal_f = _transfocator_calculate_focal_distance(deltas=deltas,
                                                     nlenses=nlenses_target,
                                                     radii=nlenses_radii)
    focal_q = 1.0 / (1.0 / focal_f - 1.0 / tf_p)

    demagnification = focal_q / tf_p
    div_q = sigmazp / demagnification

    source_demagnified = sigmaz * demagnification # focal_q / tf_p
    s_target = numpy.sqrt( (div_q * (tf_q - focal_q))**2 + (source_demagnified)**2 )

    return (s_target, focal_f, focal_q, div_q)

def _transfocator_calculate_focal_distance(deltas=[0.999998], nlenses=[1], radii=[500e-6]):

    inverse_focal_distance = 0.0
    for i, nlensesi in enumerate(nlenses):
        if nlensesi > 0:
            focal_distance_i = radii[i] / (2. * nlensesi * deltas[i])
            inverse_focal_distance += 1.0 / focal_distance_i
    if inverse_focal_distance == 0:
        return 99999999999999999999999999.
    else:
        return 1.0 / inverse_focal_distance


def transfocator_guess_configuration(
            photon_energy_ev=14000.0,
            s_target=10e-6,
            symbol=["Be","Be","Be"],
            density=[1.845,1.845,1.845],
            nlenses_max=[15,3,1],
            nlenses_radii=[500e-6,1000e-6,1500e-6],
            sigmaz=6.46e-6,
            sigmazp=6.46e-6,
            tf_p=59.60, tf_q=38.00,
            verbose=1 ):

    if s_target < sigmaz * tf_q / tf_p:
        print("Source size FWHM is: %f um" % (1e6 * 2.355 * sigmaz))
        print("Maximum Demagnifications is: %f" % (tf_p / tf_q))
        print("Minimum possible size is: %f um" % (1e6 * 2.355 * sigmaz * tf_q / tf_p))
        print("Error: redefine size")
        return None


    combinations = []
    for i in range(nlenses_max[0]+1):
        for j in range(nlenses_max[1]+1):
            for k in range(nlenses_max[2]+1):
                combinations.append([i,j,k])

    N = len(combinations)
    sizes = numpy.zeros(N, dtype=float)
    differences = numpy.zeros(N, dtype=float)

    for i, nlenses_target in enumerate(combinations):
        (S, f, q_f, div) = transfocator_compute_parameters(
            photon_energy_ev=photon_energy_ev,
            nlenses_target=nlenses_target,
            symbol=symbol,
            density=density,
            nlenses_radii=nlenses_radii,
            sigmaz=sigmaz,
            sigmazp=sigmazp,
            tf_p=tf_p,
            tf_q=tf_q,
            )
        sizes[i] = S
        differences[i] = numpy.abs(S - s_target)

    i_optimum = numpy.argsort( numpy.abs(sizes - s_target) )

    if verbose:
        print("The best 3 configs are:")
        for i in range(3):
            print("    %s got: %.3f um, target: %.3f um" % (repr(combinations[i_optimum[i]]), 1e6 * sizes[i_optimum[i]], 1e6 * s_target))

    return combinations[i_optimum[0]]


def run_S4_ideal_lens(p_slit=27.3, tf_p=100.0, tf_q=100.0, focal=1.0):
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
    coordinates = ElementCoordinates(p=p_slit, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens import S4IdealLens
    optical_element = S4IdealLens(name='Ideal Lens', focal_x=focal, focal_y=focal)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=tf_p-p_slit, q=tf_q, angle_radial=0, angle_azimuthal=0,
                                     angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_lens import S4IdealLensElement
    beamline_element = S4IdealLensElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # test plot
    if 0:
        from srxraylib.plot.gol import plot_scatter
        plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
                     title='(Intensity,Photon Energy)', plot_histograms=0)
        plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')
    return beam, light_source

def run_S4_real_lens(photon_energy=14000.0, n_lens=[1, 2, 4, 8, 1, 2, 1]):
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
        photon_energy=photon_energy,  # Photon energy (in eV)
        delta_e=0.0,  # Photon energy width (in eV)
        ng_e=100,  # Photon energy scan number of points
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread=1,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
        harmonic_number=3,  # harmonic number
        flag_autoset_flux_central_cone=1,  # value to set the flux peak
        flux_central_cone=493975722611785.4,  # value to set the flux peak
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
                                     piling_thickness=[0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003],  # syned stuff
                                     boundary_shape=boundary_shape,
                                     # syned stuff, replaces "diameter" in the shadow3 append_lens
                                     material=['Be', 'Be', 'Be', 'Be', 'Be', 'Be', 'Be'],
                                     # the material for ri_calculation_mode > 1
                                     density=[1.848, 1.848, 1.848, 1.848, 1.848, 1.848, 1.848],
                                     # the density for ri_calculation_mode > 1
                                     thickness=[4.9999999999999996e-05, 4.9999999999999996e-05, 4.9999999999999996e-05,
                                                4.9999999999999996e-05, 4.9999999999999996e-05, 4.9999999999999996e-05,
                                                4.9999999999999996e-05],
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
                                     radius=[0.0005, 0.0005, 0.0005, 0.0005, 0.001, 0.001, 0.0015],
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
                                     empty_space_after_last_interface=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
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
    return beam, light_source

def get_sigmas_all(light_source):
    E0, delta_e, npoints = light_source.get_magnetic_structure().get_energy_box()
    u = light_source.get_magnetic_structure()

    s, sp = light_source.get_undulator_photon_beam_sizes(
        undulator_E0=E0,
        undulator_length=u.length(),
    )

    if u.get_flag_emittance():
        sigma_x, sigdi_x, sigma_z, sigdi_z = light_source.get_electron_beam().get_sigmas_all()
    else:
        sigma_x, sigdi_x, sigma_z, sigdi_z = 0, 0, 0, 0

    Qa, Qs = light_source._get_q_a_and_q_s()

    Sx = numpy.sqrt((s * Qs) ** 2 + sigma_x ** 2)
    Sz = numpy.sqrt((s * Qs) ** 2 + sigma_z ** 2)
    Spx = numpy.sqrt((sp * Qa) ** 2 + sigdi_x ** 2)
    Spz = numpy.sqrt((sp * Qa) ** 2 + sigdi_z ** 2)

    return Sx, Spx, Sz, Spz

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_show

    photon_energy = 14000.0
    n_lens = [1, 2, 4, 8, 1, 2, 1]

    beam, light_source = run_S4_real_lens(photon_energy=photon_energy, n_lens=n_lens)

    Sigma_x, Sigma_xp, Sigma_z, Sigma_zp = get_sigmas_all(light_source)

    tf_p = 62.052
    tf_q = 39.378

    photon_energy_kev = photon_energy * 1e-3
    qq = numpy.linspace(-40.0, 0, 100)
    nlenses_current = [n_lens[0] + n_lens[1] + n_lens[2] + n_lens[3],
                       n_lens[4] + n_lens[5],
                       n_lens[6]]

    (size, f, q_f, div) = transfocator_compute_parameters(
        photon_energy_kev * 1e3,
        nlenses_current,
        symbol=["Be","Be","Be"],
        density=[1.845,1.845,1.845],
        nlenses_radii = [500e-6,1000e-6,1500e-6],
        sigmaz=Sigma_z,
        sigmazp=Sigma_zp,
        tf_p=tf_p,
        tf_q=tf_q + qq,
        )

    (size_H, f_H, q_f_H, div_H) = transfocator_compute_parameters(
        photon_energy_kev * 1e3,
        nlenses_current,
        symbol=["Be","Be","Be"],
        density=[1.845,1.845,1.845],
        nlenses_radii = [500e-6,1000e-6,1500e-6],
        sigmaz=Sigma_x,
        sigmazp=Sigma_xp,
        tf_p=tf_p, tf_q=tf_q + qq,
        )

    #
    #  shadow
    #

    print(">>>>>>>>>>>>tf_p, tf_q, f, q_f: ", tf_p, tf_q, f, q_f)
    # beam, light_source = run_S4_ideal_lens(p_slit=27.3, tf_p=tf_p, tf_q=tf_q, focal=f)
    # beam, light_source = run_S4_real_lens(photon_energy=14000.0, n_lens=[1, 2, 4, 8, 1, 2, 1])
    print(">>>>>>>>>>>>sigmas_all: ", get_sigmas_all(light_source))
    from shadow4.tools.beamline_tools import focnew, focnew_scan

    ticket = focnew(beam=beam, nolost=1, mode=0, center=[0,0])
    y = qq
    ylist = [focnew_scan(ticket["AX"], y) * 1e6,
             focnew_scan(ticket["AZ"], y) * 1e6,
             focnew_scan(ticket["AT"], y) * 1e6]

    FWHM_X = numpy.zeros_like(y)
    FWHM_Z = numpy.zeros_like(y)
    FWHM_FROM_STDEV_X = numpy.zeros_like(y)
    FWHM_FROM_STDEV_Z = numpy.zeros_like(y)
    for i in range(y.size):
        bb = beam.duplicate()
        bb.retrace(y[i])
        tx = bb.histo1(1, nolost=1, ref=23)
        tz = bb.histo1(3, nolost=1, ref=23)

        FWHM_X[i] = tx['fwhm']
        FWHM_Z[i] = tz['fwhm']

        FWHM_FROM_STDEV_X[i] = 2.355 * 1e6 * bb.get_standard_deviation(1, nolost=1, ref=23)  # tx['fwhm'] * 1e6
        FWHM_FROM_STDEV_Z[i] = 2.355 * 1e6 * bb.get_standard_deviation(3, nolost=1, ref=23)  # tz['fwhm'] * 1e6

    plot(
        qq, 2.355 * size_H * 1e6,
        qq, 2.355 * size * 1e6,
        y,  FWHM_X * 1e6,
        y,  FWHM_Z * 1e6,
        # y, FWHM_FROM_STDEV_X,
        # y, FWHM_FROM_STDEV_Z,
        title="FWHM for config: %s (%s) " % (nlenses_current, n_lens),
        xtitle="Y [m]",
        ytitle="size [$\mu$m]",
        legend=["X analytical","Z analytical", \
            "X ray tracing (FWHM histogram)", "Z ray tracing (FWHM histogram)", \
            "X ray tracing (2.35 * stdev)","Z ray tracing (2.35 * stdev)"],
        linestyle=[None,None,':',':','--','--'],
        color=['b', 'r', 'b', 'r', 'b', 'r',],
        grid=1, show=0)




    print("Analytical size at sample H [um]: ", 2.355 * size_H[-1] * 1e6)
    print("Analytical size at sample V [um]: ", 2.355 * size[-1]   * 1e6)
    print("Ray-tracing size at sample H [um]: ", FWHM_X[-1] * 1e6)
    print("Ray-tracing size at sample V [um]: ", FWHM_Z[-1] * 1e6)
    print("Ratio Raytracing/Analytical H [um]: ", FWHM_X[-1] / (2.355 * size_H[-1]))
    print("Ratio Raytracing/Analytical V [um]: ", FWHM_Z[-1] / (2.355 * size[-1]  ))

    #
    # Guess config
    #
    if 1:
        s_target_V = size[-1]
        s_target_H = size_H[-1]

        nlenses_optimum_V = transfocator_guess_configuration(photon_energy_ev=photon_energy_kev * 1e3,
                                                             s_target=s_target_V,
                                                             symbol=["Be","Be","Be"],
                                                             density=[1.845,1.845,1.845],
                                                             nlenses_max = [15,3,1],
                                                             nlenses_radii = [500e-6,1000e-6,1500e-6],
                                                             sigmaz=Sigma_z,
                                                             sigmazp=Sigma_zp,
                                                             tf_p=tf_p,
                                                             tf_q=tf_q,
                                                             verbose=1 )

        nlenses_optimum_H = transfocator_guess_configuration(photon_energy_ev=photon_energy_kev * 1e3,
                                                             s_target=s_target_H,
                                                             symbol=["Be","Be","Be"],
                                                             density=[1.845,1.845,1.845],
                                                             nlenses_max = [15,3,1],
                                                             nlenses_radii = [500e-6,1000e-6,1500e-6],
                                                             sigmaz=Sigma_x,
                                                             sigmazp=Sigma_xp,
                                                             tf_p=tf_p,
                                                             tf_q=tf_q,
                                                             verbose=1 )

        print("Optimum lens configuration in V is: ", nlenses_optimum_V)
        print("Optimum lens configuration in H is: ", nlenses_optimum_H)

    plot_show()