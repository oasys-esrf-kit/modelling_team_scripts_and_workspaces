import numpy
import xraylib

"""
transfocator_id30b : transfocator for id13b:
        It can:

            1) guess the lens configuration (number of lenses for each type) for a given photon energy
            and target image size. Use transfocator_compute_configuration() for this task

            2) for a given transfocator configuration, compute the main optical parameters
                (image size, focal distance, focal position and divergence).
                Use transfocator_compute_parameters() for this task

            3) Performs full ray tracing. Use id30b_ray_tracing() for this task

        Note that for the optimization and parameters calculations the transfocator configuration is
        given in keywords. For ray tracing calculations many parameters of the transfocator are hard coded
        with the values of id30b

        See main program for examples.

        Dependencies:
            Numpy
            xraylib (to compute refracion indices)
            Shadow (for ray tracing only)
            matplotlib (for some plots of ray=tracing)

        Side effects:
            When running ray tracing some files are created.

        MODIFICATION HISTORY:
           2015-03-25 srio@esrf.eu, written

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright__ = "ESRF, 2015"

def transfocator_compute_parameters(
            photon_energy_ev,
            nlenses_target,
            symbol=["Be","Be","Be"],
            density=[1.845,1.845,1.845],
            nlenses_radii = [500e-4,1000e-4,1500e-4],
            lens_diameter=0.05,
            sigmaz=6.00e-4,
            sigmazp=10.0,
            alpha = 1.0,
            tf_p=5960, tf_q=3800 ):
    """
    Computes the parameters of the optical performances of a given transfocator configuration.

    returns a l

    All length units are cm

    :param photon_energy_ev:
    :param nlenses_target: a list with the lens configuration, i.e. the number of lenses of each type.
    :param symbol:         the chemical symbol of the lens material of each type. Default symbol=["Be","Be","Be"]
    :param density:        the density of each type of lens. Default: density=[1.845,1.845,1.845]
    :param nlenses_radii:  the radii in cm of each type of lens. Default: nlenses_radii = [500e-4,1000e-4,1500e-4]
    :param lens_diameter:    the physical diameter (acceptance) in cm of the lenses. If different for each type of lens,
                            consider the smaller one. Default: lens_diameter=0.05
    :param sigmaz:         the sigma (standard deviation) of the source in cm
    :param alpha:          an adjustable parameter in [0,1](see doc). Default: 0.55 (it is 0.76 for pure Gaussian beams)
    :param tf_p:           the distance source-transfocator in cm
    :param tf_q:           the distance transfocator-image in cm
    :return:               a list with parameters (image_size (sigma), lens_focal_distance,
                           focal_position from transfocator center, divergence of beam after the transfocator)
    """

    deltas = [(1.0 - xraylib.Refractive_Index_Re(symbol[i],photon_energy_ev*1e-3,density[i])) \
              for i in range(len(symbol))]
    focal_f = _transfocator_calculate_focal_distance( deltas=deltas,\
                                                       nlenses=nlenses_target,radii=nlenses_radii)
    focal_q = 1.0 / (1.0 / focal_f - 1.0 / tf_p)

    demagnification = focal_q / tf_p
    if 0:
        div_q = alpha  * lens_diameter / focal_q
    else:
        # new
        div_p = sigmazp
        div_q = sigmazp / demagnification

    print(">>>>>>>>>>>>>>>>>> D, slit, div_p, div_q: ", demagnification, 0.5e-3 / 27.3 / 2.355, div_p, div_q)
    #corrections
    source_demagnified = sigmaz * demagnification # focal_q / tf_p
    # if source_demagnified > lens_diameter: source_demagnified = lens_diameter
    s_target = numpy.sqrt( (div_q * (tf_q - focal_q))**2 + (source_demagnified)**2 )
    return (s_target * alpha, focal_f, focal_q, div_q)

def _transfocator_calculate_focal_distance(deltas=[0.999998], nlenses=[1], radii=[500e-4]):

    inverse_focal_distance = 0.0
    for i,nlensesi in enumerate(nlenses):
        if nlensesi > 0:
            focal_distance_i = radii[i] / (2. * nlensesi*deltas[i])
            inverse_focal_distance += 1.0 / focal_distance_i
    if inverse_focal_distance == 0:
        return 99999999999999999999999999.
    else:
        return 1.0 / inverse_focal_distance



def transfocator_compute_configuration(
            photon_energy_ev,
            s_target,
            symbol=["Be","Be","Be"],
            density=[1.845,1.845,1.845],
            nlenses_max=[15,3,1],
            nlenses_radii=[500e-4,1000e-4,1500e-4],
            lens_diameter=0.05,
            sigmaz=6.46e-4,
            sigmazp=6.46e-6,
            alpha = 0.55,
            tf_p=5960, tf_q=3800,
            verbose=1 ):
    """
    Computes the optimum transfocator configuration for a given photon energy and target image size.

    All length units are cm

    :param photon_energy_ev: the photon energy in eV
    :param s_target:       the target image size in cm.
    :param symbol:         the chemical symbol of the lens material of each type. Default symbol=["Be","Be","Be"]
    :param density:        the density of each type of lens. Default: density=[1.845,1.845,1.845]
    :param nlenses_max:    the maximum allowed number of lenases for each type of lens. nlenses_max = [15,3,1]
    :param nlenses_radii:  the radii in cm of each type of lens. Default: nlenses_radii = [500e-4,1000e-4,1500e-4]
    :param lens_diameter:    the physical diameter (acceptance) in cm of the lenses. If different for each type of lens,
                            consider the smaller one. Default: lens_diameter=0.05
    :param sigmaz:         the sigma (standard deviation) of the source in cm
    :param alpha:          an adjustable parameter in [0,1](see doc). Default: 0.55 (it is 0.76 for pure Gaussian beams)
    :param tf_p:           the distance source-transfocator in cm
    :param tf_q:           the distance transfocator-image in cm
    :param:verbose:        set to 1 for verbose text output
    :return:               a list with the number of lenses of each type.

    """
    if s_target < 2.35 * sigmaz * tf_q / tf_p:
        print("Source size FWHM is: %f um" % (1e4 * 2.35 * sigmaz))
        print("Maximum Demagnifications is: %f" % (tf_p / tf_q))
        print("Minimum possible size is: %f um" % (1e4 * 2.35 * sigmaz * tf_q / tf_p))
        print("Error: redefine size")
        return None

    deltas = [(1.0 - xraylib.Refractive_Index_Re(symbol[i],photon_energy_ev*1e-3,density[i])) \
              for i in range(len(symbol))]

    focal_q_target = _tansfocator_guess_focal_position( s_target,
                                                        p=tf_p, q=tf_q,
                                                        sigmaz=sigmaz,
                                                        sigmazp=sigmazp,
                                                        alpha=alpha,
                                                        lens_diameter=lens_diameter,
                                                        method=2)

    focal_f_target = 1.0 / (1.0 / focal_q_target + 1.0 / tf_p)
    div_q_target = alpha  * lens_diameter / focal_q_target

    #corrections for extreme cases
    source_demagnified = 2.35 * sigmaz * focal_q_target / tf_p
    if source_demagnified > lens_diameter: source_demagnified = lens_diameter

    s_target_calc = numpy.sqrt( (div_q_target * (tf_q - focal_q_target))**2 + source_demagnified**2)

    nlenses_target = _transfocator_guess_configuration(focal_f_target,
                                                       deltas=deltas,
                                                       nlenses_max=nlenses_max,
                                                       radii=nlenses_radii, )
    if verbose:
        print("transfocator_compute_configuration: focal_f_target: %f"%(focal_f_target))
        print("transfocator_compute_configuration: focal_q_target: %f cm"%(focal_q_target))
        print("transfocator_compute_configuration: s_target: %f um"%(s_target_calc*1e4))
        print("transfocator_compute_configuration: nlenses_target: ",nlenses_target)

    return nlenses_target


def transfocator_nlenses_to_slots(nlenses,nlenses_max=None):
    """
    converts the transfocator configuration from a list of the number of lenses of each type,
    into a list of active (1) or inactive (0) actuators for the slots.

    :param nlenses: the list with number of lenses (e.g., [5,2,0]
    :param nlenses_max: the maximum number of lenses of each type, usually powers of two minus one.
                        E.g. [15,3,1]
    :return: a list of on (1) and off (0) slots, e.g., [1, 0, 1, 0, 0, 1, 0]
            (first type: 1*1+0*2+1*4+0*8=5, second type: 0*1+1*2=2, third type: 0*1=0)
    """
    if nlenses_max == None:
        nlenses_max = nlenses

    ss = []
    for i,iopt in enumerate(nlenses):
        if iopt > nlenses_max[i]:
            print("Error: i:%d, nlenses: %d, nlenses_max: %d"%(i,iopt,nlenses_max[i]))
        ncharacters = len("{0:b}".format(nlenses_max[i]))

        si = list( ("{0:0%db}"%(ncharacters)).format(int(iopt)) )
        si.reverse()
        ss += si

    on_off = [int(i) for i in ss]
    #print("transfocator_nlenses_to_slots: nlenses_max: ",nlenses_max," nlenses: ",nlenses," slots: ",on_off)
    return on_off

def _tansfocator_guess_focal_position( s_target, p=5960., q=3800.0, sigmaz=6.46e-4, sigmazp=6.46e-6,
                                      alpha=0.66, lens_diameter=0.05, method=2):
    x = 1e15
    if method == 1: # simple sum
        AA = 2.35*sigmaz/p
        BB = -(s_target + alpha * lens_diameter)
        CC = alpha*lens_diameter*q

        cc = numpy.roots([AA,BB,CC])
        x = cc[1]
        return x

    if method == 2: # sum in quadrature
        AA = ( (2.35 * sigmaz)**2) / (p**2)
        BB = 0.0
        CC = alpha**2 * lens_diameter**2 - s_target**2
        DD = - 2.0 * alpha**2 * lens_diameter**2 * q
        EE = alpha**2 * lens_diameter**2 * q**2
        cc = numpy.roots([AA,BB,CC,DD,EE])
        for i,cci in enumerate(cc):
            if numpy.imag(cci) == 0:
                return numpy.real(cci)

    return x

def _transfocator_guess_configuration(focal_f_target,deltas=[0.999998],nlenses_max=[15],radii=[500e-4]):

    nn = len(nlenses_max)

    ncombinations = (1+nlenses_max[0]) * (1+nlenses_max[1]) * (1+nlenses_max[2])

    icombinations = 0
    aa = numpy.zeros((3,ncombinations),dtype=int)
    bb = numpy.zeros(ncombinations)
    for i0 in range(1+nlenses_max[0]):
        for i1 in range(1+nlenses_max[1]):
            for i2 in range(1+nlenses_max[2]):
                aa[0,icombinations] = i0
                aa[1,icombinations] = i1
                aa[2,icombinations] = i2
                bb[icombinations] = focal_f_target - _transfocator_calculate_focal_distance(deltas=deltas,nlenses=[i0,i1,i2],radii=radii)
                icombinations += 1
    bb1 = numpy.abs(bb)
    ibest = bb1.argmin()

    return (aa[:,ibest]).tolist()

#
#########################
#
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

def run_S4_real_lens(p_slit=27.3, tf_p=100.0, tf_q=100.0, focal=1.0):
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
                                     n_lens=[1, 2, 4, 8, 1, 2, 1],
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
    from srxraylib.plot.gol import plot

    Sigma_x  =  28.2381 # um
    Sigma_z  =  6.03339 # um
    Sigma_xp =  7.20068 # urad
    Sigma_zp =  5.51623 # urad

    alpha = 0.55

    tf_p = 6205.2
    tf_q = 3937.8

    photon_energy_kev = 14.0  # float(input("Enter photon energy in keV: "))

    #
    # Guess config
    #
    if 0:

        s_target_um       = 100.0   # float(input("Enter target focal dimension in microns: "))

        # slit_position = 27.3
        # tf_p = 62.052 # m
        # tf_q = 39.378 # m


        nlenses_optimum = transfocator_compute_configuration(photon_energy_kev * 1e3,
                                                             s_target_um * 1e-4,
                                                             symbol=["Be","Be","Be"],
                                                             density=[1.845,1.845,1.845],
                                                             nlenses_max = [15,3,1],
                                                             nlenses_radii = [500e-4,1000e-4,1500e-4],
                                                             lens_diameter=0.05,
                                                             sigmaz=Sigma_z * 1e-4,
                                                             alpha = alpha,
                                                             tf_p=tf_p,
                                                             tf_q=tf_q,
                                                             verbose=1 )

        print("Optimum lens configuration is: ",nlenses_optimum)

        if nlenses_optimum == None:
            raise Exception("Cannot find solution.")
        print("Guessed configuration: ", nlenses_optimum)


        #
        # compute size for given solution
        #

        # print("Activate slots: ",transfocator_nlenses_to_slots(nlenses_optimum,nlenses_max=[15,3,1]))
        # this calculates the parameters (image size, etc) for a given lens configuration


        for nlenses_current in [nlenses_optimum, [15,3,1]]:
            (size, f, q_f, div) = transfocator_compute_parameters(
                photon_energy_kev*1e3, nlenses_current,\
                symbol=["Be","Be","Be"],
                density=[1.845,1.845,1.845],
                nlenses_radii = [500e-4,1000e-4,1500e-4],
                lens_diameter=0.05,
                sigmaz=sigmaz, alpha = alpha,
                tf_p=tf_p, tf_q=tf_q,
                )

            # print("For configuration ",nlenses_current," we get: ")
            # print("  size: %f cm, focal length: %f cm, focal distance: %f cm, divergence: %f rad: "%(size, f, q_f, div))
            print("Config:  %s size: %f um, focal length: %f m, focal distance: %f m, divergence: %f u rad: "% \
                  (repr(nlenses_current), size * 1e4, f * 1e-2, q_f * 1e-2, div * 1e6))


    for nlenses_current in [[15,3,1]]:
        alpha = 1
        qq = numpy.linspace(-4000.0, 0, 100)
        (size, f, q_f, div) = transfocator_compute_parameters(
            photon_energy_kev*1e3, nlenses_current,\
            symbol=["Be","Be","Be"],
            density=[1.845,1.845,1.845],
            nlenses_radii = [500e-4,1000e-4,1500e-4],
            lens_diameter=0.05,
            sigmaz=Sigma_z * 1e-4,
            alpha = alpha,
            sigmazp=Sigma_zp * 1e-6,
            tf_p=tf_p, tf_q=tf_q + qq,
            )

        (size_H, f_H, q_f_H, div_H) = transfocator_compute_parameters(
            photon_energy_kev*1e3, nlenses_current,\
            symbol=["Be","Be","Be"],
            density=[1.845,1.845,1.845],
            nlenses_radii = [500e-4,1000e-4,1500e-4],
            lens_diameter=0.05,
            sigmaz=Sigma_x * 1e-4,
            alpha = alpha,
            sigmazp=Sigma_xp * 1e-6,
            tf_p=tf_p, tf_q=tf_q + qq,
            )

        #
        #  shadow
        #
        print(">>>>>>>>>>>>tf_p, tf_q, f: ", tf_p*1e-2, tf_q*1e-2, f*1e-2)
        beam, light_source = run_S4_ideal_lens(p_slit=27.3, tf_p=tf_p*1e-2, tf_q=tf_q*1e-2, focal=f*1e-2)
        # beam, light_source = run_S4_real_lens(p_slit=27.3, tf_p=tf_p * 1e-2, tf_q=tf_q * 1e-2, focal=f * 1e-2)
        print(">>>>>>>>>>>>sigmas_all: ", get_sigmas_all(light_source))
        from shadow4.tools.beamline_tools import focnew, focnew_scan

        ticket = focnew(beam=beam, nolost=1, mode=0, center=[0,0])
        y = qq * 1e-2
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

            FWHM_X[i] = tx['fwhm'] * 1e6
            FWHM_Z[i] = tz['fwhm'] * 1e6

            FWHM_FROM_STDEV_X[i] = 2.355 * 1e6 * bb.get_standard_deviation(1, nolost=1, ref=23)  # tx['fwhm'] * 1e6
            FWHM_FROM_STDEV_Z[i] = 2.355 * 1e6 * bb.get_standard_deviation(3, nolost=1, ref=23)  # tz['fwhm'] * 1e6



        print("Analytical size at sample H [um]: ", 2.355 * size_H[-1] * 1e4)
        print("Analytical size at sample V [um]: ", 2.355 * size[-1] * 1e4)
        #
        # recompute guess configuration
        #
        # nlenses_optimum = transfocator_compute_configuration(photon_energy_kev * 1e3,
        #                                                      size[-1],
        #                                                      symbol=["Be","Be","Be"],
        #                                                      density=[1.845,1.845,1.845],
        #                                                      nlenses_max = [15,3,1],
        #                                                      nlenses_radii = [500e-4,1000e-4,1500e-4],
        #                                                      lens_diameter=0.05,
        #                                                      sigmaz=Sigma_z * 1e-4,
        #                                                      sigmazp=Sigma_zp,
        #                                                      alpha = 1.0,
        #                                                      tf_p=tf_p,
        #                                                      tf_q=tf_q,
        #                                                      verbose=1 )
        #
        # print("To obtain ", size[-1], "cm in V we use: ", nlenses_optimum)

        #
        # (size, f, q_f, div) = transfocator_compute_parameters(
        #     photon_energy_kev*1e3,
        #     nlenses_optimum,
        #     symbol=["Be","Be","Be"],
        #     density=[1.845,1.845,1.845],
        #     nlenses_radii = [500e-4,1000e-4,1500e-4],
        #     lens_diameter=0.05,
        #     sigmaz=Sigma_x * 1e-4,
        #     alpha = 1.0,
        #     sigmazp=Sigma_zp * 1e-6,
        #     tf_p=tf_p,
        #     tf_q=tf_q,
        #     )
        #
        # print("To obtain ", size[-1], " in V we use: ", nlenses_optimum, "and got: ", size)


        #
        # plot
        #
        plot(
            qq * 1e-2, 2.355 * size_H * 1e4,
            qq * 1e-2, 2.355 * size * 1e4,
            y,  FWHM_X, # ylist[0],
            y,  FWHM_Z, # ylist[1],
            y, FWHM_FROM_STDEV_X,  # ylist[0],
            y, FWHM_FROM_STDEV_Z,  # ylist[1],
            # y, 2.355 * ylist[0],
            # y, 2.355 * ylist[1],
            title="FWHM",
            xtitle="Y [m]",
            ytitle="size [$\mu$m]",
            legend=["X analytical","Z analytical", \
                "X ray tracing from histo", "Z ray tracing from histo", \
                "X ray tracing from stdev","Z ray tracing from stdev"],
            linestyle=[None,None,':',':','--','--'],
            color=['b', 'r', 'b', 'r', 'b', 'r',],
            grid=1)