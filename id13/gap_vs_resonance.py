

def get_E_K(do_plot=1, energies=None, ns=None):
    import numpy
    import scipy.constants as codata
    from srxraylib.plot.gol import plot

    #
    # inputs
    #
    gap_min                = 4 # mm
    gap_max                = 30 # mm
    period_length          = 0.0164 # m
    number_of_periods      = 121.951
    ring_current           = 0.2 # A
    electron_energy_in_GeV = 6 # GeV
    pow_dens_screen        = 31.1 # distance to screen in m
    auto_harmonic_number   = 1
    A = [3.1777, 0.0, 0.9846, 1.0621, 1.0, 1.0992]

    #
    # calculations
    #
    gap_mm = numpy.linspace(gap_min * 0.9, gap_max * 1.1, 1000)
    i_half = len(A) // 2

    # get K vs gap
    Bmax = numpy.zeros_like(gap_mm)
    for i in range(i_half):
        Bmax += A[i] * numpy.exp(-numpy.pi * (i + 1) * A[i + i_half] * gap_mm / (period_length * 1e3))

    Karray = Bmax * period_length * codata.e / (2 * numpy.pi * codata.m_e * codata.c)

    print(">>>> Karray: ", Karray)
    # resonance energy
    gamma1 = 1e9 * electron_energy_in_GeV / (codata.m_e *  codata.c**2 / codata.e)
    Karray_horizontal = numpy.zeros_like(Karray)
    Bfield = Bmax # Karray / (period_length * codata.e / (2 * numpy.pi * codata.m_e * codata.c))

    theta_x = 0.0
    theta_z = 0.0
    wavelength = (period_length / (2.0 * gamma1 ** 2)) * (1 + Karray ** 2 / 2.0 + Karray_horizontal ** 2 / 2.0 + gamma1 ** 2 * (theta_x ** 2 + theta_z ** 2))
    wavelength /= auto_harmonic_number
    frequency = codata.c / wavelength * auto_harmonic_number
    energy_in_ev = codata.h * frequency / codata.e
    E1_array = energy_in_ev * auto_harmonic_number

    # power
    ptot = (number_of_periods / 6) * codata.value('characteristic impedance of vacuum') * ring_current * codata.e * 2 * numpy.pi * codata.c * gamma1 ** 2 * (Karray ** 2 + Karray_horizontal ** 2) / period_length

    ### From: Undulators, Wigglers and their Applications - H. Onuki & P Elleaume ###
    ### Chapter 3: Undulator radiation eqs. 68 and 69 ###
    g_k = Karray * ((Karray ** 6) + (24 / 7) * (Karray ** 4) + 4 * (Karray ** 2) + (16 / 7)) / ((1 + (Karray ** 2)) ** (7 / 2))

    p_dens_peak = ((21 * (gamma1 ** 2)) / (16 * numpy.pi * Karray) * ptot * g_k) / ((pow_dens_screen * 1e3) ** 2)


    #
    # inrerpolation
    ENERGIES = []
    KS = []
    if energies is None:
        energies = [5800, 7800, 10000, 11000, 12700, 17000, 20000, 20000, 23500, 25000, 30000, 35000, 40000]
        ns = [1,1,1,1,1,1,1,3,3,3,3,3,3]

    for i, energy in enumerate(energies):
        gap = numpy.interp(energy / ns[i], E1_array, gap_mm)
        print("interpolated (n=%d) energy: %g keV, gap: %.2g mm, K: %.5g, Ptot: %.3g kW, Ppeak@screen: %.3g W/mm2: " %
              (ns[i], energy *1e-3,
              gap,
              numpy.interp(gap, gap_mm, Karray),
              1e-3 * numpy.interp(gap, gap_mm, ptot),
              numpy.interp(gap, gap_mm, p_dens_peak),
              ))

        ENERGIES.append(energy)
        KS.append(numpy.interp(gap, gap_mm, Karray))


    #
    # plots
    #
    if do_plot:
        plot(gap_mm, Karray, title="K vs Gap", xtitle="Gap [mm]", ytitle="K", show=0)
        plot(gap_mm, Bfield, title="B vs Gap", xtitle="Gap [mm]", ytitle="B [T]", show=0)
        plot(E1_array, gap_mm,
                  E1_array * 3, gap_mm,
                  E1_array * 5, gap_mm,
                  title="Gap vs resonance energy", xtitle="Photon energy [eV]", ytitle="Gap [mm]",
                  legend=['harmonic 1', 'harmonic 3', 'harmonic 5'], show=0)
        plot(gap_mm, ptot,        title="Power vs Gap", xtitle="Gap [mm]", ytitle="Power [W]", show=0)
        plot(gap_mm, p_dens_peak, title="Power density peak at screen vs Gap", xtitle="Gap [mm]", ytitle="Power density peak at screen [W/mm2]", show=1)

    return numpy.array(ENERGIES), numpy.array(KS)

def run_shadow4(Ei, Ki):
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam
    electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
    electron_beam.set_sigmas_all(sigma_x=3.33432e-05, sigma_y=7.28361e-06, sigma_xp=4.52892e-06, sigma_yp=1.93718e-06)

    # magnetic structure
    from shadow4.sources.undulator.s4_undulator import S4Undulator
    source = S4Undulator(
        K_vertical=Ki, #2.277,  # syned Undulator parameter
        period_length=0.0164,  # syned Undulator parameter
        number_of_periods=121.951,  # syned Undulator parameter
        emin=Ei, # 5800.0,  # Photon energy scan from energy (in eV)
        emax=Ei, # 5800.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=2.5e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='internal',  # internal, pysru, srw
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=1,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.05,  # for SRW backpropagation, the range factor
        srw_resolution=50,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=0,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.05,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=0,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=0,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.5,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    beamline.set_light_source(light_source)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.001, x_right=0.001, y_bottom=-0.001, y_top=0.001)

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator PRIMARY SLITS',
                               boundary_shape=boundary_shape,
                               i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                               i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=27.2, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.00035, x_right=0.00035, y_bottom=-0.00035, y_top=0.00035)

    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator MONO SLITS',
                               boundary_shape=boundary_shape,
                               i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                               i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=3.9, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
                                     boundary_shape=None, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=Ei, # 5802.7595568424,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.223060463, angle_azimuthal=0,
                                     angle_radial_out=1.223060463)
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Plane Crystal',
                                     boundary_shape=None, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=Ei, # 5802.7595568424,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0,
                                     # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0,  # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.223060417, angle_azimuthal=3.141592654,
                                     angle_radial_out=1.223060417)
    movements = None
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                             movements=movements, input_beam=beam)

    beam, mirr = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # test plot
    if 0:
        from srxraylib.plot.gol import plot_scatter
        plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
                     title='(Intensity,Photon Energy)', plot_histograms=0)
        plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')

    return beam

if __name__ == "__main__":
    E, K = get_E_K(do_plot=0)

    GEOMET_RATIO = []
    INTENS_RATIO = []

    for i in range(E.size):
        beam = run_shadow4(Ei=E[i], Ki=K[i])
        GEOMET_RATIO.append(beam.get_number_of_rays(nolost=1) / beam.get_number_of_rays(nolost=0))
        INTENS_RATIO.append(beam.get_intensity(nolost=1) / beam.get_intensity(nolost=0)          )


    print(">>>> E[i], K[i], GEOMET_RATIO[i], INTENS_RATIO[i] <<<<<")
    for i in range(E.size):
        print("%.3f %.3f  %.3f  %.3f" % (E[i], K[i], GEOMET_RATIO[i], INTENS_RATIO[i]))

