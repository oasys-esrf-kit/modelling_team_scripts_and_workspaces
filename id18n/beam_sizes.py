


def run_S4_source(photon_energy=7000.0, harmonic_number=1):
    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam
    electron_beam = S4ElectronBeam(energy_in_GeV=6,energy_spread=0.001,current=0.2)
    electron_beam.set_sigmas_all(sigma_x=2.87429e-05, sigma_y=5.15531e-06, sigma_xp=4.18025e-06, sigma_yp=1.93975e-06)

    # magnetic structure
    from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
    source = S4UndulatorGaussian(
        period_length     = 0.0186,     # syned Undulator parameter (length in m)
        number_of_periods = 134.5, # syned Undulator parameter
        photon_energy     = photon_energy, # Photon energy (in eV)
        delta_e           = 0.0, # Photon energy width (in eV)
        ng_e              = 100, # Photon energy scan number of points
        flag_emittance    = 1, # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread = 1, # when sampling rays: Use e- energy spread (0=No, 1=Yes)
        harmonic_number    = harmonic_number, # harmonic number
        flag_autoset_flux_central_cone  = 1, # value to set the flux peak
        flux_central_cone  = 3230278823662982.0, # value to set the flux peak
        )


    # light source
    from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
    light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam, magnetic_structure=source,nrays=5000,seed=5676561)
    # beam = light_source.get_beam()
    #
    # # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')

    return light_source

if __name__ == "__main__":
    import numpy

    ENERGY = [7000, 17000, 35000, 40000]
    N = [1, 3, 5, 5]

    fwhm = True

    for i in range(len(N)):
        ls = run_S4_source(photon_energy=ENERGY[i], harmonic_number=N[i])
        sigma_x, sigdi_x, sigma_z, sigdi_z = ls.get_electron_beam().get_sigmas_all()
        us, ua = ls.get_undulator_photon_beam_sizes(undulator_E0=ENERGY[i], undulator_length=2.5)
        x = ls.norm_energ_spr(harmonic_number=N[i])
        qs = ls.q_s(x, factor=0.5)
        qa = ls.q_a(x)
        uus = us * numpy.sqrt(qs**2 - 1.0)
        uua = ua * numpy.sqrt(qa**2 - 1.0)


        if fwhm:
            ff = 2.355
        else:
            ff = 1.00

        print("\n=============== photon energy: ", ENERGY[i], "==========================")
        if fwhm:
            print("Values are FWHM")
        else:
            print("Values are SIGMA")
        print("Hsize %.3f %.3f %.3f %.3f" % ( ff * 1e6 * us, ff * 1e6 * sigma_x, ff * 1e6 * uus, ff * 1e6 * numpy.sqrt(us ** 2 + sigma_x ** 2 + uus ** 2)))
        print("Vsize %.3f %.3f %.3f %.3f" % ( ff * 1e6 * us, ff * 1e6 * sigma_z, ff * 1e6 * uus, ff * 1e6 * numpy.sqrt(us ** 2 + sigma_z ** 2 + uus ** 2)))
        print("Hangle %.3f %.3f %.3f %.3f" % ( ff * 1e6 * ua, ff * 1e6 * sigdi_x, ff * 1e6 * uua, ff * 1e6 * numpy.sqrt(ua ** 2 + sigdi_x ** 2 + uua ** 2)))
        print("Vangle %.3f %.3f %.3f %.3f" % ( ff * 1e6 * ua, ff * 1e6 * sigdi_z, ff * 1e6 * uua, ff * 1e6 * numpy.sqrt(ua ** 2 + sigdi_z ** 2 + uua ** 2)))

