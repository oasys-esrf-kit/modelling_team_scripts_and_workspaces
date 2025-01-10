
import numpy
from crystalpy.util.calc_xcrystal import calc_xcrystal_angular_scan, calc_xcrystal_energy_scan, calc_xcrystal_alphazachariasen_scan

def get_fwhm(histogram, bins, ret0=None):
    fwhm = ret0
    quote = ret0
    coordinates = None

    if histogram.size > 1:
        quote = numpy.max(histogram)*0.5
        cursor = numpy.where(histogram >= quote)

        if histogram[cursor].size > 1:
            bin_size    = bins[1]-bins[0]
            fwhm        = bin_size*(cursor[0][-1]-cursor[0][0])
            coordinates = (bins[cursor[0][0]], bins[cursor[0][-1]])

    return fwhm, quote, coordinates

def calculate(energy=7000.0, nreflections=1):
    bunch_out_dict, diffraction_setup, deviations = calc_xcrystal_angular_scan(
        # material_constants_library_flag=self.material_constants_library_flag,
        crystal_name              = 'Si',
        thickness                 = 0.007,
        miller_h                  = 1,
        miller_k                  = 1,
        miller_l                  = 1,
        asymmetry_angle           = 0.0,
        energy                    = energy,
        angle_deviation_min       = -0.0001 * 0.35, # 0.5, 1.0
        angle_deviation_max       = 0.00015 * 0.35, # 0.5, 1.0
        angle_deviation_points    = 1000,
        angle_center_flag         = 2,
        calculation_method        = 1, # 0=Zachariasen, 1=Guigay
        is_thick                  = 0,
        use_transfer_matrix       = 0,
        geometry_type_index       = 0,
        calculation_strategy_flag = 0, # 0=mpmath 1=numpy 2=numpy-truncated
                )

    tmp = numpy.zeros((bunch_out_dict["energies"].size,7))
    tmp[:, 0] = deviations / 1e-06
    tmp[:, 1] = 7000.0
    tmp[:, 2] = bunch_out_dict["phaseP"] * nreflections
    tmp[:, 3] = bunch_out_dict["phaseS"] * nreflections
    # tmp[:, 4] = circular polarization
    tmp[:, 5] = bunch_out_dict["intensityP"] ** nreflections
    tmp[:, 6] = bunch_out_dict["intensityS"] ** nreflections


    return tmp, diffraction_setup.angleBragg(energy)

if __name__ == "__main__":
    from srxraylib.plot.gol import plot

    energy = 40000.0
    tmp, thetaB = calculate(energy=energy, nreflections=1)

    fwhm_s, _, _ = get_fwhm(tmp[:, 6], tmp[:, 0])
    fwhm_p, _, _ = get_fwhm(tmp[:, 5], tmp[:, 0])
    peak_s = tmp[:, 6].max()
    peak_p = tmp[:, 5].max()


    peak_s = tmp[:, 6].max()
    peak_p = tmp[:, 5].max()

    integral_s = numpy.trapz( tmp[:, 6], tmp[:, 0])
    integral_p = numpy.trapz( tmp[:, 5], tmp[:, 0])

    deltaE_s = energy * 1e-6 * fwhm_s / numpy.tan(thetaB)
    deltaE_p = energy * 1e-6 * fwhm_p / numpy.tan(thetaB)

    print("Photon energy: %.1f keV, Bragg angle: %.1f deg" % (1e-3 * energy, numpy.degrees(thetaB)))

    print("DE[eV] %.3f %.3f"                               % ( deltaE_s, deltaE_p))
    print("FWHM[urad] %.1f %.1f"                           % ( fwhm_s,     fwhm_p))
    print("PEAK %.2f %.2f"                                 % ( peak_s,     peak_p))
    print("INTEGRAL[urad] %.3f %.3f"                       % ( integral_s, integral_p))

    print("\n%.3f %.3f" % ( deltaE_s, deltaE_p))
    print("%.1f %.1f" % ( fwhm_s,     fwhm_p))
    print("%.2f %.2f" % ( peak_s,     peak_p))
    print("%.3f %.3f" % ( integral_s, integral_p))


    plot(tmp[:,0], tmp[:,6], tmp[:,0], tmp[:,5], xtitle=r"$\theta$-$\theta_B$ [$\mu$ rad]", legend=["S-pol","P-pol"])

