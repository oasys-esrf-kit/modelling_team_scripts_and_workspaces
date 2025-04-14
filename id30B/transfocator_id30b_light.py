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

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_show

    Sigma_x  =  28.2381e-6 # m
    Sigma_z  =  6.03339e-6 # m
    Sigma_xp =  7.20068e-6 # rad
    Sigma_zp =  5.51623e-6 # rad

    tf_p = 62.052
    tf_q = 39.378

    photon_energy_kev = float(input("Enter photon energy in keV: "))
    s_target_fwhm_um = float(input("Enter target beam size at sample in um FWHM: "))



    #
    # Guess config
    #
    if 1:
        s_target_V = s_target_fwhm_um / 2.355 * 1e-6
        s_target_H = s_target_fwhm_um / 2.355 * 1e-6

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


    #
    # plot
    #
    qq = numpy.linspace(-40.0, 0, 100)

    (size_VV, f, q_f, div) = transfocator_compute_parameters(
        photon_energy_kev * 1e3,
        nlenses_optimum_V,
        symbol=["Be","Be","Be"],
        density=[1.845,1.845,1.845],
        nlenses_radii = [500e-6,1000e-6,1500e-6],
        sigmaz=Sigma_z,
        sigmazp=Sigma_zp,
        tf_p=tf_p,
        tf_q=tf_q + qq,
        )

    (size_VH, f_H, q_f_H, div_H) = transfocator_compute_parameters(
        photon_energy_kev * 1e3,
        nlenses_optimum_V,
        symbol=["Be","Be","Be"],
        density=[1.845,1.845,1.845],
        nlenses_radii = [500e-6,1000e-6,1500e-6],
        sigmaz=Sigma_x,
        sigmazp=Sigma_xp,
        tf_p=tf_p, tf_q=tf_q + qq,
        )

    (size_HV, f, q_f, div) = transfocator_compute_parameters(
        photon_energy_kev * 1e3,
        nlenses_optimum_H,
        symbol=["Be","Be","Be"],
        density=[1.845,1.845,1.845],
        nlenses_radii = [500e-6,1000e-6,1500e-6],
        sigmaz=Sigma_z,
        sigmazp=Sigma_zp,
        tf_p=tf_p,
        tf_q=tf_q + qq,
        )

    (size_HH, f_H, q_f_H, div_H) = transfocator_compute_parameters(
        photon_energy_kev * 1e3,
        nlenses_optimum_H,
        symbol=["Be","Be","Be"],
        density=[1.845,1.845,1.845],
        nlenses_radii = [500e-6,1000e-6,1500e-6],
        sigmaz=Sigma_x,
        sigmazp=Sigma_xp,
        tf_p=tf_p, tf_q=tf_q + qq,
        )



    print("Analytical size optimized for V at sample H [um]: ", 2.355 * size_VH[-1] * 1e6)
    print("Analytical size optimized for V at sample V [um]: ", 2.355 * size_VV[-1] * 1e6)

    print("Analytical size optimized for H at sample H [um]: ", 2.355 * size_HH[-1] * 1e6)
    print("Analytical size optimized for H at sample V [um]: ", 2.355 * size_HV[-1] * 1e6)

    plot(
        qq, 2.355 * size_VH * 1e6,
        qq, 2.355 * size_VV * 1e6,
        qq, 2.355 * size_HH * 1e6,
        qq, 2.355 * size_HV * 1e6,
        title="FWHM",
        xtitle="Y [m]",
        ytitle="size FWHM [$\mu$m]",
        legend=["optimized for V %s, V analytical" % repr(nlenses_optimum_V), "optimized for V %s, H analytical" % repr(nlenses_optimum_V),
                "optimized for H %s, V analytical" % repr(nlenses_optimum_H), "optimized for H %s, H analytical" % repr(nlenses_optimum_H)],
        linestyle=[None, None, ':', ':'],
        color=['b', 'r', 'b', 'r'],
        grid=1, show=0)

    plot_show()