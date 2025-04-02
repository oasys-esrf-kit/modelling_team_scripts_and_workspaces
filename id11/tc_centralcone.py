
import numpy
import scipy.constants as codata
# from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource


def Fn(K, n):
    """
    Calculate the universal function Fn(x).

    Parameters
    ----------
    K : float
        The K value
    n : int, optional
        The harmonic number (it must be odd).

    Returns
    -------
    float
        Fn(x)

    References
    ----------
    see fig 2.4 in X-RAY DATA BOOKLET https://xdb.lbl.gov/

    """
    from scipy.special import jn, yn, jv, yv
    if n % 2 == 1:
        cst1 = ((n * K) / (1. + (K ** 2) / 2.)) ** 2
        cst2 = (n * K ** 2) / (4. + (2. * K ** 2))
        Fn = cst1 * (jn(0.5 * (n - 1), cst2) - jn(0.5 * (n + 1), cst2)) ** 2
    else:
        Fn = 0.0
    return Fn


def Qn(K, n):
    """
    Calculate the universal function Qn(x).

    Parameters
    ----------
    K : float
        The K value
    n : int, optional
        The harmonic number (it must be odd).

    Returns
    -------
    float
        Qn(x)

    References
    ----------
    see fig 2.6 in X-RAY DATA BOOKLET https://xdb.lbl.gov/

    """

    if n == 0:
        raise Exception(' the harmonic number can not be 0')
    res = (1. + 0.5 * K ** 2) * Fn(K, n) / n
    return res


def calculate_tc_centralcone(electron_energy_in_GeV=6.0, current=0.2,
    kmin=0.001, kmax=1.832, kpoints=10, period_length=0.018, number_of_periods=111.111, harmonics=[1,3,5,7]):

    karray = numpy.linspace(kmin, kmax, kpoints)
    gamma1 = 1e9 * electron_energy_in_GeV / (codata.m_e * codata.c ** 2 / codata.e)
    wavelength = (period_length / (2.0 * gamma1 ** 2)) * (1 + karray ** 2 / 2.0)
    energy1_in_ev = codata.h * codata.c / wavelength / codata.e

    out = []
    for n in harmonics:
        F = (numpy.pi * codata.alpha * 1e-3 / codata.e * number_of_periods * current * Qn(karray, n))
        out.append(F)

    return karray, energy1_in_ev, out

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.size': 22})


    kmin=0.001
    kpoints=100
    electron_energy_in_GeV=6.0
    current=0.2
    harmonics=[1,3,5,7,9,11]

    #
    # u18
    #
    period_length=0.018
    number_of_periods=111.111
    kmax_g5p5=1.832
    kmax_g6=1.656

    _, energy1_in_ev_u18_g5p5, u18_g5p5 = calculate_tc_centralcone(kmin=kmax_g6, kmax=kmax_g5p5, kpoints=kpoints,
                                                     electron_energy_in_GeV=electron_energy_in_GeV, current=current,
                                                     harmonics=harmonics,
                                                     period_length=period_length, number_of_periods=number_of_periods)

    _, energy1_in_ev_u18_g6, u18_g6 = calculate_tc_centralcone(kmin=kmin, kmax=kmax_g6, kpoints=kpoints,
                                                     electron_energy_in_GeV=electron_energy_in_GeV, current=current,
                                                     harmonics=harmonics,
                                                     period_length=period_length, number_of_periods=number_of_periods)

    #
    # u20.5
    #
    period_length=0.0205
    number_of_periods=97.561
    kmax_g5p5=2.557
    kmax_g6=2.334

    _, energy1_in_ev_u20p5_g5p5, u20p5_g5p5 = calculate_tc_centralcone(kmin=kmax_g6, kmax=kmax_g5p5, kpoints=kpoints,
                                                     electron_energy_in_GeV=electron_energy_in_GeV, current=current,
                                                     harmonics=harmonics,
                                                     period_length=period_length, number_of_periods=number_of_periods)

    _, energy1_in_ev_u20p5_g6, u20p5_g6 = calculate_tc_centralcone(kmin=kmin, kmax=kmax_g6, kpoints=kpoints,
                                                     electron_energy_in_GeV=electron_energy_in_GeV, current=current,
                                                     harmonics=harmonics,
                                                     period_length=period_length, number_of_periods=number_of_periods)


    from srxraylib.plot.gol import plot
    plot(
         1e-3 * energy1_in_ev_u18_g6 * 1, u18_g6[0],
         1e-3 * energy1_in_ev_u18_g6 * 3, u18_g6[1],
         1e-3 * energy1_in_ev_u18_g6 * 5, u18_g6[2],
         1e-3 * energy1_in_ev_u18_g6 * 7, u18_g6[3],
         1e-3 * energy1_in_ev_u18_g6 * 9, u18_g6[4],
         1e-3 * energy1_in_ev_u18_g6 * 11, u18_g6[5],
         1e-3 * energy1_in_ev_u18_g5p5 * 1, u18_g5p5[0],
         1e-3 * energy1_in_ev_u18_g5p5 * 3, u18_g5p5[1],
         1e-3 * energy1_in_ev_u18_g5p5 * 5, u18_g5p5[2],
         1e-3 * energy1_in_ev_u18_g5p5 * 7, u18_g5p5[3],
         1e-3 * energy1_in_ev_u18_g5p5 * 9, u18_g5p5[4],
         1e-3 * energy1_in_ev_u18_g5p5 * 11, u18_g5p5[5],
        ############
         1e-3 * energy1_in_ev_u20p5_g6 * 1, u20p5_g6[0],
         1e-3 * energy1_in_ev_u20p5_g6 * 3, u20p5_g6[1],
         1e-3 * energy1_in_ev_u20p5_g6 * 5, u20p5_g6[2],
         1e-3 * energy1_in_ev_u20p5_g6 * 7, u20p5_g6[3],
         1e-3 * energy1_in_ev_u20p5_g6 * 9, u20p5_g6[4],
         1e-3 * energy1_in_ev_u20p5_g6 * 11, u20p5_g6[5],
         1e-3 * energy1_in_ev_u20p5_g5p5 * 1, u20p5_g5p5[0],
         1e-3 * energy1_in_ev_u20p5_g5p5 * 3, u20p5_g5p5[1],
         1e-3 * energy1_in_ev_u20p5_g5p5 * 5, u20p5_g5p5[2],
         1e-3 * energy1_in_ev_u20p5_g5p5 * 7, u20p5_g5p5[3],
         1e-3 * energy1_in_ev_u20p5_g5p5 * 9, u20p5_g5p5[4],
         1e-3 * energy1_in_ev_u20p5_g5p5 * 11, u20p5_g5p5[5],
         color=['b','b','b','b','b','b','r','r','r','r','r','r','k','k','k','k','k','k','r','r','r','r','r','r',],
         ytitle="Flux [photons/s/0.1%bw] in central cone", xtitle="Photon Energy [keV]",
        ylog=1, yrange=[1e11,1e16], xrange=[0,150], figsize=(12,8), show=0)

    plt.grid()

    plot(
         1e-3 * energy1_in_ev_u18_g6 * 1, u18_g6[0],
         1e-3 * energy1_in_ev_u18_g6 * 3, u18_g6[1],
         1e-3 * energy1_in_ev_u18_g6 * 5, u18_g6[2],
         1e-3 * energy1_in_ev_u18_g5p5 * 1, u18_g5p5[0],
         1e-3 * energy1_in_ev_u18_g5p5 * 3, u18_g5p5[1],
         1e-3 * energy1_in_ev_u18_g5p5 * 5, u18_g5p5[2],
        ############
         1e-3 * energy1_in_ev_u20p5_g6 * 1, u20p5_g6[0],
         1e-3 * energy1_in_ev_u20p5_g6 * 3, u20p5_g6[1],
         1e-3 * energy1_in_ev_u20p5_g6 * 5, u20p5_g6[2],
         1e-3 * energy1_in_ev_u20p5_g5p5 * 1, u20p5_g5p5[0],
         1e-3 * energy1_in_ev_u20p5_g5p5 * 3, u20p5_g5p5[1],
         1e-3 * energy1_in_ev_u20p5_g5p5 * 5, u20p5_g5p5[2],
         color=['b','b','b','r','r','r','k','k','k','r','r','r',],
         ytitle="Flux [photons/s/0.1%bw] in central cone", xtitle="Photon Energy [keV]",
        ylog=1, yrange=[1e13,1e16], xrange=[0,50], figsize=(12,8), show=0)

    plt.grid()

    plt.show()