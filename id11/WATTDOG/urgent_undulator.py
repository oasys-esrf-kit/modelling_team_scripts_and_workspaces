# from srxraylib.plot.gol import plot, plot_image
import scipy.constants as codata
import numpy

"""
urgent_undulator.py
Exact implementation of Walker & Diviacco, Rev. Sci. Instrum. 63, 392 (1992)

PLANE UNDULATOR (Kx=0), Section VI Power Density.

Amplitude (from paper, Kx=0):
    A = xi * ( 2*alpha_x*S0 - Ky*(S1 + S_minus1),
               2*alpha_y*S0,
               0 )

where:
    S_q = sum_{p=-inf}^{inf} J_p(Y) * J_{2p+q+n}(X)    q = 0, +1, -1

    xi    = n / (1 + Ky^2/2 + alpha^2)
    X     = 2 * xi * alpha_x * Ky          (alpha_x = gamma*theta_x)
    Y     = xi * Ky^2 / 4                  (Kx=0 -> Y = xi*(Ky^2-0)/4)

Power density (Section VI):
    dI/dOmega = (e*gamma^4*Ib) / (eps0*lambda0*N)
                * |A_n|^2 / (1 + Ky^2/2 + alpha^2)
"""

import numpy as np
from scipy.special import jv

ECHARGE = 1.60217663e-19
EPS0    = 8.85418782e-12
E0_EV   = 0.51099895e6


def harmonic_energy(n, K, period_m, gamma):
    hc_eVm = 1239.84193e-9
    return n * 2.0 * gamma**2 * hc_eVm / (period_m * (1.0 + K**2 / 2.0))


def _Sq(X, Y, n, q, p_max=20):
    """
    S_q = sum_{p=-pmax}^{pmax} J_p(Y) * J_{2p+q+n}(X)
    X, Y are 2D arrays; returns 2D array.
    """
    S = np.zeros_like(X, dtype=float)
    for p in range(-p_max, p_max + 1):
        S += jv(p, Y) * jv(2*p + q + n, X)
    return S


def power_density_harmonic(n, Ky, N_periods, period_m, gamma, current_A,
                            theta_x_rad, theta_y_rad, p_max=20):
    """Power density [W/mrad^2] for harmonic n. Walker 1992, plane undulator."""
    Ky2     = Ky**2
    alpha_x = gamma * theta_x_rad
    alpha_y = gamma * theta_y_rad
    alpha2  = alpha_x**2 + alpha_y**2

    denom = 1.0 + Ky2 / 2.0 + alpha2
    xi    = n / denom
    X     = 2.0 * xi * alpha_x * Ky
    Y     = xi * Ky2 / 4.0               # positive: Y = xi*(Ky^2-Kx^2)/4, Kx=0

    S0  = _Sq(X, Y, n,  0, p_max)
    S1  = _Sq(X, Y, n,  1, p_max)
    Sm1 = _Sq(X, Y, n, -1, p_max)

    # A components (Kx=0)
    Ax = xi * (2.0 * alpha_x * S0 - Ky * (S1 + Sm1))
    Ay = xi * (2.0 * alpha_y * S0)

    An2 = Ax**2 + Ay**2

    # Physical prefactor (Walker 1992, Section VI)
    # prefactor = ECHARGE * gamma**4 * current_A / (EPS0 * period_m * N_periods)
    # Try this instead:
    # prefactor = ECHARGE * gamma ** 4 * current_A / (4 * np.pi * EPS0 * period_m * N_periods)
    prefactor1 = gamma ** 4 * current_A / (period_m * N_periods) # from parameters
    prefactor2 = 1.0 # constants, TO DO

    return prefactor1 * prefactor2 * An2 / denom


def power_density_all_harmonics(Ky, N_periods, period_m, gamma, current_A,
                                 theta_x_rad, theta_y_rad,
                                 n_harmonics=200, p_max=20,
                                 return_individual=True):
    total = np.zeros_like(theta_x_rad, dtype=float)
    indiv = {}
    for n in range(1, n_harmonics + 1):
        pd_n = power_density_harmonic(n, Ky, N_periods, period_m,
                                       gamma, current_A,
                                       theta_x_rad, theta_y_rad, p_max)
        total += pd_n
        if return_individual:
            indiv[n] = pd_n
        if n % 20 == 0:
            print(f"  harmonic {n}/{n_harmonics}")
    return (total, indiv) if return_individual else total

def run_urgent(h5_parameters, n_harmonics=10, do_plot=0, show=1, ):

    from xoppylib.sources.srundplug import calc2d_from_harmonics_urgent

    bl = dict()
    bl['ElectronBeamDivergenceH'] = h5_parameters["ELECTRONBEAMDIVERGENCEH"]
    bl['ElectronBeamDivergenceV'] = h5_parameters["ELECTRONBEAMDIVERGENCEV"]
    bl['ElectronBeamSizeH']       = h5_parameters["ELECTRONBEAMSIZEH"]
    bl['ElectronBeamSizeV']       = h5_parameters["ELECTRONBEAMSIZEV"]
    bl['ElectronCurrent']         = h5_parameters["ELECTRONCURRENT"]
    bl['ElectronEnergy']          = h5_parameters["ELECTRONENERGY"]
    bl['ElectronEnergySpread']    = h5_parameters["ELECTRONENERGYSPREAD"]
    bl['Kv']                      = h5_parameters["KV"]
    bl['Kh']                      = h5_parameters["KH"]
    bl['Kphase']                  = h5_parameters["KPHASE"]
    bl['NPeriods']                = h5_parameters["NPERIODS"]
    bl['PeriodID']                = h5_parameters["PERIODID"]
    bl['distance']                = h5_parameters["DISTANCE"]
    bl['gapH']                    = h5_parameters["GAPH"]
    bl['gapV']                    = h5_parameters["GAPV"]

    # horizontal, vertical, power_density = calc2d_urgent(bl, zero_emittance=False, fileName=None, fileAppend=False,
    #                                     hSlitPoints=21, vSlitPoints=51)


    # X1, Y1, POWER_DENSITY1 = calc2d_urgent(bl,
    #                             zero_emittance=False, fileName=None, fileAppend=False,
    #                             hSlitPoints=h5_parameters["HSLITPOINTS"],
    #                             vSlitPoints=h5_parameters["VSLITPOINTS"],)

    X, Y, POWER_DENSITY, POWER_DENSITY_HARMONICS, ENERGY_HARMONICS, FLUX = calc2d_from_harmonics_urgent(bl,
                                zero_emittance=False, fileName=None, fileAppend=False,
                                hSlitPoints=h5_parameters["HSLITPOINTS"],
                                vSlitPoints=h5_parameters["VSLITPOINTS"],
                                harmonic_max=n_harmonics,
                                return_flux=1)

    print(POWER_DENSITY_HARMONICS.shape, ENERGY_HARMONICS.shape, FLUX.shape, X.shape, Y.shape)
    dx = X[1] - X[0]
    dy = Y[1] - Y[0]
    print("\nharm, power, power-from-flux, power-density-peak, flux: ")
    for i in range(POWER_DENSITY_HARMONICS.shape[0]):
        print("%d  %10.3f  %10.3f  %10.3f  %g" % (i,
              POWER_DENSITY_HARMONICS[i].sum() * dx * dy,
              (FLUX * ENERGY_HARMONICS * codata.e)[i].sum() * dx * dy, # inefficient, but correct
              POWER_DENSITY_HARMONICS[i].max() * dx * dy,
              FLUX[i].sum() * dx * dy,)
              )

    # read_urgent("urgent.out")
    # example plot
    if do_plot:
        from srxraylib.plot.gol import plot_image, plot_show

        print(POWER_DENSITY.max(), POWER_DENSITY.shape, X.shape, Y.shape)
        plot_image(POWER_DENSITY, X, Y, xtitle="H [mm]", ytitle="V [mm]", title="Power density W/mm2 **FROM HARMONICS**", show=0)

        plot_image(FLUX.sum(axis=0), X, Y, xtitle="H [mm]", ytitle="V [mm]", title="Flux  photons/mm2/? **FROM HARMONICS**", show=0)


        # print(POWER_DENSITY1.max(), POWER_DENSITY1.shape, X1.shape, Y1.shape)
        # plot_image(POWER_DENSITY1, X1, Y1, xtitle="H [mm]", ytitle="V [mm]", title="Power density W/mm2", show=0)

        if do_plot:
            if show: plot_show()

    return X, Y, POWER_DENSITY, POWER_DENSITY_HARMONICS, ENERGY_HARMONICS, FLUX


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import time

    n_harmonics = 10 # 24
    calculate_python  = 1
    calculate_fortran = 1

##

    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"] = 6.0
    h5_parameters["ELECTRONENERGYSPREAD"] = 0.001
    h5_parameters["ELECTRONCURRENT"] = 0.2
    h5_parameters["ELECTRONBEAMSIZEH"] = 3.34281e-05
    h5_parameters["ELECTRONBEAMSIZEV"] = 7.28139e-06
    h5_parameters["ELECTRONBEAMDIVERGENCEH"] = 4.51097e-06
    h5_parameters["ELECTRONBEAMDIVERGENCEV"] = 1.94034e-06
    h5_parameters["PERIODID"] = 0.018
    h5_parameters["NPERIODS"] = 111
    h5_parameters["KV"] = 1.358
    h5_parameters["KH"] = 0.0
    h5_parameters["KPHASE"] = 0.0
    h5_parameters["DISTANCE"] = 30.0
    h5_parameters["GAPH"] = 0.01
    h5_parameters["GAPV"] = 0.01
    h5_parameters["HSLITPOINTS"] = 48
    h5_parameters["VSLITPOINTS"] = 30
    h5_parameters["METHOD"] = 1
    h5_parameters["USEEMITTANCES"] = 1
    h5_parameters["MASK_FLAG"] = 0
    h5_parameters["MASK_ROT_H_DEG"] = 0.0
    h5_parameters["MASK_ROT_V_DEG"] = 0.0
    h5_parameters["MASK_H_MIN"] = -1000.0
    h5_parameters["MASK_H_MAX"] = 1000.0
    h5_parameters["MASK_V_MIN"] = -1000.0
    h5_parameters["MASK_V_MAX"] = 1000.0

    ##
    period_m  = h5_parameters["PERIODID"]
    Ky        = h5_parameters["KV"]
    N_periods = h5_parameters["NPERIODS"]
    current_A = h5_parameters["ELECTRONCURRENT"]
    E_GeV     = h5_parameters["ELECTRONENERGY"]

    #
    # direct calculation
    #
    Z_m       = h5_parameters["DISTANCE"]
    B0 = Ky / (0.9336 * period_m * 1e2)
    gamma = 1e9 * E_GeV / (codata.m_e *  codata.c**2 / codata.e)

    print(f"K  = {Ky:.4f}")
    print(f"E1 = {harmonic_energy(1, Ky, period_m, gamma)/1e3:.3f} keV")
    print(f"B0 = {B0:.3f} T")

    if calculate_python:
        x_m = np.linspace(-0.5 * h5_parameters["GAPH"], 0.5 * h5_parameters["GAPH"], 2 * h5_parameters["HSLITPOINTS"])
        y_m = np.linspace(-0.5 * h5_parameters["GAPV"], 0.5 * h5_parameters["GAPV"], 2 * h5_parameters["VSLITPOINTS"])
        XX, YY  = np.meshgrid(x_m, y_m)
        theta_x = XX / h5_parameters["DISTANCE"]
        theta_y = YY / h5_parameters["DISTANCE"]

        t0 = time.time()
        total_pd, indiv = power_density_all_harmonics(
            Ky, N_periods, period_m, gamma, current_A,
            theta_x, theta_y, n_harmonics=n_harmonics, p_max=20)
        print(f"Done in {time.time()-t0:.2f} s")

    if calculate_fortran:
        U_X, U_Y, U_POWER_DENSITY, U_POWER_DENSITY_HARMONICS, U_ENERGY_HARMONICS, U_FLUX = (
            run_urgent(h5_parameters, n_harmonics=n_harmonics, do_plot=0, show=0))




    #
    # plots
    #


    if calculate_python:
        #
        # PYTHON CALCULATION
        #
        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        fig.suptitle(f"URGENT-PYTHON (Walker 1992) — K={Ky:.3f}, N={N_periods}, "
                     f"E={gamma*E0_EV/1e9:.1f} GeV, I={current_A*1e3:.0f} mA")
        panels = [('Total H1-%d' % n_harmonics, total_pd)] + \
                 [(f'H{n}', indiv[n]) for n in [1, 2, 3, 4, 5, 6, 7]]
        extent = [x_m[0]*1e3, x_m[-1]*1e3, y_m[0]*1e3, y_m[-1]*1e3]
        for ax, (title, data) in zip(axes.flat, panels):
            im = ax.imshow(data, extent=extent, origin='lower',
                           aspect='equal', cmap='inferno')
            ax.set_title(title, fontsize=9)
            ax.set_xlabel('x [mm]'); ax.set_ylabel('y [mm]')
            plt.colorbar(im, ax=ax, label='W/mrad^2', fraction=0.046)
        plt.tight_layout()

    if calculate_fortran:
        #
        # FORTRAN CODE
        #
        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        fig.suptitle(f"URGENT-FORTRAN (Walker 1992) — K={Ky:.3f}, N={N_periods}, "
                     f"E={gamma*E0_EV/1e9:.1f} GeV, I={current_A*1e3:.0f} mA")
        panels = [('Total H1-%d' % n_harmonics, (U_FLUX.sum(axis=0)).T)] + \
                 [(f'H{n+1}', U_FLUX[n].T) for n in [0, 1, 2, 3, 4, 5, 6]]
        extent = [U_X[0], U_X[-1], U_Y[0], U_Y[-1]]
        for ax, (title, data) in zip(axes.flat, panels):
            im = ax.imshow(data, extent=extent, origin='lower',
                           aspect='equal', cmap='inferno')
            ax.set_title(title, fontsize=9)
            ax.set_xlabel('x [mm]'); ax.set_ylabel('y [mm]')
            plt.colorbar(im, ax=ax, label='W/mrad^2', fraction=0.046)
        plt.tight_layout()

    plt.show()



    #
    # check total intensity
    #
    if 0:
        integ_total = 0.0

        # Numerical integral of your power density map over solid angle
        dtheta_x = theta_x[0, 1] - theta_x[0, 0]  # rad
        dtheta_y = theta_y[1, 0] - theta_y[0, 0]  # rad
        d_omega = dtheta_x * dtheta_y  # sr = rad^2
        d_omega_mrad2 = dtheta_x * dtheta_y * 1e6  # sr = rad^2

        flux_ok = [
        1.06305e+18,
        5.578e+17,
        3.02725e+17,
        1.71251e+17,
        1.0002e+17,
        5.97779e+16,
        3.63149e+16,
        2.23241e+16,
        1.38411e+16,
        8.63484e+15]

        print(f"\n{'n':>4}  {'E [keV]':>10}  {'Peak [????]':>16}  {'Integ [????]':>16} {'Factor ':>16}")
        print("  " + "-"*36)
        for n in range(1, n_harmonics + 1):
            peak = indiv[n].max()
            integ = indiv[n].sum()
            integ_total += integ
            print(f"  {n:2d}    {harmonic_energy(n,Ky,period_m,gamma)/1e3:8.3f}"
                  f"    {peak:14.6g}    {integ:14.6g}  {flux_ok[n-1]/integ:14.6g}")


    # Total power radiated by undulator [W]
    # P = C_gamma * E[GeV]^4 * I[A] * (K^2 * N) / (lambda0 * (1+K^2/2))  -- CHECK
    # Standard formula from Kim 1989:
    # P [kW] = 0.6328 * E[GeV]^2 * B0[T]^2 * N * I[A]  -- wiggler/undulator

    # P_kim
    P_kim = (N_periods / 6) * codata.value('characteristic impedance of vacuum') * \
        current_A * codata.e * 2 * numpy.pi * codata.c * gamma ** 2 * Ky ** 2 / period_m

    print(f"E_GeV     = {E_GeV:.4f} GeV")
    print(f"K         = {Ky:.4f}")
    print(f"Total power (Kim): {P_kim:.2f} W")



"""
harm, power, power-from-flux, power-density-peak, flux: 
0    1110.815    1108.777       1.652  1.06305e+18
1    1154.743    1152.928       1.516  5.578e+17
2     944.851     943.444       1.394  3.02725e+17
3     713.451     712.429       1.023  1.71251e+17
4     518.770     518.051       0.743  1.0002e+17
5     369.225     368.731       0.527  5.97779e+16
6     259.166     258.831       0.380  3.63149e+16
7     180.166     179.942       0.270  2.23241e+16
8     124.332     124.183       0.193  1.38411e+16
9      85.296      85.198       0.139  8.63484e+15

"""



