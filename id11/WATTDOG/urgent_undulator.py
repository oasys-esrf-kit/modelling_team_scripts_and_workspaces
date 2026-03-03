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
    prefactor = gamma ** 4 * current_A / (period_m * N_periods)

    return prefactor * An2 / denom * 1e-6   # W/mrad^2


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



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import time

    n_harmonics = 16 # 24
    period_m  = 0.018
    Ky        = 1.0 # 1.358
    B0 = Ky / (0.9336 * period_m * 1e2)
    N_periods = 111
    gamma     = 11741.707085460033
    current_A = 0.2
    Z_m       = 30.0

    print(f"K  = {Ky:.4f}")
    print(f"E1 = {harmonic_energy(1, Ky, period_m, gamma)/1e3:.3f} keV")
    print(f"B0 = {B0:.3f} T")

    # Total power radiated by undulator [W]
    # P = C_gamma * E[GeV]^4 * I[A] * (K^2 * N) / (lambda0 * (1+K^2/2))  -- CHECK
    # Standard formula from Kim 1989:
    # P [kW] = 0.6328 * E[GeV]^2 * B0[T]^2 * N * I[A]  -- wiggler/undulator
    E_GeV = 6.0
    # P_kim
    P_kim = (N_periods / 6) * codata.value('characteristic impedance of vacuum') * \
        current_A * codata.e * 2 * numpy.pi * codata.c * gamma ** 2 * Ky ** 2 / period_m

    print(f"E_GeV     = {E_GeV:.4f} GeV")
    print(f"K         = {Ky:.4f}")
    print(f"Total power: {P_kim:.2f} W")



    x_m = np.linspace(-0.005, 0.005, 150)
    y_m = np.linspace(-0.005, 0.005, 150)
    XX, YY  = np.meshgrid(x_m, y_m)
    theta_x = XX / Z_m
    theta_y = YY / Z_m

    t0 = time.time()
    total_pd, indiv = power_density_all_harmonics(
        Ky, N_periods, period_m, gamma, current_A,
        theta_x, theta_y, n_harmonics=n_harmonics, p_max=20)
    print(f"Done in {time.time()-t0:.2f} s")
    print(f"Peak total: {total_pd.max():.4f} W/mrad^2")

    print(f"\n{'n':>4}  {'E [keV]':>10}  {'Peak [W/mrad^2]':>16}  {'Peak [W/mm^2]':>16} {'Integ [??]':>16}")
    print("  " + "-"*36)
    integ_total = 0.0

    # Numerical integral of your power density map over solid angle
    dtheta_x = theta_x[0, 1] - theta_x[0, 0]  # rad
    dtheta_y = theta_y[1, 0] - theta_y[0, 0]  # rad
    d_omega = dtheta_x * dtheta_y  # sr = rad^2
    d_omega_mrad2 = dtheta_x * dtheta_y * 1e6  # sr = rad^2


    for n in range(1, n_harmonics + 1):
        peak = indiv[n].max() / (Z_m ** 2)
        integ = indiv[n].sum()
        integ_total += integ
        print(f"  {n:2d}    {harmonic_energy(n,Ky,period_m,gamma)/1e3:8.3f}"
              f"    {indiv[n].max():14.6g}    {peak:14.6g}    {integ:14.6g}")

    P_numerical = total_pd.sum() * d_omega_mrad2  # W

    ratio = 0.00023356698217765035 # P_kim / P_numerical

    P_corrected = P_numerical * ratio
    print("P_numerical, ratio, P_corrected = : ", P_numerical, ratio, P_corrected)


    # plot

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    fig.suptitle(f"URGENT (Walker 1992) — K={Ky:.3f}, N={N_periods}, "
                 f"E={gamma*E0_EV/1e9:.1f} GeV, I={current_A*1e3:.0f} mA")
    panels = [('Total H1-200', total_pd)] + \
             [(f'H{n}', indiv[n]) for n in [1, 2, 3, 4, 5, 6, 7]]
    extent = [x_m[0]*1e3, x_m[-1]*1e3, y_m[0]*1e3, y_m[-1]*1e3]
    for ax, (title, data) in zip(axes.flat, panels):
        im = ax.imshow(data, extent=extent, origin='lower',
                       aspect='equal', cmap='inferno')
        ax.set_title(title, fontsize=9)
        ax.set_xlabel('x [mm]'); ax.set_ylabel('y [mm]')
        plt.colorbar(im, ax=ax, label='W/mrad^2', fraction=0.046)
    plt.tight_layout()
    # plt.savefig("power_density_urgent.png", dpi=150, bbox_inches='tight')
    plt.show()




