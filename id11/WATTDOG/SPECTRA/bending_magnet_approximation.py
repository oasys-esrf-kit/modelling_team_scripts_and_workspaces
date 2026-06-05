"""
CPMU18 spectrum — Bending Magnet (wiggler) approximation via srfunc.

get_bending_magnet_spectrum() has exactly the same signature and returns
exactly the same four arrays as get_undulator_spectrum():
    energy [eV], flux [ph/s/0.1%bw], spectral_power [W/eV], cumulated_power [W]

The slit-integrated result is the primary (returned) output, mirroring
what xoppy_calc_undulator_spectrum produces with GAPH/GAPV set.

The full-vertical result (hdiv = 2K/gamma) is computed internally and
also returned via a second helper for comparison plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from srxraylib.sources.srfunc import sync_ene

# -----------------------------------------------------------------------
# Hard-coded machine / energy-range parameters
# -----------------------------------------------------------------------
_E_GEV       = 6.0       # electron energy [GeV]
_I_A         = 0.2       # beam current [A]
_DISTANCE    = 23.0      # source-to-slit distance [m]
_E_MIN_EV    = 1000.0
_E_MAX_EV    = 200000.0
_N_POINTS    = 2000
_PSI_NPOINTS = 50        # vertical integration points (slit case)

_UNDULATOR_DB = {
    'u18': dict(period_m=0.018,  n_periods=111),
    'u20': dict(period_m=0.0205, n_periods=98),
    'u22': dict(period_m=0.022,  n_periods=91),
}

# -----------------------------------------------------------------------
def _bm_flux_and_power(ec_ev, n_periods, f_psi,
                       hdiv_mrad, psi_min=0.0, psi_max=0.0):
    """Internal helper: compute flux, spectral_power, cumulated_power for one case."""
    energy = np.linspace(_E_MIN_EV, _E_MAX_EV, _N_POINTS)

    flux = np.array(sync_ene(
        f_psi       = f_psi,
        energy_ev   = energy,
        ec_ev       = ec_ev,
        e_gev       = _E_GEV,
        i_a         = _I_A,
        hdiv_mrad   = hdiv_mrad,
        psi_min     = psi_min,
        psi_max     = psi_max,
        psi_npoints = _PSI_NPOINTS,
    )).flatten() * 2 * n_periods

    eV_to_J        = 1.60218e-19
    spectral_power = flux * eV_to_J / 1e-3           # W/eV
    de             = energy[1] - energy[0]
    cumulated_power = np.cumsum(spectral_power) * de  # W

    return energy, flux, spectral_power, cumulated_power


# -----------------------------------------------------------------------
def get_bending_magnet_spectrum(K=1.0, u='u18', slit_h_mm=1.8, slit_v_mm=1.0):
    """
    Approximate the undulator spectrum using the 2N × BM (wiggler) formula.

    Identical signature and return arrays to get_undulator_spectrum().
    Returns the slit-integrated spectrum (through GAPH × GAPV aperture).

    Parameters
    ----------
    K         : deflection parameter
    u         : undulator tag ('u18', 'u20', 'u22')
    slit_h_mm : slit horizontal full width [mm]
    slit_v_mm : slit vertical full height [mm]  (TOTAL angle)

    Returns
    -------
    energy          : photon energy [eV]               shape (N,)
    flux            : flux [ph/s/0.1%bw]               shape (N,)
    spectral_power  : spectral power [W/eV]            shape (N,)
    cumulated_power : cumulated power [W]               shape (N,)
    """
    db        = _UNDULATOR_DB[u]
    period_m  = db['period_m']
    n_periods = db['n_periods']

    gamma  = _E_GEV * 1e3 / 0.51099895          # Lorentz factor
    B0_T   = K / (0.9337 * period_m * 100.0)    # peak field [T]
    EC_EV  = 0.6650 * _E_GEV**2 * B0_T * 1e3   # critical energy [eV]

    psi_half_mrad  = (slit_v_mm * 1e-3 / _DISTANCE) * 1e3 / 2.0  # half of total angle
    hdiv_full_mrad = 2.0 * K / gamma * 1e3
    hdiv_slit_mrad = (slit_h_mm * 1e-3 / _DISTANCE) * 1e3
    hdiv_used_mrad = min(hdiv_slit_mrad, hdiv_full_mrad)

    return _bm_flux_and_power(
        ec_ev     = EC_EV,
        n_periods = n_periods,
        f_psi     = 2,
        hdiv_mrad = hdiv_used_mrad,
        psi_min   = -psi_half_mrad,
        psi_max   = +psi_half_mrad,
    )


def get_bending_magnet_spectrum_full_vertical(K=1.0, u='u18',
                                              slit_h_mm=1.8, slit_v_mm=1.0):
    """
    Same as get_bending_magnet_spectrum() but integrates over the full vertical
    emission with hdiv = 2K/gamma (no slit clipping in V).
    slit_h_mm / slit_v_mm kept in signature for interface compatibility but
    only slit_h_mm is used (horizontal clipping still applies).
    """
    db        = _UNDULATOR_DB[u]
    period_m  = db['period_m']
    n_periods = db['n_periods']

    gamma  = _E_GEV * 1e3 / 0.51099895
    B0_T   = K / (0.9337 * period_m * 100.0)
    EC_EV  = 0.6650 * _E_GEV**2 * B0_T * 1e3

    hdiv_full_mrad = 2.0 * K / gamma * 1e3

    return _bm_flux_and_power(
        ec_ev     = EC_EV,
        n_periods = n_periods,
        f_psi     = 0,
        hdiv_mrad = hdiv_full_mrad,
    )


# -----------------------------------------------------------------------
# Main script
# -----------------------------------------------------------------------
if __name__ == '__main__':

    K         = 1.6563
    u         = 'u18'
    slit_h_mm = 1.8
    slit_v_mm = 1.0

    # --- primary result: slit (same as get_undulator_spectrum) ---
    energy, flux, spectral_power, cumulated_power = get_bending_magnet_spectrum(
        K=K, u=u, slit_h_mm=slit_h_mm, slit_v_mm=slit_v_mm,
    )

    # --- full-vertical result for comparison ---
    _, flux_full, sp_full, cp_full = get_bending_magnet_spectrum_full_vertical(
        K=K, u=u,
    )

    # Labels / derived values for plot annotation
    db        = _UNDULATOR_DB[u]
    period_m  = db['period_m']
    n_periods = db['n_periods']
    gamma     = _E_GEV * 1e3 / 0.51099895
    B0_T      = K / (0.9337 * period_m * 100.0)
    EC_EV     = 0.6650 * _E_GEV**2 * B0_T * 1e3
    L_m       = n_periods * period_m
    hdiv_full = 2.0 * K / gamma * 1e3
    psi_half  = (slit_v_mm * 1e-3 / _DISTANCE) * 1e3 / 2.0
    hdiv_used = min((slit_h_mm * 1e-3 / _DISTANCE) * 1e3, hdiv_full)
    P_total_W = 0.6328e3 * _E_GEV**2 * B0_T**2 * L_m * _I_A

    print(f"  B0={B0_T:.4f} T  Ec={EC_EV/1e3:.3f} keV  2K/gamma={hdiv_full:.4f} mrad")
    print(f"  Power full V : {cp_full[-1]:.1f} W")
    print(f"  Power slit   : {cumulated_power[-1]:.1f} W")
    print(f"  P_total anal : {P_total_W:.1f} W")

    label_full = f"Full V, hdiv=2K/γ={hdiv_full:.3f} mrad  →  {cp_full[-1]:.0f} W"
    label_slit = (f"Slit {slit_h_mm}×{slit_v_mm} mm² @{_DISTANCE} m  "
                  f"(hdiv={hdiv_used:.3f} mrad, ψ=±{psi_half:.3f} mrad)"
                  f"  →  {cumulated_power[-1]:.0f} W")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        f"{u.upper()} BM approximation  |  "
        f"E={_E_GEV} GeV, K={K}, B0={B0_T:.3f} T, "
        f"Ec={EC_EV/1e3:.1f} keV, N={n_periods}, I={_I_A} A",
        fontsize=10,
    )

    ax = axes[0]
    ax.semilogy(energy / 1e3, flux_full,       label=label_full)
    ax.semilogy(energy / 1e3, flux,       '--', label=label_slit)
    ax.set_xlabel("Photon energy [keV]")
    ax.set_ylabel("Flux [ph/s/0.1%bw]")
    ax.set_title("Flux  (2N × BM)")
    ax.legend(fontsize=7)
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1]
    ax.plot(energy / 1e3, sp_full * 1e3,              label=label_full)
    ax.plot(energy / 1e3, spectral_power * 1e3,  '--', label=label_slit)
    ax.set_xlabel("Photon energy [keV]")
    ax.set_ylabel("Spectral power [mW/eV]")
    ax.set_title("Spectral power  (2N × BM)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.plot(energy / 1e3, cp_full,                  label=label_full)
    ax.plot(energy / 1e3, cumulated_power,      '--', label=label_slit)
    ax.axhline(P_total_W, color='gray', ls=':', lw=1.5,
               label=f"Analytical total  {P_total_W:.0f} W")
    ax.set_xlabel("Photon energy [keV]")
    ax.set_ylabel("Cumulated power [W]")
    ax.set_title("Cumulated power  (2N × BM)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    # outpng = "/mnt/user-data/outputs/cpmu18_bm_approximation.png"
    # plt.savefig(outpng, dpi=150, bbox_inches='tight')
    plt.show()
    print(f"\nPlot saved to {outpng}")