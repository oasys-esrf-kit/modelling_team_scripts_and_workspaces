"""

Power density calculations and total power on slit after Kim (1989)

Main routines:
calc2d_kim: Power density on a slit from an insertion device (Kim 1989 analytical).
id_power_on_slit: total power [W] through a slit from an insertion device.

Implements Kim (1989) directly:

  Angular power density (Eq. 4.61):
      d2P/dOmega(phi, psi) = (d2P/dOmega)|_0 * fK(gamma*phi, gamma*psi)

  On-axis density, practical units (after Eq. 4.63):
      (d2P/dOmega)|_0 [W/mrad^2] = 10.84 * B0[T] * E^4[GeV] * I[A] * N * G(K)
  (the factor 2N and the interference factor G(K) are already contained here:
   10.84 = 2 * 5.42, and N and G(K) appear explicitly).

  Normalized angular function (Eq. 4.62), valid for ANY K:
      fK(g*phi, g*psi) = 16K/(7*pi*G(K)) * INT_{-pi}^{pi} dxi
                           sin^2(xi)/D^5 * [ (1+X^2-Y^2)^2 + 4 X^2 Y^2 ]
      X = gamma*psi,  Y = K*cos(xi) - gamma*phi,  D = 1 + X^2 + Y^2
  fK is normalized so fK(0,0) = 1.

  Interference factor (Eq. 4.60):
      G(K) = K*(K^6 + (24/7)K^4 + 4K^2 + 16/7) / (1+K^2)^(7/2)   -> 1 for K>>1

Power through the slit is the on-axis density times the integral of fK over
the angular window subtended by the slit (the integral is carried out in
mrad so it matches the W/mrad^2 units of the on-axis density):

      P = (d2P/dOmega)|_0 * INT INT fK(gamma*phi, gamma*psi) d(phi) d(psi)

The slit at distance D subtends:
      phi in [-w_h/(2D), +w_h/(2D)]   horizontal  [rad]
      psi in [-w_v/(2D), +w_v/(2D)]   vertical    [rad]

No bending-magnet approximation and no separate 2N*G(K) multiplication are
needed: the full fK already carries the correct horizontal and vertical angle
dependence, and 2N*G(K) lives in the on-axis density.

References:
    K.-J. Kim, AIP Conf. Proc. 184, 565 (1989). Eqs. 4.60-4.63, 4.66
        https://doi.org/10.1063/1.38046
"""
import numpy
import numpy as np
from scipy.integrate import trapezoid
import scipy.constants as codata
import os
from scipy.ndimage import gaussian_filter

#from id_power_on_slit import fK, G_K, on_axis_density, ELECTRON_REST_GEV



ELECTRON_REST_GEV = codata.value('electron mass energy equivalent in MeV') * 1e-3 # 0.51099895e-3   # m_e c^2 [GeV]


# -- G(K): interference / wiggler factor (Kim 1989 eq. 4.60) ------------------

def G_K(K):
    """
    Interference factor G(K), Kim (1989) eq. (4.60).

        G(K) = K * [K^6 + (24/7)K^4 + 4K^2 + 16/7] / (1 + K^2)^(7/2)

    G(K) -> 1 for K >> 1 (pure wiggler limit); G(K) < 1 at moderate K.
    """
    K2 = K ** 2
    num = K * (K2 ** 3 + (24. / 7.) * K2 ** 2 + 4. * K2 + 16. / 7.)
    den = (1.0 + K2) ** (7. / 2.)
    return num / den


# -- fK: normalized angular power function (Kim 1989 eq. 4.62) ----------------

def fK(gamma_phi, gamma_psi, K, GK=None, n_xi=512):
    """
    Normalized angular power function fK(gamma*phi, gamma*psi), Kim eq. (4.62).

    Parameters
    ----------
    gamma_phi : array_like   gamma * phi  (phi = horizontal angle [rad])
    gamma_psi : array_like   gamma * psi  (psi = vertical   angle [rad])
                             must broadcast against gamma_phi
    K         : float        deflection parameter
    GK        : float        precomputed G(K) (optional)
    n_xi      : int          phase-integration points over xi in [-pi, pi]

    Returns
    -------
    fK : ndarray   dimensionless, with fK(0,0) = 1
    """
    if GK is None:
        GK = G_K(K)

    gphi = np.atleast_1d(np.asarray(gamma_phi, dtype=float))
    gpsi = np.atleast_1d(np.asarray(gamma_psi, dtype=float))
    gphi, gpsi = np.broadcast_arrays(gphi, gpsi)

    xi = np.linspace(-np.pi, np.pi, n_xi)

    # add a trailing xi axis for vectorized integration
    X = gpsi[..., None]                       # X = gamma*psi
    Y = K * np.cos(xi) - gphi[..., None]      # Y = K cos(xi) - gamma*phi
    D = 1.0 + X ** 2 + Y ** 2

    integrand = (np.sin(xi) ** 2 / D ** 5
                 * ((1.0 + X ** 2 - Y ** 2) ** 2 + 4.0 * X ** 2 * Y ** 2))
    val = trapezoid(integrand, xi, axis=-1)

    return (16.0 * K / (7.0 * np.pi * GK)) * val


# -- on-axis power density (Kim 1989 eq. 4.63, practical units) ----------------

def on_axis_density(K, period_m, n_periods, energy_GeV, current, B0_T=None):
    """
    Forward power density (d2P/dOmega)|_0 in W/mrad^2, Kim eq. (4.63):

        (d2P/dOmega)|_0 = 10.84 * B0 * E^4 * I * N * G(K)   [W/mrad^2]
    """
    if B0_T is None:
        # Kim eq. (4.3): K = 0.934 * lambda_u[cm] * B0[T]
        B0_T = K / (0.9337 * period_m * 100.0)
    return 10.84 * B0_T * energy_GeV ** 4 * current * n_periods * G_K(K)


# -- Main functions -------------------------------------------------------------

def id_power_on_slit(
        K=1.6563,
        period_m=0.018,
        n_periods=111,
        energy_GeV=6.0,
        current=0.2,
        slit_h_mm=1.8,
        slit_v_mm=1.0,
        distance=23.0,
        n_phi=101,
        n_psi=101,
):
    """
    Total power [W] through a rectangular slit from an insertion device,
    using Kim (1989) eqs. (4.60), (4.62), (4.63).

    Parameters
    ----------
    K          : float   deflection parameter
    period_m   : float   undulator period [m]
    n_periods  : int     number of periods N
    energy_GeV : float   electron beam energy [GeV]
    current    : float   beam current [A]
    slit_h_mm  : float   slit horizontal full width  [mm]
    slit_v_mm  : float   slit vertical   full height [mm]
    distance   : float   source-to-slit distance [m]
    n_phi      : int     horizontal angle grid points
    n_psi      : int     vertical   angle grid points

    Returns
    -------
    power : float   [W]
    """
    gamma = energy_GeV / ELECTRON_REST_GEV
    GK = G_K(K)

    # on-axis power density [W/mrad^2]  (already contains 2N * G(K))
    dP0 = on_axis_density(K, period_m, n_periods, energy_GeV, current)

    # slit half-angles [rad]
    phi_half = slit_h_mm * 1e-3 / (2.0 * distance)
    psi_half = slit_v_mm * 1e-3 / (2.0 * distance)

    phi = np.linspace(-phi_half, phi_half, n_phi)   # [rad]
    psi = np.linspace(-psi_half, psi_half, n_psi)   # [rad]
    PHI, PSI = np.meshgrid(phi, psi, indexing='ij')

    f = fK(gamma * PHI, gamma * PSI, K, GK=GK)

    # integrate fK over the slit; use mrad to match W/mrad^2
    integral = trapezoid(trapezoid(f, psi * 1e3, axis=1), phi * 1e3)

    return dP0 * integral


"""
calc2d_kim: 2D power density on a slit from an insertion device, computed
analytically with the Kim (1989) formulas, as a drop-in analog of
calc2d_urgent().

The angular power density is (Kim 1989 eqs. 4.61-4.63):

    d2P/dOmega(theta, psi) = (d2P/dOmega)|_0 * fK(gamma*theta, gamma*psi)
    (d2P/dOmega)|_0 = 10.84 * B0[T] * E^4[GeV] * I[A] * N * G(K)  [W/mrad^2]

projected onto the slit plane at distance D:

    x[mm] = theta[mrad]*D[m],  y[mm] = psi[mrad]*D[m]
    dP/dA[W/mm^2] = (d2P/dOmega)[W/mrad^2] / D^2[m^2]

Finite electron-beam emittance (zero_emittance=False) is included by convolving
the single-electron density with a Gaussian whose rms widths at the slit are

    sigma_x = sqrt(sizeH^2 + (divH*D)^2),   sigma_y = sqrt(sizeV^2 + (divV*D)^2)

which broadens the density but conserves the total power (Kim 1989 sec. 4.2.4).

Requires fK, G_K, on_axis_density from id_power_on_slit.py.
"""

scanCounter = 0
def calc2d_kim(bl, zero_emittance=False, fileName=None, fileAppend=False,
               hSlitPoints=21, vSlitPoints=51, n_xi=256):
    r"""
        Power density on a slit from an insertion device (Kim 1989 analytical).

        input:  a dictionary 'bl' with the beamline (same keys as calc2d_urgent)
        output: (hhh, vvv, int_mesh2)
                hhh [mm] horizontal coordinates, length 2*hSlitPoints-1
                vvv [mm] vertical   coordinates, length 2*vSlitPoints-1
                int_mesh2 [W/mm^2]  power density map
                (and, if fileName is given, a SPEC file is written/appended)
    """
    global scanCounter
    print("Inside calc2d_kim")

    # ---- parameters --------------------------------------------------------
    E_GeV   = bl['ElectronEnergy']
    I_A     = bl['ElectronCurrent']
    lam_u   = bl['PeriodID']          # m
    N       = bl['NPeriods']
    K       = bl['Kv']                # planar undulator: vertical field -> use Kv
    Kh      = bl.get('Kh', 0.0)
    D       = bl['distance']          # m
    gapH    = bl['gapH']              # m, full horizontal aperture
    gapV    = bl['gapV']              # m, full vertical   aperture

    if Kh != 0.0:
        print("  WARNING: Kh != 0; calc2d_kim assumes a planar device and uses K = Kv.")

    gamma = E_GeV / ELECTRON_REST_GEV
    B0    = K / (0.9337 * lam_u * 100.0)            # peak field [T], Kim eq. 4.3
    GK    = G_K(K)
    dP0   = on_axis_density(K, lam_u, N, E_GeV, I_A)  # W/mrad^2 (on-axis)
    dA0   = dP0 / D ** 2                              # W/mm^2  (on-axis, at slit)

    # ---- slit grid (mirror a quadrant, exactly like calc2d_urgent) ---------
    hh = numpy.linspace(0.0, gapH * 1e3 / 2.0, hSlitPoints)   # mm (quadrant)
    vv = numpy.linspace(0.0, gapV * 1e3 / 2.0, vSlitPoints)   # mm (quadrant)
    hhh = numpy.concatenate((-hh[::-1], hh[1:]))              # full, 2*hSlitPoints-1
    vvv = numpy.concatenate((-vv[::-1], vv[1:]))              # full, 2*vSlitPoints-1
    dh = hh[1] - hh[0]
    dv = vv[1] - vv[0]

    # convolution widths at the slit plane [mm]
    sigx = numpy.sqrt(bl['ElectronBeamSizeH'] ** 2 +
                      (bl['ElectronBeamDivergenceH'] * D) ** 2) * 1e3
    sigy = numpy.sqrt(bl['ElectronBeamSizeV'] ** 2 +
                      (bl['ElectronBeamDivergenceV'] * D) ** 2) * 1e3

    # ---- evaluate the density ----------------------------------------------
    def density_on(Hcoords, Vcoords):
        """W/mm^2 on the (Hcoords x Vcoords) grid, single electron."""
        gth = gamma * (Hcoords * 1e-3) / D            # gamma*theta (horizontal)
        gps = gamma * (Vcoords * 1e-3) / D            # gamma*psi   (vertical)
        GTH, GPS = numpy.meshgrid(gth, gps, indexing='ij')
        return dA0 * fK(GTH, GPS, K, GK=GK, n_xi=n_xi)

    if zero_emittance:
        int_mesh2 = density_on(hhh, vvv)
    else:
        # pad by 4 sigma so the convolution does not lose the tails, then crop
        npx = int(numpy.ceil(4.0 * sigx / dh))
        npy = int(numpy.ceil(4.0 * sigy / dv))
        Hext = dh * numpy.arange(-(hSlitPoints - 1 + npx), hSlitPoints + npx)
        Vext = dv * numpy.arange(-(vSlitPoints - 1 + npy), vSlitPoints + npy)
        dens_ext = density_on(Hext, Vext)
        dens_ext = gaussian_filter(dens_ext, sigma=(sigx / dh, sigy / dv),
                                   mode='constant', cval=0.0)
        int_mesh2 = dens_ext[npx:npx + (2 * hSlitPoints - 1),
                             npy:npy + (2 * vSlitPoints - 1)]

    totPower = float(int_mesh2.sum() * dh * dv)      # W, power inside the slit

    # ---- write SPEC file (same layout as calc2d_urgent) --------------------
    if fileName is not None:
        if fileAppend:
            f = open(fileName, "a")
        else:
            scanCounter = 0
            f = open(fileName, "w")
            f.write("#F " + fileName + "\n")

        def _header(title, ncol, label, total=None):
            global scanCounter
            scanCounter += 1
            f.write("\n#S %d %s\n" % (scanCounter, title))
            for i, j in bl.items():
                f.write("#UD %s = %s\n" % (i, j))
            f.write("#UD hSlitPoints =  %f\n" % (hSlitPoints))
            f.write("#UD vSlitPoints =  %f\n" % (vSlitPoints))
            f.write("#UD B0 [T] =  %f\n" % (B0))
            f.write("#UD G(K) =  %f\n" % (GK))
            f.write("#UD zero_emittance =  %s\n" % (zero_emittance))
            if total is not None:
                f.write("#UD Total power [W]: " + repr(total) + "\n")
            f.write("#N %d\n" % ncol)
            f.write("#L  " + label + "\n")

        # whole slit (2D map)
        _header("Undulator power density calculation using Kim (whole slit)",
                3, "H[mm]  V[mm]  PowerDensity[W/mm^2]", total=totPower)
        for i in range(len(hhh)):
            for j in range(len(vvv)):
                f.write("%f  %f  %f\n" % (hhh[i], vvv[j], int_mesh2[i, j]))

        # H profile (vertical centre)
        _header("Undulator power density calculation using Kim: H profile",
                2, "H[mm]  PowerDensity[W/mm2]", total=totPower)
        jc = len(vvv) // 2
        for i in range(len(hhh)):
            f.write("%f  %f\n" % (hhh[i], int_mesh2[i, jc]))

        # V profile (horizontal centre)
        _header("Undulator power density calculation using Kim: V profile",
                2, "V[mm]  PowerDensity[W/mm2]", total=totPower)
        ic = len(hhh) // 2
        for j in range(len(vvv)):
            f.write("%f  %f\n" % (vvv[j], int_mesh2[ic, j]))

        f.close()
        where = os.path.join(os.getcwd(), fileName)
        print("Data appended to file: %s" % where if fileAppend
              else "File written to disk: %s" % where)

    print("Power density peak KIM: [W/mm2]: " + repr(float(int_mesh2.max())))
    print("Total power on slit KIM [W]: " + repr(float(totPower)))
    print("\n--------------------------------------------------------\n\n")

    return (hhh, vvv, int_mesh2)


# ── self-test ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    if 0:
        # (1) normalization: fK(0,0) must equal 1
        print("fK normalization check:")
        for K_val in [0.5, 1.0, 1.6563, 3.0]:
            print(f"  K={K_val:6.4f}:  fK(0,0) = {fK(0.0, 0.0, K_val)[0]:.6f}")


    if 0: # power in an opened slit
        # (2) total power: integrating fK over all angles must reproduce the total
        #     wiggler power. NOTE: Kim's *symbolic* eq. (4.66) is correct, but its
        #     *practical* coefficient as printed ("6.4") is ~10x too large; the
        #     internally consistent / standard literature value is 0.633:
        #         P[kW] = 0.633 * E^2[GeV] * B0^2[T] * I[A] * L[m]
        #     (verified: on-axis density eq.4.63 + the fK integral reproduce 0.633,
        #      and symbolic eq.4.66 evaluates to the same total power.)
        print("\nTotal-power check (full-angle integral vs corrected eq. 4.66):")
        print(f"  {'K':>7} {'P_integral[kW]':>15} {'P_eq4.66[kW]':>14} {'ratio':>7}")
        E_GeV, I_A = 6.0, 0.2
        lam, N = 0.018, 111
        for K_val in [0.5, 1.0, 1.6563, 3.0]:
            gamma = E_GeV / ELECTRON_REST_GEV
            B0 = K_val / (0.9337 * lam * 100.0)
            dP0 = on_axis_density(K_val, lam, N, E_GeV, I_A)
            # generous angular window: horizontal ~ K/gamma, vertical ~ 1/gamma
            phi = np.linspace(-(K_val + 6.0) / gamma, (K_val + 6.0) / gamma, 241)
            psi = np.linspace(-8.0 / gamma, 8.0 / gamma, 161)
            PHI, PSI = np.meshgrid(phi, psi, indexing='ij')
            f = fK(gamma * PHI, gamma * PSI, K_val, n_xi=256)
            integ = trapezoid(trapezoid(f, psi * 1e3, axis=1), phi * 1e3)
            P_int = dP0 * integ / 1e3  # kW
            P_ref = 0.633 * E_GeV ** 2 * B0 ** 2 * I_A * (N * lam)  # kW
            print(f"  {K_val:>7.4f} {P_int:>15.3f} {P_ref:>14.3f} {P_int / P_ref:>7.3f}")

    if 0: # power in a 1.8 x 1.0 mm2 slit at 23 m
        cases = [
            ("CPMU18", 1.6563, 0.018, 111, 1000),
            ("CPMU20", 2.334, 0.0205, 98, 1080),
            ("IVU22", 1.543, 0.022, 91, 620),
        ]
        print("\nPower through 1.8 x 1.0 mm slit at 23 m, 6 GeV, 0.2 A:")
        print(f"  {'Undulator':<10} {'K':>7} {'P_slit[W]':>11} {'SRW ref[W]':>11} {'ratio':>7}")
        print("  " + "-" * 50)
        for name, K_val, lam_, N_, srw_ref in cases:
            P = id_power_on_slit(K=K_val, period_m=lam_, n_periods=N_,
                                     energy_GeV=6.0, current=0.2,
                                     slit_h_mm=1.8, slit_v_mm=1.0, distance=23.0)
            print(f"  {name:<10} {K_val:>7.4f} {P:>11.1f} {srw_ref:>11.1f} {P / srw_ref:>7.3f}")

    #
    #
    #
    if 0: # K-scan. importing for external SRW calculation
        # import sys
        # sys.path.append('/modelling_team_scripts_and_workspaces/id11/WATTDOG/SPECTRA/')
        from wattdog_xoppylib import get_id_spectrum

        method = 0  # SRW

        KK = []
        PP=[]
        SS=[]
        for K_val in [0.5, 1.0, 1.2, 1.4, 1.6563]:
            KK.append(K_val)
            P = id_power_on_slit(K=K_val, period_m=0.018, n_periods=111,
                                      energy_GeV=6.0, current=0.2,
                                      slit_h_mm=1.8, slit_v_mm=1.0, distance=23.0)
            PP.append(P)
            s0 = get_id_spectrum(K_val, u='u18', method=method, slit_h_mm=1.8, slit_v_mm=1.0)
            SS.append(s0[3][-1])

        print()
        print("CPMU18 K-scan:")
        print(f"  {'K':>6}  {'P_slit [W]':>12}  {'P_numeric [W]':>12}")
        for i in range(len(KK)):
            print(f"  {KK[i]:>6.4f}  {PP[i]:>12.1f}  {SS[i]:>12.1f}")


    if 0: # scan of slit aperture
        import numpy
        slit_direction = 1 # 0=H, 1=V
        slit_h_mm = 1.8
        slit_v_mm = 1.0

        if slit_direction == 0:
            slit_scan = numpy.linspace(0.001, slit_h_mm, num=11)
        else:
            slit_scan = numpy.linspace(0.001, slit_v_mm, num=11)

        PP = []
        for i, slit_value in enumerate(slit_scan):
            print(i, "from", slit_scan.size)
            if slit_direction == 0:
                P = id_power_on_slit(K=1.6563, period_m=0.018, n_periods=111,
                                      energy_GeV=6.0, current=0.2,
                                      slit_h_mm=slit_value, slit_v_mm=slit_v_mm, distance=23.0)
            else:
                P = id_power_on_slit(K=1.6563, period_m=0.018, n_periods=111,
                                      energy_GeV=6.0, current=0.2,
                                      slit_h_mm=slit_h_mm, slit_v_mm=slit_value, distance=23.0)

            PP.append(P)


        from srxraylib.plot.gol import plot

        plot(slit_scan, PP, xtitle="slit scan", ytitle="power [W]", title="direction = " + "%s"%(['H','V'][slit_direction]))

    if 1: # power density map
        bl = dict(
            ElectronEnergy          = 6.0,
            ElectronCurrent         = 0.2,
            ElectronBeamSizeH       = 3.34281e-05,
            ElectronBeamSizeV       = 7.28139e-06,
            ElectronBeamDivergenceH = 4.51097e-06,
            ElectronBeamDivergenceV = 1.94034e-06,
            ElectronEnergySpread    = 0.001,
            PeriodID                = 0.018,
            NPeriods                = 111,
            Kv                      = 1.6563,
            distance                = 23.0,
            gapH                    = 0.0018,
            gapV                    = 0.001,
        )

        print("=== zero emittance ===")
        h0, v0, m0 = calc2d_kim(bl, zero_emittance=True)
        print("=== with emittance ===")
        h1, v1, m1 = calc2d_kim(bl, zero_emittance=False)

        # cross-check the zero-emittance total against get_id_power_on_slit
        from scipy.integrate import trapezoid
        P_ref = id_power_on_slit(K=bl['Kv'], period_m=bl['PeriodID'],
                                     n_periods=bl['NPeriods'],
                                     energy_GeV=bl['ElectronEnergy'],
                                     current=bl['ElectronCurrent'],
                                     slit_h_mm=bl['gapH'] * 1e3,
                                     slit_v_mm=bl['gapV'] * 1e3,
                                     distance=bl['distance'])
        dh, dv = h0[1] - h0[0], v0[1] - v0[0]
        P_rect = m0.sum() * dh * dv                              # URGENT-style sum
        P_trap = trapezoid(trapezoid(m0, v0, axis=1), h0)        # trapezoid on same grid
        print("cross-check (zero emittance, 1x1 mm slit at 30 m):")
        print(f"  calc2d_kim totPower (sum*dh*dv, URGENT rule) = {P_rect:.2f} W")
        print(f"  same grid, trapezoid rule                    = {P_trap:.2f} W")
        print(f"  get_id_power_on_slit (fine trapezoid)        = {P_ref:.2f} W")
        print(f"  -> trapezoid/reference ratio = {P_trap/P_ref:.4f} "
              f"(rule difference, not a model difference)")

        print(h1.shape, v1.shape, m1.shape)
        # example plot
        if True:
            from srxraylib.plot.gol import plot_image
            plot_image(m1, h1, v1,xtitle="H [mm]",ytitle="V [mm]",title="Power density W/mm2")
        #
        # end script