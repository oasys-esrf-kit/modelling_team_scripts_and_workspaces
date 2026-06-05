import numpy
import scipy.constants as codata
from fontTools.ttLib import xmlToTag


def K_vs_gap(gap_mm=6.0, u='u18', ):
    Bmax = 0.0
    if u == 'u18':
        period_length = 0.018  # m
        A = [2.978, 0.0, 0.988, 1.0826, 1.0, 1.1456]
        i_half = len(A) // 2
        for i in range(i_half): Bmax += A[i] * numpy.exp(-numpy.pi * (i + 1) * A[i + i_half] * gap_mm / (period_length * 1e3))
    elif u == 'u20':
        period_length = 0.0205  # m
        A = [3.2379, 0.0, 1.005, 1.1008, 1.0, 1.1468]
        i_half = len(A) // 2
        for i in range(i_half): Bmax += A[i] * numpy.exp(-numpy.pi * (i + 1) * A[i + i_half] * gap_mm / (period_length * 1e3))
    elif u == 'u22':
        period_length = 0.022  # m
        A = [2.2342, 1.1228]
        i_half = len(A) // 2
        for i in range(i_half): Bmax += A[i] * numpy.exp(-numpy.pi * (i + 1) * A[i + i_half] * gap_mm / (period_length * 1e3))

    return Bmax * period_length * codata.e / (2 * numpy.pi * codata.m_e * codata.c)


def get_id_spectrum( K=1.0, u='u18', method=0, slit_h_mm=1.8, slit_v_mm=1.0):
    if method == 0:
        return get_undulator_spectrum(K, u, slit_h_mm=slit_h_mm, slit_v_mm=slit_v_mm)
    elif method == 1:
        return get_ws_spectrum(K, u, slit_h_mm=slit_h_mm, slit_v_mm=slit_v_mm)
    elif method == 2:
        return get_undulator_power_density_from_harmonics(K, u, slit_h_mm=slit_h_mm, slit_v_mm=slit_v_mm)
    elif method == 3:
        return get_bending_magnet_spectrum(K, u, slit_h_mm=slit_h_mm, slit_v_mm=slit_v_mm)

def get_undulator_spectrum( K=1.0, u='u18', slit_h_mm=1.8, slit_v_mm=1.0):

    if u == 'u18':
        #
        # script to make the calculations (created by XOPPY:undulator_spectrum)
        #
        from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_spectrum
        energy, flux, spectral_power, cumulated_power = xoppy_calc_undulator_spectrum(
            ELECTRONENERGY=6.0,
            ELECTRONENERGYSPREAD=0.001,
            ELECTRONCURRENT=0.2,
            ELECTRONBEAMSIZEH=3.34281e-05,
            ELECTRONBEAMSIZEV=7.28139e-06,
            ELECTRONBEAMDIVERGENCEH=4.51097e-06,
            ELECTRONBEAMDIVERGENCEV=1.94034e-06,
            PERIODID=0.018,
            NPERIODS=111,
            KV=K,
            KH=0.0,
            KPHASE=0.0,
            DISTANCE=23.0,
            GAPH=1e-3 * slit_h_mm,
            GAPV=1e-3 * slit_v_mm,
            GAPH_CENTER=0.0,
            GAPV_CENTER=0.0,
            PHOTONENERGYMIN=1000.0,
            PHOTONENERGYMAX=200000.0,
            PHOTONENERGYPOINTS=2000,
            METHOD=2,
            USEEMITTANCES=1)
    elif u == 'u20':
        #
        # script to make the calculations (created by XOPPY:undulator_spectrum)
        #
        from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_spectrum
        energy, flux, spectral_power, cumulated_power = xoppy_calc_undulator_spectrum(
            ELECTRONENERGY=6.0,
            ELECTRONENERGYSPREAD=0.001,
            ELECTRONCURRENT=0.2,
            ELECTRONBEAMSIZEH=3.34281e-05,
            ELECTRONBEAMSIZEV=7.28139e-06,
            ELECTRONBEAMDIVERGENCEH=4.51097e-06,
            ELECTRONBEAMDIVERGENCEV=1.94034e-06,
            PERIODID=0.0205,
            NPERIODS=98,
            KV=K,
            KH=0.0,
            KPHASE=0.0,
            DISTANCE=23.0,
            GAPH=1e-3 * slit_h_mm,
            GAPV=1e-3 * slit_v_mm,
            GAPH_CENTER=0.0,
            GAPV_CENTER=0.0,
            PHOTONENERGYMIN=1000.0,
            PHOTONENERGYMAX=200000.0,
            PHOTONENERGYPOINTS=2000,
            METHOD=2,
            USEEMITTANCES=1)
    elif u == 'u22':
        #
        # script to make the calculations (created by XOPPY:undulator_spectrum)
        #
        from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_spectrum
        energy, flux, spectral_power, cumulated_power = xoppy_calc_undulator_spectrum(
            ELECTRONENERGY=6.0,
            ELECTRONENERGYSPREAD=0.001,
            ELECTRONCURRENT=0.2,
            ELECTRONBEAMSIZEH=3.34281e-05,
            ELECTRONBEAMSIZEV=7.28139e-06,
            ELECTRONBEAMDIVERGENCEH=4.51097e-06,
            ELECTRONBEAMDIVERGENCEV=1.94034e-06,
            PERIODID=0.022,
            NPERIODS=91,
            KV=K,
            KH=0.0,
            KPHASE=0.0,
            DISTANCE=23.0,
            GAPH=1e-3 * slit_h_mm,
            GAPV=1e-3 * slit_v_mm,
            GAPH_CENTER=0.0,
            GAPV_CENTER=0.0,
            PHOTONENERGYMIN=1000.0,
            PHOTONENERGYMAX=200000.0,
            PHOTONENERGYPOINTS=2000,
            METHOD=2,
            USEEMITTANCES=1)

    return energy, flux, spectral_power, cumulated_power

def get_ws_spectrum( K=1.0, u='u18', slit_h_mm=1.8, slit_v_mm=1.0):

    if u == 'u18':
        #
        # script to make the calculations (created by XOPPY:WS)
        #
        import numpy
        from xoppylib.xoppy_run_binaries import xoppy_calc_ws

        out_file = xoppy_calc_ws(
            ENERGY=6.0,
            CUR=200.0,
            PERIOD=1.7999999999999998,
            N=111.0,
            KX=0.0,
            KY=K,
            EMIN=1000.0,
            EMAX=200000.0,
            NEE=2000,
            D=23.0,
            XPC=0.0,
            YPC=0.0,
            XPS=slit_h_mm,
            YPS=slit_v_mm,
            NXP=10,
            NYP=10,
        )

        # data to pass to power
        data = numpy.loadtxt(out_file)
        energy = data[:, 0]
        flux = data[:, 1]
        spectral_power = data[:, 2]
        cumulated_power = data[:, 3]

    elif u == 'u20':

        #
        # script to make the calculations (created by XOPPY:WS)
        #
        import numpy
        from xoppylib.xoppy_run_binaries import xoppy_calc_ws

        out_file = xoppy_calc_ws(
            ENERGY=6.0,
            CUR=200.0,
            PERIOD=2.0500000000000003,
            N=98.0,
            KX=0.0,
            KY=K,
            EMIN=1000.0,
            EMAX=200000.0,
            NEE=2000,
            D=23.0,
            XPC=0.0,
            YPC=0.0,
            XPS=slit_h_mm,
            YPS=slit_v_mm,
            NXP=10,
            NYP=10,
        )

        # data to pass to power
        data = numpy.loadtxt(out_file)
        energy = data[:, 0]
        flux = data[:, 1]
        spectral_power = data[:, 2]
        cumulated_power = data[:, 3]

    elif u == 'u22':
        #
        # script to make the calculations (created by XOPPY:WS)
        #
        import numpy
        from xoppylib.xoppy_run_binaries import xoppy_calc_ws

        out_file = xoppy_calc_ws(
            ENERGY=6.0,
            CUR=200.0,
            PERIOD=2.1999999999999997,
            N=91.0,
            KX=0.0,
            KY=K,
            EMIN=1000.0,
            EMAX=200000.0,
            NEE=2000,
            D=23.0,
            XPC=0.0,
            YPC=0.0,
            XPS=slit_h_mm,
            YPS=slit_v_mm,
            NXP=10,
            NYP=10,
        )

        # data to pass to power
        data = numpy.loadtxt(out_file)
        energy = data[:, 0]
        flux = data[:, 1]
        spectral_power = data[:, 2]
        cumulated_power = data[:, 3]

    return energy, flux, spectral_power, cumulated_power

def get_undulator_power_density_from_harmonics( K=1.0, u='u18', slit_h_mm=1.8, slit_v_mm=1.0, use_python=0):

    if u == 'u18':
        #
        # script to make the calculations (created by XOPPY:undulator_spectrum)
        #
        from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_power_density_from_harmonics

        # define inputs here to be written in h5 file
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
        h5_parameters["KV"] = K
        h5_parameters["KH"] = 0.0
        h5_parameters["KPHASE"] = 0.0
        h5_parameters["DISTANCE"] = 23.0
        h5_parameters["GAPH"] = slit_h_mm * 1e-3
        h5_parameters["GAPV"] = slit_v_mm * 1e-3
        h5_parameters["HSLITPOINTS"] = 41
        h5_parameters["VSLITPOINTS"] = 41
        h5_parameters["METHOD"] = use_python  # 0=urgent (fortran), 1=urgentpy (python)
        h5_parameters["USEEMITTANCES"] = 1
        h5_parameters["MASK_FLAG"] = 1
        h5_parameters["MASK_ROT_H_DEG"] = 0.0
        h5_parameters["MASK_ROT_V_DEG"] = 0.0
        h5_parameters["MASK_H_MIN"] = -0.5 * slit_h_mm
        h5_parameters["MASK_H_MAX"] =  0.5 * slit_h_mm
        h5_parameters["MASK_V_MIN"] = -0.5 * slit_v_mm
        h5_parameters["MASK_V_MAX"] =  0.5 * slit_v_mm
        h5_parameters["harmonic_max"] = 25  # maximum harmonic to calculate
        h5_parameters["photon_energy_bin"] = 500.0  # in eV, for calculating Spectral Power

        h, v, power_density, code, power_density_harmonics, energy_harmonics, spectral_power, spectral_power_energy = xoppy_calc_undulator_power_density_from_harmonics(
            ELECTRONENERGY=h5_parameters["ELECTRONENERGY"],
            ELECTRONENERGYSPREAD=h5_parameters["ELECTRONENERGYSPREAD"],
            ELECTRONCURRENT=h5_parameters["ELECTRONCURRENT"],
            ELECTRONBEAMSIZEH=h5_parameters["ELECTRONBEAMSIZEH"],
            ELECTRONBEAMSIZEV=h5_parameters["ELECTRONBEAMSIZEV"],
            ELECTRONBEAMDIVERGENCEH=h5_parameters["ELECTRONBEAMDIVERGENCEH"],
            ELECTRONBEAMDIVERGENCEV=h5_parameters["ELECTRONBEAMDIVERGENCEV"],
            PERIODID=h5_parameters["PERIODID"],
            NPERIODS=h5_parameters["NPERIODS"],
            KV=h5_parameters["KV"],
            KH=h5_parameters["KH"],
            KPHASE=h5_parameters["KPHASE"],
            DISTANCE=h5_parameters["DISTANCE"],
            GAPH=h5_parameters["GAPH"],
            GAPV=h5_parameters["GAPV"],
            HSLITPOINTS=h5_parameters["HSLITPOINTS"],
            VSLITPOINTS=h5_parameters["VSLITPOINTS"],
            METHOD=h5_parameters["METHOD"],
            harmonic_max=h5_parameters["harmonic_max"],
            USEEMITTANCES=h5_parameters["USEEMITTANCES"],
            MASK_FLAG=h5_parameters["MASK_FLAG"],
            MASK_ROT_H_DEG=h5_parameters["MASK_ROT_H_DEG"],
            MASK_ROT_V_DEG=h5_parameters["MASK_ROT_V_DEG"],
            MASK_H_MIN=h5_parameters["MASK_H_MIN"],
            MASK_H_MAX=h5_parameters["MASK_H_MAX"],
            MASK_V_MIN=h5_parameters["MASK_V_MIN"],
            MASK_V_MAX=h5_parameters["MASK_V_MAX"],
            h5_file="undulator_power_density_from_harmonics.h5",
            h5_entry_name="XOPPY_POWERDENSITY_FROM_HARMONICS",
            h5_initialize=True,
            h5_parameters=h5_parameters,
            photon_energy_bin=h5_parameters["photon_energy_bin"],
        )
    elif u == 'u20':

        #
        # script to make the calculations (created by XOPPY:undulator_spectrum)
        #
        from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_power_density_from_harmonics

        # define inputs here to be written in h5 file
        h5_parameters = dict()
        h5_parameters["ELECTRONENERGY"] = 6.0
        h5_parameters["ELECTRONENERGYSPREAD"] = 0.001
        h5_parameters["ELECTRONCURRENT"] = 0.2
        h5_parameters["ELECTRONBEAMSIZEH"] = 3.34281e-05
        h5_parameters["ELECTRONBEAMSIZEV"] = 7.28139e-06
        h5_parameters["ELECTRONBEAMDIVERGENCEH"] = 4.51097e-06
        h5_parameters["ELECTRONBEAMDIVERGENCEV"] = 1.94034e-06
        h5_parameters["PERIODID"] = 0.0205
        h5_parameters["NPERIODS"] = 98
        h5_parameters["KV"] = K
        h5_parameters["KH"] = 0.0
        h5_parameters["KPHASE"] = 0.0
        h5_parameters["DISTANCE"] = 23.0
        h5_parameters["GAPH"] = slit_h_mm * 1e-3
        h5_parameters["GAPV"] = slit_v_mm * 1e-3
        h5_parameters["HSLITPOINTS"] = 41
        h5_parameters["VSLITPOINTS"] = 41
        h5_parameters["METHOD"] = use_python  # 0=urgent (fortran), 1=urgentpy (python)
        h5_parameters["USEEMITTANCES"] = 1
        h5_parameters["MASK_FLAG"] = 1
        h5_parameters["MASK_ROT_H_DEG"] = 0.0
        h5_parameters["MASK_ROT_V_DEG"] = 0.0
        h5_parameters["MASK_H_MIN"] = -0.5 * slit_h_mm
        h5_parameters["MASK_H_MAX"] =  0.5 * slit_h_mm
        h5_parameters["MASK_V_MIN"] = -0.5 * slit_v_mm
        h5_parameters["MASK_V_MAX"] =  0.5 * slit_v_mm
        h5_parameters["harmonic_max"] = 45  # maximum harmonic to calculate
        h5_parameters["photon_energy_bin"] = 500.0  # in eV, for calculating Spectral Power

        h, v, power_density, code, power_density_harmonics, energy_harmonics, spectral_power, spectral_power_energy = xoppy_calc_undulator_power_density_from_harmonics(
            ELECTRONENERGY=h5_parameters["ELECTRONENERGY"],
            ELECTRONENERGYSPREAD=h5_parameters["ELECTRONENERGYSPREAD"],
            ELECTRONCURRENT=h5_parameters["ELECTRONCURRENT"],
            ELECTRONBEAMSIZEH=h5_parameters["ELECTRONBEAMSIZEH"],
            ELECTRONBEAMSIZEV=h5_parameters["ELECTRONBEAMSIZEV"],
            ELECTRONBEAMDIVERGENCEH=h5_parameters["ELECTRONBEAMDIVERGENCEH"],
            ELECTRONBEAMDIVERGENCEV=h5_parameters["ELECTRONBEAMDIVERGENCEV"],
            PERIODID=h5_parameters["PERIODID"],
            NPERIODS=h5_parameters["NPERIODS"],
            KV=h5_parameters["KV"],
            KH=h5_parameters["KH"],
            KPHASE=h5_parameters["KPHASE"],
            DISTANCE=h5_parameters["DISTANCE"],
            GAPH=h5_parameters["GAPH"],
            GAPV=h5_parameters["GAPV"],
            HSLITPOINTS=h5_parameters["HSLITPOINTS"],
            VSLITPOINTS=h5_parameters["VSLITPOINTS"],
            METHOD=h5_parameters["METHOD"],
            harmonic_max=h5_parameters["harmonic_max"],
            USEEMITTANCES=h5_parameters["USEEMITTANCES"],
            MASK_FLAG=h5_parameters["MASK_FLAG"],
            MASK_ROT_H_DEG=h5_parameters["MASK_ROT_H_DEG"],
            MASK_ROT_V_DEG=h5_parameters["MASK_ROT_V_DEG"],
            MASK_H_MIN=h5_parameters["MASK_H_MIN"],
            MASK_H_MAX=h5_parameters["MASK_H_MAX"],
            MASK_V_MIN=h5_parameters["MASK_V_MIN"],
            MASK_V_MAX=h5_parameters["MASK_V_MAX"],
            h5_file="undulator_power_density_from_harmonics.h5",
            h5_entry_name="XOPPY_POWERDENSITY_FROM_HARMONICS",
            h5_initialize=True,
            h5_parameters=h5_parameters,
            photon_energy_bin=h5_parameters["photon_energy_bin"],
        )
    elif u == 'u22':

        #
        # script to make the calculations (created by XOPPY:undulator_spectrum)
        #
        from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_power_density_from_harmonics

        # define inputs here to be written in h5 file
        h5_parameters = dict()
        h5_parameters["ELECTRONENERGY"] = 6.0
        h5_parameters["ELECTRONENERGYSPREAD"] = 0.001
        h5_parameters["ELECTRONCURRENT"] = 0.2
        h5_parameters["ELECTRONBEAMSIZEH"] = 3.34281e-05
        h5_parameters["ELECTRONBEAMSIZEV"] = 7.28139e-06
        h5_parameters["ELECTRONBEAMDIVERGENCEH"] = 4.51097e-06
        h5_parameters["ELECTRONBEAMDIVERGENCEV"] = 1.94034e-06
        h5_parameters["PERIODID"] = 0.022
        h5_parameters["NPERIODS"] = 91
        h5_parameters["KV"] = K
        h5_parameters["KH"] = 0.0
        h5_parameters["KPHASE"] = 0.0
        h5_parameters["DISTANCE"] = 23.0
        h5_parameters["GAPH"] = slit_h_mm * 1e-3
        h5_parameters["GAPV"] = slit_v_mm * 1e-3
        h5_parameters["HSLITPOINTS"] = 41
        h5_parameters["VSLITPOINTS"] = 41
        h5_parameters["METHOD"] = use_python  # 0=urgent (fortran), 1=urgentpy (python)
        h5_parameters["USEEMITTANCES"] = 1
        h5_parameters["MASK_FLAG"] = 1
        h5_parameters["MASK_ROT_H_DEG"] = 0.0
        h5_parameters["MASK_ROT_V_DEG"] = 0.0
        h5_parameters["MASK_H_MIN"] = -0.5 * slit_h_mm
        h5_parameters["MASK_H_MAX"] =  0.5 * slit_h_mm
        h5_parameters["MASK_V_MIN"] = -0.5 * slit_v_mm
        h5_parameters["MASK_V_MAX"] =  0.5 * slit_v_mm
        h5_parameters["harmonic_max"] = 30  # maximum harmonic to calculate
        h5_parameters["photon_energy_bin"] = 500.0  # in eV, for calculating Spectral Power

        h, v, power_density, code, power_density_harmonics, energy_harmonics, spectral_power, spectral_power_energy = xoppy_calc_undulator_power_density_from_harmonics(
            ELECTRONENERGY=h5_parameters["ELECTRONENERGY"],
            ELECTRONENERGYSPREAD=h5_parameters["ELECTRONENERGYSPREAD"],
            ELECTRONCURRENT=h5_parameters["ELECTRONCURRENT"],
            ELECTRONBEAMSIZEH=h5_parameters["ELECTRONBEAMSIZEH"],
            ELECTRONBEAMSIZEV=h5_parameters["ELECTRONBEAMSIZEV"],
            ELECTRONBEAMDIVERGENCEH=h5_parameters["ELECTRONBEAMDIVERGENCEH"],
            ELECTRONBEAMDIVERGENCEV=h5_parameters["ELECTRONBEAMDIVERGENCEV"],
            PERIODID=h5_parameters["PERIODID"],
            NPERIODS=h5_parameters["NPERIODS"],
            KV=h5_parameters["KV"],
            KH=h5_parameters["KH"],
            KPHASE=h5_parameters["KPHASE"],
            DISTANCE=h5_parameters["DISTANCE"],
            GAPH=h5_parameters["GAPH"],
            GAPV=h5_parameters["GAPV"],
            HSLITPOINTS=h5_parameters["HSLITPOINTS"],
            VSLITPOINTS=h5_parameters["VSLITPOINTS"],
            METHOD=h5_parameters["METHOD"],
            harmonic_max=h5_parameters["harmonic_max"],
            USEEMITTANCES=h5_parameters["USEEMITTANCES"],
            MASK_FLAG=h5_parameters["MASK_FLAG"],
            MASK_ROT_H_DEG=h5_parameters["MASK_ROT_H_DEG"],
            MASK_ROT_V_DEG=h5_parameters["MASK_ROT_V_DEG"],
            MASK_H_MIN=h5_parameters["MASK_H_MIN"],
            MASK_H_MAX=h5_parameters["MASK_H_MAX"],
            MASK_V_MIN=h5_parameters["MASK_V_MIN"],
            MASK_V_MAX=h5_parameters["MASK_V_MAX"],
            h5_file="undulator_power_density_from_harmonics.h5",
            h5_entry_name="XOPPY_POWERDENSITY_FROM_HARMONICS",
            h5_initialize=True,
            h5_parameters=h5_parameters,
            photon_energy_bin=h5_parameters["photon_energy_bin"],
        )

    de = spectral_power_energy[1] - spectral_power_energy[0]
    energy = spectral_power_energy
    flux = energy * 0  # TODO
    # spectral_power = spectral_power
    cumulated_power = (spectral_power * de).cumsum()

    return energy, flux, spectral_power, cumulated_power



###############################################

# import numpy as np


# -----------------------------------------------------------------------
# Hard-coded machine / energy-range parameters
# -----------------------------------------------------------------------
# _E_GEV       = 6.0       # electron energy [GeV]
# _I_A         = 0.2       # beam current [A]
# _DISTANCE    = 23.0      # source-to-slit distance [m]
# _E_MIN_EV    = 1000.0
# _E_MAX_EV    = 200000.0
# _N_POINTS    = 2000
# _PSI_NPOINTS = 50        # vertical integration points (slit case)



# # -----------------------------------------------------------------------
# def _bm_flux_and_power(ec_ev, n_periods, f_psi,
#                        hdiv_mrad, psi_min=0.0, psi_max=0.0):
#     from srxraylib.sources.srfunc import sync_ene
#     """Internal helper: compute flux, spectral_power, cumulated_power for one case."""
#     energy = np.linspace(_E_MIN_EV, _E_MAX_EV, _N_POINTS)
#
#     flux = np.array(sync_ene(
#         f_psi       = f_psi,
#         energy_ev   = energy,
#         ec_ev       = ec_ev,
#         e_gev       = _E_GEV,
#         i_a         = _I_A,
#         hdiv_mrad   = hdiv_mrad,
#         psi_min     = psi_min,
#         psi_max     = psi_max,
#         psi_npoints = _PSI_NPOINTS,
#     )).flatten() * 2 * n_periods
#
#     eV_to_J        = 1.60218e-19
#     spectral_power = flux * eV_to_J / 1e-3           # W/eV
#     de             = energy[1] - energy[0]
#     cumulated_power = np.cumsum(spectral_power) * de  # W
#
#     return energy, flux, spectral_power, cumulated_power



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
    from srxraylib.sources.srfunc import sync_ene

    if u == 'u18':
        period_m=0.018
        n_periods=111
    elif u == 'u20':
        period_m=0.0205
        n_periods=98
    elif u == 'u22':
        period_m=0.022
        n_periods=91


    _E_GEV = 6.0
    _I_A         = 0.2       # beam current [A]
    _DISTANCE    = 23.0      # source-to-slit distance [m]
    _E_MIN_EV    = 1000.0
    _E_MAX_EV    = 200000.0
    _N_POINTS    = 2000
    _PSI_NPOINTS = 50        # vertical integration points (slit case)

    gamma  = _E_GEV * 1e3 / 0.51099895          # Lorentz factor
    B0_T   = K / (0.9337 * period_m * 100.0)    # peak field [T]
    EC_EV  = 0.6650 * _E_GEV**2 * B0_T * 1e3   # critical energy [eV]

    psi_half_mrad  = (slit_v_mm * 1e-3 / _DISTANCE) * 1e3 / 2.0  # half of total angle
    hdiv_full_mrad = 2.0 * K / gamma * 1e3
    hdiv_slit_mrad = (slit_h_mm * 1e-3 / _DISTANCE) * 1e3
    hdiv_used_mrad = min(hdiv_slit_mrad, hdiv_full_mrad)

    energy = numpy.linspace(_E_MIN_EV, _E_MAX_EV, _N_POINTS)

    flux = numpy.array(sync_ene(
        f_psi       = 2,
        energy_ev   = energy,
        ec_ev       = EC_EV,
        e_gev       = _E_GEV,
        i_a         = _I_A,
        hdiv_mrad   = hdiv_used_mrad,
        psi_min     = -psi_half_mrad,
        psi_max     = psi_half_mrad,
        psi_npoints = 50,
    )).flatten() * 2 * n_periods

    eV_to_J        = 1.60218e-19
    spectral_power = flux * eV_to_J / 1e-3           # W/eV
    de             = energy[1] - energy[0]
    cumulated_power = numpy.cumsum(spectral_power) * de  # W

    return energy, flux, spectral_power, cumulated_power

#######################################################################################



def get_att_transmission(attenuators_json, attenuators_up_to_axis=None, energy=None, verbose=1, do_plot=0):
    from syned.util.json_tools import load_from_json_file, load_from_json_url

    syned_filterbox = load_from_json_file(attenuators_json)

    if verbose:
        print("**** attenuators ****")
        print("Using file: ", attenuators_json)
        print("Number of blocks in file: ", syned_filterbox.get_n())

    if attenuators_up_to_axis is not None: syned_filterbox.set_n(attenuators_up_to_axis + 1)
    materials, thicknesses, densities = syned_filterbox.get_lists_materials_thicknesses_densities(cumulate=1)

    if verbose:
        print("Number of blocks used: ", syned_filterbox.get_n())
        for i in range(len(materials)):
            print(i, materials[i], thicknesses[i], densities[i])


    #
    # script to make the calculations (created by XOPPY:xpower)
    #
    import numpy
    from xoppylib.power.xoppy_calc_power import xoppy_calc_power
    try: import xraylib
    except: print("xraylib not available")
    from dabax.dabax_xraylib import DabaxXraylib

    out_dictionary = xoppy_calc_power(
            energy,
            numpy.ones_like(energy),
            substance = materials,
            thick     = thicknesses, # in mm (for filters)
            angle     = [6,] * len(materials), # in mrad (for mirrors)
            dens      = densities,
            roughness = [0,]  * len(materials), # in A (for mirrors)
            flags     = [0,]  * len(materials), # 0=Filter, 1=Mirror
            nelements = len(materials),
            FILE_DUMP = 0,
            material_constants_library = DabaxXraylib(file_f1f2="f1f2_Windt.dat",file_CrossSec="CrossSec_EPDL97.dat"),
            )

    # data to pass
    energy = out_dictionary["data"][0,:]
    spectral_power = out_dictionary["data"][-1,:]

    numpy.savetxt("cumulated_filters.txt", numpy.column_stack([energy, spectral_power]))
    print("File cumulated_filters.txt written to disk.")

    #
    # example plots
    #
    if do_plot:
        from srxraylib.plot.gol import plot
        plot(out_dictionary["data"][0,:], out_dictionary["data"][1,:],
            out_dictionary["data"][0,:], out_dictionary["data"][-1,:],
            xtitle=out_dictionary["labels"][0],
            legend=[out_dictionary["labels"][1],out_dictionary["labels"][-1]],
            title='Spectral Power [W/eV]')

    return energy, spectral_power



def get_mono_absorption(energy, thick=2.0, file_CrossSec="CrossSec_NIST_MassEnergyAbsorption.dat"):
    import numpy
    spectral_power = numpy.ones_like(energy)

    #
    # script to make the calculations (created by XOPPY:xpower)
    #

    import numpy
    from xoppylib.power.xoppy_calc_power import xoppy_calc_power
    try:
        import xraylib
    except:
        print("xraylib not available")
    from dabax.dabax_xraylib import DabaxXraylib

    out_dictionary = xoppy_calc_power(
        energy,
        spectral_power,
        substance=['Si', ],
        thick=[thick, ],  # in mm (for filters)
        angle=[3, ],  # in mrad (for mirrors)
        dens=['?', ],
        roughness=[0, ],  # in A (for mirrors)
        flags=[0, ],  # 0=Filter, 1=Mirror
        nelements=1,
        FILE_DUMP=0,
        material_constants_library=DabaxXraylib(file_f1f2="f1f2_Windt.dat",
                                                file_CrossSec=file_CrossSec),
    )

    return 1.0 - spectral_power


def print_power(spectral_power, energy, title="Power [W]"):
    print(title, get_power(spectral_power, energy))

def get_power(spectral_power, energy):
    return (numpy.cumulative_sum(spectral_power) * (energy[1] - energy[0]))[-1]

_spectrum_cache = {}

def calculate_power(method=1, # 0=SRW, 1=WS
                    u='u20',  # 'u20', 'u18' 'u22'
                    K=1.0,
                    slit_h_mm=1.8,
                    slit_v_mm=1.0,
                    attenuators_json='id11_wattdog_attenuators_2028_syned.json',
                    attenuators_up_to_axis=7,
                    file_CrossSec="CrossSec_NIST_MassEnergyAbsorption.dat",  # For Laue and Bragg
                    ):
    #
    # spectra
    #
    POWER = [] # to store results

    cache_key = (method, u, round(K, 6), round(slit_h_mm, 3), round(slit_v_mm, 3))
    if cache_key in _spectrum_cache:
        energy, flux, spectral_power, cumulated_power = _spectrum_cache[cache_key]
        print(f"  [cache hit] method={method} u={u} K={K:.4f}")
    else:
        energy, flux, spectral_power, cumulated_power = get_id_spectrum(K, u=u, method=method, slit_h_mm=slit_h_mm, slit_v_mm=slit_v_mm)
        _spectrum_cache[cache_key] = (energy, flux, spectral_power, cumulated_power)
    POWER.append(get_power(spectral_power, energy))

    #
    # attenuators
    #
    energy, transmission = get_att_transmission(attenuators_json,  attenuators_up_to_axis=attenuators_up_to_axis,
                                                energy=energy, verbose=1, do_plot=0)

    spectral_power_attenuated = spectral_power * transmission

    POWER.append(get_power(spectral_power_attenuated, energy)) # after attenuators

    spectral_power_absorbed_in_bragg = spectral_power_attenuated * get_mono_absorption(energy, thick=90.0, file_CrossSec=file_CrossSec)
    POWER.append(get_power(spectral_power_absorbed_in_bragg, energy)) # absorbed in laue


    spectral_power_absorbed_in_laue = spectral_power_attenuated * get_mono_absorption(energy, thick=2.0, file_CrossSec=file_CrossSec)
    POWER.append(get_power(spectral_power_absorbed_in_laue, energy)) # absorbed in laue

    return POWER


#
#
#

if '__main__' == __name__:

    if 0:
        #
        # settings
        #
        method = 3 # 0=SRW, 1=WS, 2=URGENT-HARMONICS, 3=BM approx
        file_CrossSec = "CrossSec_NIST_MassEnergyAbsorption.dat" # For Laue and Bragg
        # file_CrossSec = "CrossSec_EPDL97.dat"  # For Laue and Bragg


        #
        # parameters
        #
        u = 'u18' # 'u20' # 'u18' 'u22'
        gap_mm = 6.0
        K = K_vs_gap(gap_mm=gap_mm, u=u)
        slit_h_mm = 1.8
        slit_v_mm = 1.0
        attenuators_json = 'id11_wattdog_attenuators_2028_syned.json'
        attenuators_up_to_axis = 0  # 0 gets FE window, 6 uses all 6 attenuator axes


        #
        # results
        #
        if 0:
            for i in range(7):
                POWER = calculate_power(method=method,
                                        u=u,
                                        K=K,
                                        slit_h_mm=slit_h_mm,
                                        slit_v_mm=slit_v_mm,
                                        attenuators_json=attenuators_json,
                                        attenuators_up_to_axis=i,
                                        file_CrossSec=file_CrossSec,
                                        )

                print("\n\n---------------------------------------------------------------")
                print("\n\n%s(%s), K=%.3f (gap=%.1f mm) %.1fx%.1f mm^2@23m, atts up to: %d, mono CS from: %s" % \
                      (u, ["SRW", "WS"][method], K, gap_mm, slit_h_mm, slit_v_mm, attenuators_up_to_axis, file_CrossSec))
                POWER_TXT = ["\nPower at FE: %.3f W", "Power attenuated: %.3f W", "Absorbed in Bragg: %.3f W", "Absorbed in Laue: %.3f W"]
                for i in range(len(POWER)):
                    print(POWER_TXT[i] % (POWER[i]))



        POWER = calculate_power(method=method,
                                u=u,
                                K=K,
                                slit_h_mm=slit_h_mm,
                                slit_v_mm=slit_v_mm,
                                attenuators_json=attenuators_json,
                                attenuators_up_to_axis=attenuators_up_to_axis,
                                file_CrossSec=file_CrossSec,
                                )

        print("\n\n---------------------------------------------------------------")
        print("\n\n%s(%s), K=%.3f (gap=%.1f mm) %.1fx%.1f mm^2@23m, atts up to(0=FE window, 6=all): %d, mono CS from: %s" % \
              (u, ["SRW", "WS", "URGENT", "BM"][method], K, gap_mm, slit_h_mm, slit_v_mm, attenuators_up_to_axis, file_CrossSec))
        POWER_TXT = ["\nPower at FE (no window): %.3f W", "Power attenuated: %.3f W", "Absorbed in Bragg: %.3f W", "Absorbed in Laue: %.3f W"]
        for i in range(len(POWER)):
            print(POWER_TXT[i] % (POWER[i]))
        print("\n\n")


    #
    # make attenuators plot
    #
    if 1:
        from srxraylib.plot.gol import plot
        e = numpy.linspace(1000.0, 200000, 200)
        RESULT=[]
        for i in range(7):
            e, r = get_att_transmission('id11_wattdog_attenuators_2028_syned.json',
                                        attenuators_up_to_axis=i, energy=e, verbose=0)
            RESULT.append(r)

        e, r = get_att_transmission('id11_wattdog_attenuators_2026_syned.json',
                                    attenuators_up_to_axis=i, energy=e, verbose=0)
        RESULT.append(r)

        plot(e * 1e-3, RESULT[0],
             e * 1e-3, RESULT[1],
             e * 1e-3, RESULT[2],
             e * 1e-3, RESULT[3],
             e * 1e-3, RESULT[4],
             e * 1e-3, RESULT[5],
             e * 1e-3, RESULT[6],
             e * 1e-3, RESULT[7],
             legend=["2028.0","2028.1","2028.2","2028.3","2028.4","2028.5","2028.6","2026"],
             xtitle="Energy [keV]",ytitle="Attenuators transmission", grid=1,
             linestyle=[None,None,None,None,None,None,None,":"],
             )

    from syned.util.json_tools import load_from_json_file, load_from_json_url
    syned_filterbox = load_from_json_file('id11_wattdog_attenuators_2028_syned.json')


    for i in range(7):
        syned_filterbox_i = syned_filterbox.duplicate()
        syned_filterbox_i.set_n(i + 1)
        materials, thicknesses, densities = syned_filterbox_i.get_lists_materials_thicknesses_densities(cumulate=1)
        print("2028.%d " % i, materials, thicknesses)

    #
    # make spectra plot
    #
    if 0:
        # energy, flux, spectral_power, cumulated_power
        u='u18'
        Kmax = {'u18': 1.6563, 'u20': 2.334, 'u22': 1.5426}
        N = {'u18': 111, 'u20': 98, 'u22': 91}
        PERIOD = {'u18': 0.018, 'u20': 0.0205, 'u22': 0.022}
        s0 = get_id_spectrum(K=Kmax[u], u=u, method=0, slit_h_mm=1.8 * 10, slit_v_mm=1.0 * 5)
        s1 = get_id_spectrum(K=Kmax[u], u=u, method=1, slit_h_mm=1.8 * 10, slit_v_mm=1.0 * 5)
        s2 = get_id_spectrum(K=Kmax[u], u=u, method=2, slit_h_mm=1.8 * 10, slit_v_mm=1.0 * 5)
        s3 = get_id_spectrum(K=Kmax[u], u=u, method=3, slit_h_mm=1.8 * 10, slit_v_mm=1.0 * 5)
        # s0 = get_id_spectrum(K=Kmax[u], u=u, method=0, slit_h_mm=1.8, slit_v_mm=1.0)
        # s1 = get_id_spectrum(K=Kmax[u], u=u, method=1, slit_h_mm=1.8, slit_v_mm=1.0)
        # s2 = get_id_spectrum(K=Kmax[u], u=u, method=2, slit_h_mm=1.8, slit_v_mm=1.0)
        # s3 = get_id_spectrum(K=Kmax[u], u=u, method=3, slit_h_mm=1.8, slit_v_mm=1.0)

        # total_power = 7.257e-2 * 6**2 * Kmax[u]**2 * N[u] * 0.2 / PERIOD[u]

        # Formula 1: via K
        total_power = 7.257e-5 * 6 ** 2 * Kmax[u] ** 2 * N[u] * 0.2 / PERIOD[u]
        print(f"P (K formula) = {total_power:.2f} kW")

        # Formula 2: 2N bending magnets
        B = Kmax[u] / 0.9337 / (PERIOD[u] * 1e2)  # peak field [T]
        gamma = 6000 / 0.511  # Lorentz factor
        theta = numpy.pi * Kmax[u] / gamma  # half-period bending angle [rad]
        total_power_bm = 2.110 * 6 ** 3 * B * theta * 2 * N[u] * 0.2
        print(f"P (BM formula) = {total_power_bm:.2f} kW")


        from srxraylib.plot.gol import plot

        legend = ["SRW", "WS", "URGENT", "BM"]
        plot(s0[0] * 1e-3, s0[2],
             s1[0] * 1e-3, s1[2],
             s2[0] * 1e-3, s2[2],
             s3[0] * 1e-3, s3[2],
             ylog=1,
             legend=legend, title="%s K=%.3f" % (u, Kmax[u]),
             xtitle="Energy [keV]",ytitle="Spectral Power [W/eV]", grid=1, show=0,
             )

        plot(s0[0] * 1e-3, s0[3],
             s1[0] * 1e-3, s1[3],
             s2[0] * 1e-3, s2[3],
             s3[0] * 1e-3, s3[3],
             ylog=0,
             legend=legend, title="%s K=%.3f Total power %.3f kW" % (u, Kmax[u], total_power),
             xtitle="Energy [keV]",ytitle="Cumulated Power [W]", grid=1, show=1,
             )