import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import scipy.constants as codata
import csv
import io
import os
import contextlib

from xoppylib.sources.calc2d_kim import calc2d_kim

#
# inputs (same storage-ring / beam parameters as example.py and example_power_on_slit.py)
#
ring_current            = 0.2     # A
electron_energy_in_GeV  = 6.0     # GeV
electron_energy_spread  = 0.001
electron_beam_size_h    = 3.335657e-05   # m
electron_beam_size_v    = 7.28023e-06    # m
electron_beam_div_h     = 4.52532e-06    # rad
electron_beam_div_v     = 1.93834e-06    # rad

undulator_length        = 2.0     # m   (used for the fixed-L case)
number_of_periods_fixed = 138.889 # periods (used for the fixed-N case)

# front-end slit, as in example_power_on_slit.py
slit_distance_m = 23.0   # m
slit_gap_h_m    = 0.002  # m (2.0 mm)
slit_gap_v_m    = 0.001  # m (1.0 mm)
slit_h_points   = 21
slit_v_points   = 11

period_min_mm, period_max_mm = 3, 30    # mm
B_min, B_max = 0.0, 3.0                 # T
n_period, n_B = 40, 40

period_mm = numpy.linspace(period_min_mm, period_max_mm, n_period)
B         = numpy.linspace(B_min, B_max, n_B)
period_m  = period_mm * 1e-3

gamma1 = 1e9 * electron_energy_in_GeV / (codata.m_e * codata.c**2 / codata.e)


def power_on_slit_grid(n_periods_of_period):
    """n_periods_of_period(period_m) -> number of periods for that period.
    Returns a (n_period, n_B) array of power [W] through the slit."""
    out = numpy.zeros((n_period, n_B))
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        for ip, p_m in enumerate(period_m):
            N = n_periods_of_period(p_m)
            for ib, Bv in enumerate(B):
                K = Bv * p_m * codata.e / (2 * numpy.pi * codata.m_e * codata.c)
                bl = dict(ElectronBeamDivergenceH=electron_beam_div_h,
                          ElectronBeamDivergenceV=electron_beam_div_v,
                          ElectronBeamSizeH=electron_beam_size_h,
                          ElectronBeamSizeV=electron_beam_size_v,
                          ElectronCurrent=ring_current,
                          ElectronEnergy=electron_energy_in_GeV,
                          ElectronEnergySpread=electron_energy_spread,
                          Kv=max(K, 1e-6), Kh=0.0, Kphase=0.0,
                          NPeriods=N, PeriodID=p_m,
                          distance=slit_distance_m, gapH=slit_gap_h_m, gapV=slit_gap_v_m)
                h, v, p = calc2d_kim(bl, zero_emittance=False, fileName=None,
                                      hSlitPoints=slit_h_points, vSlitPoints=slit_v_points)
                out[ip, ib] = p.sum() * (h[1] - h[0]) * (v[1] - v[0])
    return out


CACHE_FILE = "power_on_slit_grid_cache_%dx%d.npz" % (n_period, n_B)
if os.path.exists(CACHE_FILE):
    _cache = numpy.load(CACHE_FILE)
    power_fixedL, power_fixedN = _cache['power_fixedL'], _cache['power_fixedN']
else:
    # case 1: fixed undulator length L -> N = L / period varies with period
    power_fixedL = power_on_slit_grid(lambda p_m: undulator_length / p_m)
    # case 2: fixed number of periods N
    power_fixedN = power_on_slit_grid(lambda p_m: number_of_periods_fixed)
    numpy.savez(CACHE_FILE, power_fixedL=power_fixedL, power_fixedN=power_fixedN)

vmin, vmax = 0, 5000

#
# resonance energy hyperbolas: on-axis fundamental E1 = (h*c/e) / [period/(2*gamma1^2) * (1+K^2/2)].
# Solve for K(period) such that E1 = target_energy_eV, for each target energy.
#
target_energies_eV = [5000, 10000, 15000, 20000]  # 5, 10, 15, 20 keV
hc_eVm = codata.h * codata.c / codata.e

resonance_lines = {}
for target_energy_eV in target_energies_eV:
    wavelength_target = hc_eVm / target_energy_eV
    K2 = 2 * (wavelength_target * 2 * gamma1 ** 2 / period_m - 1)
    K_line = numpy.where(K2 >= 0, numpy.sqrt(numpy.clip(K2, 0, None)), numpy.nan)
    B_line = K_line * 2 * numpy.pi * codata.m_e * codata.c / (codata.e * period_m)
    resonance_lines[target_energy_eV] = B_line

#
# ESRF undulators (Name, Period, Bmax at Gap min) from esrf_undulators_kmax.csv,
# split into CPMU vs other devices
#
esrf_periods, esrf_Bmax, esrf_colors = [], [], []
with open("esrf_undulators_kmax.csv", newline='') as f:
    for r in csv.DictReader(f):
        if r['Bmax [T]'] == '':
            continue
        esrf_periods.append(float(r['Period [mm]']))
        esrf_Bmax.append(float(r['Bmax [T]']))
        esrf_colors.append('#C0504D' if 'CPMU' in r['Name'] else '#4F81BD')

fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)

for ax, data, title in zip(
        axes,
        [power_fixedL, power_fixedN],
        ["Fixed undulator length\n(L = %.1f m)" % undulator_length,
         "Fixed number of periods\n(N = %.0f)" % number_of_periods_fixed]):
    im = ax.imshow(data.T, origin='lower', aspect='auto', cmap='viridis', vmin=vmin, vmax=vmax,
                    extent=[period_mm[0], period_mm[-1], B[0], B[-1]])
    for target_energy_eV in target_energies_eV:
        line = resonance_lines[target_energy_eV]
        ax.plot(period_mm, line, color='white', linestyle='--', linewidth=1.2)
        finite = numpy.isfinite(line) & (line <= B_max)
        if finite.any():
            i_lbl = numpy.where(finite)[0][-1]
            ax.annotate("%d keV" % (target_energy_eV / 1000), (period_mm[i_lbl], line[i_lbl]),
                        color='white', fontsize=8, ha='right', va='bottom')
    ax.scatter(esrf_periods, esrf_Bmax, s=35, color=esrf_colors, edgecolor='black',
               linewidth=0.5, zorder=4)
    ax.set_xlim(period_min_mm, period_max_mm)
    ax.set_ylim(B_min, B_max)
    ax.set_xlabel("Period [mm]")
    ax.set_title(title)

axes[0].scatter([], [], s=35, color='#C0504D', edgecolor='black', linewidth=0.5, label='CPMU')
axes[0].scatter([], [], s=35, color='#4F81BD', edgecolor='black', linewidth=0.5, label='other ESRF undulator')
axes[0].legend(loc='upper left', fontsize=8, framealpha=0.9)
axes[0].set_ylabel("Magnetic field [T]")
fig.suptitle("Power in front-end slit (%.1f x %.1f mm HxV @ %.0f m) vs Period and Magnetic field\n"
             "(E = %.0f GeV, I = %.2f A)" %
             (slit_gap_h_m * 1e3, slit_gap_v_m * 1e3, slit_distance_m,
              electron_energy_in_GeV, ring_current), y=1.08)
fig.colorbar(im, ax=axes, label="Power in slit [W]", fraction=0.046, pad=0.04)

plt.savefig("power_on_slit_vs_period_and_B_compare.png", dpi=150, bbox_inches='tight')
plt.close()
print("saved power_on_slit_vs_period_and_B_compare.png")
