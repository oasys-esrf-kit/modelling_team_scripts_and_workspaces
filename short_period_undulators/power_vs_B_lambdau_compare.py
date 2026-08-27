import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import scipy.constants as codata
import csv

#
# inputs (same storage-ring / beam parameters as example.py)
#
ring_current            = 0.2     # A
electron_energy_in_GeV  = 6       # GeV
undulator_length        = 2.0     # m   (used for the fixed-L case)
number_of_periods_fixed = 138.889 # periods (used for the fixed-N case)

period_min_mm, period_max_mm = 3, 30    # mm
B_min, B_max = 0.0, 3.0                 # T
n_period, n_B = 200, 200

period_mm = numpy.linspace(period_min_mm, period_max_mm, n_period)
B         = numpy.linspace(B_min, B_max, n_B)
period_m  = period_mm * 1e-3

PERIOD, BFIELD = numpy.meshgrid(period_m, B, indexing='ij')
K = BFIELD * PERIOD * codata.e / (2 * numpy.pi * codata.m_e * codata.c)
gamma1 = 1e9 * electron_energy_in_GeV / (codata.m_e * codata.c**2 / codata.e)
Z0 = codata.value('characteristic impedance of vacuum')

def ptot_from_N(N):
    return (N / 6) * Z0 * ring_current * codata.e * 2 * numpy.pi * codata.c * gamma1 ** 2 * K ** 2 / PERIOD

# case 1: fixed undulator length L -> N = L / period varies with period
ptot_fixedL = ptot_from_N(undulator_length / PERIOD) / 1e3   # kW

# case 2: fixed number of periods N -> undulator length L = N * period varies with period
ptot_fixedN = ptot_from_N(number_of_periods_fixed) / 1e3     # kW

vmin, vmax = 0, 40

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
        [ptot_fixedL, ptot_fixedN],
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
fig.suptitle("Total power vs Period and Magnetic field (E = %.0f GeV, I = %.2f A)" %
             (electron_energy_in_GeV, ring_current), y=1.05)
fig.colorbar(im, ax=axes, label="Total power [kW]", fraction=0.046, pad=0.04)

plt.savefig("power_vs_period_and_B_compare.png", dpi=150, bbox_inches='tight')
plt.close()
print("saved power_vs_period_and_B_compare.png")
