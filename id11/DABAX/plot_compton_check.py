#
# Visual check of the Compton correction applied by add_compton_all.py:
#   - ratio = f2_corrected / f2_original for Cu and Ge (from the actual
#     Windt table, capped at its 100 keV upper limit, plus a curve computed
#     directly from xraylib extending up to 300 keV)
#   - mu (linear attenuation coefficient, cm^-1): Photoelectric and Compton
#     contributions plotted separately, computed directly from xraylib
#     (CS_Photo, CS_Compt), independent of the Windt f2 table's 100 keV limit
#
import re
import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xraylib


R_E_CM = 2.8179403262e-13  # classical electron radius, cm
N_A = 6.02214076e23        # 1/mol
HC_EV_ANGSTROM = 12398.42  # eV * Angstrom

ELEMENTS = ["Cu", "Ge"]
EMAX_EXT_KEV = 300.0


def get_section(text, symbol):
    pat = re.compile(r"(#S\s+\d+\s+" + symbol + r"\b.*?)(?=\n#S\s+\d+\s+|\Z)", re.S)
    return pat.search(text).group(1)


def parse_data(section_text):
    lines = section_text.splitlines()
    data_lines = [l for l in lines if l.strip() and not l.startswith("#")]
    data = np.array([[float(v) for v in line.split()] for line in data_lines])
    return data[:, 0], data[:, 1], data[:, 2]  # E_eV, f1, f2


def xraylib_curves(Z, rho, E_keV):
    """Photoelectric and Compton mu [cm^-1], plus the ratio used for the
    f2 correction, evaluated directly from xraylib (valid well beyond the
    100 keV limit of the Windt f2 table)."""
    ratio = np.full_like(E_keV, np.nan)
    mu_photo = np.full_like(E_keV, np.nan)
    mu_compt = np.full_like(E_keV, np.nan)
    last_ratio = 1.0
    for i, e in enumerate(E_keV):
        try:
            photo = xraylib.CS_Photo(Z, e)
            compt = xraylib.CS_Compt(Z, e)
        except Exception:
            continue
        last_ratio = (photo + compt) / photo if photo != 0 else last_ratio
        ratio[i] = last_ratio
        mu_photo[i] = photo * rho
        mu_compt[i] = compt * rho
    return ratio, mu_photo, mu_compt


orig_text = open("f1f2_Windt.dat").read()
comp_text = open("f1f2_WindtWithCompton.dat").read()

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

for row, symbol in enumerate(ELEMENTS):
    Z = xraylib.SymbolToAtomicNumber(symbol)
    rho = xraylib.ElementDensity(Z)
    A = xraylib.AtomicWeight(Z)

    Eo, f1o, f2o = parse_data(get_section(orig_text, symbol))
    Ec, f1c, f2c = parse_data(get_section(comp_text, symbol))
    assert np.allclose(Eo, Ec)

    ratio_table = f2c / f2o

    # xraylib-direct curves, independent of the Windt table's 100 keV limit
    E_ext_keV = np.linspace(Eo[0] * 1e-3, EMAX_EXT_KEV, 1000)
    ratio_xrl, mu_photo_xrl, mu_compt_xrl = xraylib_curves(Z, rho, E_ext_keV)

    ax_ratio = axes[row, 0]
    ax_ratio.plot(E_ext_keV, ratio_xrl, "-", color="gray", label="xraylib (photo+compt)/photo, to 300 keV")
    ax_ratio.plot(Eo * 1e-3, ratio_table, "o", ms=2, color="C1", label="applied ratio (f2_corr/f2_orig, table)")
    ax_ratio.set_xlim(0, EMAX_EXT_KEV)
    ax_ratio.set_xlabel("Photon energy [keV]")
    ax_ratio.set_ylabel("ratio")
    ax_ratio.set_title(f"{symbol}: Compton correction ratio")
    ax_ratio.legend(fontsize=8)
    ax_ratio.grid(True, alpha=0.3)

    ax_mu = axes[row, 1]
    ax_mu.plot(E_ext_keV, mu_photo_xrl, "-", color="C0", label=f"{symbol} Photoelectric")
    ax_mu.plot(E_ext_keV, mu_compt_xrl, "-", color="C1", label=f"{symbol} Compton")
    ax_mu.set_xlim(0, EMAX_EXT_KEV)
    ax_mu.set_ylim(0, 3.0)
    ax_mu.set_xlabel("Photon energy [keV]")
    ax_mu.set_ylabel("Attenuation Coefficient [cm$^{-1}$]")
    ax_mu.set_title(f"{symbol}: linear attenuation coefficient")
    ax_mu.legend(fontsize=8)
    ax_mu.grid(True, alpha=0.3)

plt.tight_layout()
outfile = "compton_check_Cu_Ge.png"
plt.savefig(outfile, dpi=150)
print("Written:", outfile)
plt.show()
