import numpy as np
import pylab as plt
import json
import os
from wattdog_xoppylib import calculate_power

RESULTS_FILE = "plot_wattdog_xoppylib_slit_v_linearity.json"

devices     = ['u18', 'u20', 'u22']
attenuators = ['2028.0', '2028.1', '2028.2', '2028.3', '2028.4', '2028.5', '2028.6', '2026']
methods     = [0, 1]          # 0=SRW, 1=WS
method_labels = {0: 'SRW', 1: 'WS'}

Kmax     = {'u18': 1.6563, 'u20': 2.334, 'u22': 1.5426}
slit_h_mm = 1.8
slit_v_vals = np.linspace(0.001, 1.001, 5)

file_CrossSec = "CrossSec_NIST_MassEnergyAbsorption.dat"

# results[device][method][attenuator][slit_v] = (power_ll, power_bb)
# power_ll : power absorbed in LL crystal (2.5 mm Si)
# power_bb : power absorbed in Bragg crystal (90 mm Si)

# --------------------------------------------------------------------------
# load or calculate
# --------------------------------------------------------------------------
if os.path.exists(RESULTS_FILE):
    print(f"Loading results from {RESULTS_FILE}")
    with open(RESULTS_FILE) as f:
        flat = json.load(f)
    results = {}
    for key, val in flat.items():
        device, method_str, att, sh_str = key.split("|")
        method  = int(method_str)
        slit_v  = float(sh_str)
        if device not in results:
            results[device] = {}
        if method not in results[device]:
            results[device][method] = {a: {} for a in attenuators}
        results[device][method][att][round(slit_v, 4)] = tuple(val)
else:
    results = {}
    for method in methods:
        for device in devices:
            if device not in results:
                results[device] = {}
            results[device][method] = {att: {} for att in attenuators}
            for slit_v in slit_v_vals:
                for att in attenuators:
                    if att == '2026':
                        attenuators_json      = 'id11_wattdog_attenuators_2026_syned.json'
                        attenuators_up_to_axis = None
                    else:
                        attenuators_json      = 'id11_wattdog_attenuators_2028_syned.json'
                        attenuators_up_to_axis = int(att.split(".")[1])

                    POWER = calculate_power(method=method, u=device, K=Kmax[device],
                                            slit_h_mm=slit_h_mm, slit_v_mm=slit_v,
                                            attenuators_json=attenuators_json,
                                            attenuators_up_to_axis=attenuators_up_to_axis,
                                            file_CrossSec=file_CrossSec)
                    # POWER = [total, attenuated, absorbed_bragg, absorbed_laue]
                    # results[device][method][att][round(slit_v, 4)] = (POWER[3], POWER[2])
                    results[device][method][att][round(slit_v, 4)] = (POWER[3], POWER[1]) # BRAGG absorbes everything!!

    # dump to file
    flat = {}
    for device in results:
        for method in results[device]:
            for att in results[device][method]:
                for slit_v, val in results[device][method][att].items():
                    flat[f"{device}|{method}|{att}|{slit_v}"] = list(val)
    with open(RESULTS_FILE, "w") as f:
        json.dump(flat, f, indent=2)
    print(f"Results saved to {RESULTS_FILE}")

# --------------------------------------------------------------------------
# linearity test: R² of linear fit y = m*x + b
# --------------------------------------------------------------------------
def r2(x, y):
    p = np.polyfit(x, y, 1)
    y_fit = np.polyval(p, x)
    ss_res = np.sum((np.asarray(y) - y_fit) ** 2)
    ss_tot = np.sum((np.asarray(y) - np.mean(y)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

# --------------------------------------------------------------------------
# plot
# --------------------------------------------------------------------------
f, a = plt.subplots(len(devices), len(attenuators),
                    constrained_layout=True, figsize=(18, 8))

for j, device in enumerate(devices):
    for i, att in enumerate(attenuators):
        r2_lines = []
        for method in methods:
            lbl = method_labels[method]
            ll_vals = [results[device][method][att][round(sh, 4)][0] for sh in slit_v_vals]
            bb_vals = [results[device][method][att][round(sh, 4)][1] for sh in slit_v_vals]
            r2_ll = r2(slit_v_vals, ll_vals)
            r2_bb = r2(slit_v_vals, bb_vals)
            r2_lines.append(f'll {lbl} R²={r2_ll:.4f}')
            r2_lines.append(f'bb {lbl} R²={r2_bb:.4f}')
            a[j, i].plot(slit_v_vals, ll_vals, 'o-',  label=f'll {lbl}')
            a[j, i].plot(slit_v_vals, bb_vals, 's--', label=f'bb {lbl}')

        a[j, i].set(title=f'{device} K={Kmax[device]}  {att}', ylim=(0, None), xlim=(0, None))
        a[j, i].grid(True, linestyle='--', alpha=0.4, color='gray', linewidth=0.5)
        a[j, i].text(0.02, 0.98, '\n'.join(r2_lines), transform=a[j, i].transAxes,
                     fontsize=7, va='top', ha='left', family='monospace')
        if i == 0:
            a[j, i].legend(loc='lower right', fontsize=6)

for j in range(len(devices)):
    a[j, 0].set(ylabel='Power / W')
for i in range(len(attenuators)):
    a[-1, i].set(xlabel='slit V / mm')

plt.suptitle('WattDog — slit V linearity at Kmax, power absorbed in LL and Bragg crystals')

plt.tight_layout()
outpng = "plot_wattdog_xoppylib_slit_v_linearity"
plt.savefig(outpng, dpi=150, bbox_inches='tight')


plt.show()
