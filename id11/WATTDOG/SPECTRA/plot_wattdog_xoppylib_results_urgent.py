import numpy as np
import pylab as plt
import json
import os
from wattdog_xoppylib import calculate_power

RESULTS_FILE = "plot_wattdog_xoppylib_results_urgent.json"

devices     = ['u18', 'u20', 'u22']
attenuators = ['2028.0', '2028.1', '2028.2', '2028.3', '2028.4', '2028.5', '2028.6', '2026']
methods     = [0, 2]          # 0=SRW, 1=WS, 2=URGENT, 3=BM
method_labels = {0: 'SRW', 1: 'WS', 2: 'URGENT', 3: 'BM'}

# # srio: cpmu20 of id16 has K=2.334 at 6mm gap
# IDS = {'cpmu18': {'KV': 1.6563, 'gap_min': 6.0, 'nperiods': 111, 'period': 0.018},
#        'u22': {'KV': 1.5426, 'gap_min': 6.8, 'nperiods': 91, 'period': 0.022},
#        'cpmu20': {'KV': 2.334, 'gap_min': 6.0, 'nperiods': 98, 'period': 0.0205}}


Kmax  = {'u18': 1.6563, 'u20': 2.334, 'u22': 1.5426}
Kvs   = {u: np.append(np.arange(0.5, Kmax[u], 0.2), Kmax[u]) for u in devices}

file_CrossSec = "CrossSec_NIST_MassEnergyAbsorption.dat"

# results[device][method][attenuator][K] = (power_ll, power_bb)
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
        device, method_str, att, K_str = key.split("|")
        method = int(method_str)
        K = float(K_str)
        if device not in results:
            results[device] = {}
        if method not in results[device]:
            results[device][method] = {a: {} for a in attenuators}
        results[device][method][att][round(K, 1)] = tuple(val)
else:
    results = {}
    for method in methods:
        for device in devices:
            if device not in results:
                results[device] = {}
            results[device][method] = {att: {} for att in attenuators}
            for K in Kvs[device]:
                for att in attenuators:
                    if att == '2026':
                        attenuators_json = 'id11_wattdog_attenuators_2026_syned.json'
                        attenuators_up_to_axis = None
                    else:
                        attenuators_json = 'id11_wattdog_attenuators_2028_syned.json'
                        attenuators_up_to_axis = int(att.split(".")[1])

                    POWER = calculate_power(method=method, u=device, K=K,
                                            attenuators_json=attenuators_json,
                                            attenuators_up_to_axis=attenuators_up_to_axis,
                                            file_CrossSec=file_CrossSec)
                    # POWER[0] source POWER[1] after att POWER[2] absorbed BRAGG, POWER[3] absorbed LA
                    # POWER = [total, attenuated, absorbed_bragg, absorbed_laue]
                    # results[device][method][att][round(K, 1)] = (POWER[3], POWER[2])
                    results[device][method][att][round(K, 1)] = (POWER[3], POWER[1]) # SUPPOSE BRAGG ABSORBS EVERYTHING!

    # dump to file
    flat = {}
    for device in results:
        for method in results[device]:
            for att in results[device][method]:
                for K, val in results[device][method][att].items():
                    flat[f"{device}|{method}|{att}|{K}"] = list(val)
    with open(RESULTS_FILE, "w") as f:
        json.dump(flat, f, indent=2)
    print(f"Results saved to {RESULTS_FILE}")

# --------------------------------------------------------------------------
# plot
# --------------------------------------------------------------------------
f, a = plt.subplots(len(devices), len(attenuators),
                    constrained_layout=True, figsize=(18, 8))

for j, device in enumerate(devices):
    for i, att in enumerate(attenuators):
        for method in methods:
            lbl = method_labels[method]
            ll_vals = [results[device][method][att][round(K, 1)][0] for K in Kvs[device]]
            bb_vals = [results[device][method][att][round(K, 1)][1] for K in Kvs[device]]
            a[j, i].plot(Kvs[device], ll_vals, 'o-',  label=f'll {lbl}')
            a[j, i].plot(Kvs[device], bb_vals, 's--', label=f'bb {lbl}')

        a[j, i].set(title=f'{device}  {att}', ylim=(0, None), xlim=(0, None))
        # Add customized grid
        a[j, i].grid(True, linestyle='--', alpha=0.4, color='gray', linewidth=0.5)
        if i == 0:
            a[j, i].legend(loc='upper left', fontsize=6)

for j in range(len(devices)):
    a[j, 0].set(ylabel='Power / W')
for i in range(len(attenuators)):
    a[-1, i].set(xlabel='K')

plt.suptitle('WattDog — power absorbed in LL and BB crystals')

plt.tight_layout()
outpng = "plot_wattdog_xoppylib_results_urgent.png"
plt.savefig(outpng, dpi=150, bbox_inches='tight')
plt.show()
