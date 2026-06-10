import numpy as np
import pylab as plt
import json
import os
from wattdog_light import calculate_power

RESULTS_FILE = "plot_wattdog_light_results_bm.json"

devices     = ['u18', 'u20', 'u22']
attenuators = ['2028.0', '2028.1', '2028.2', '2028.3', '2028.4', '2028.5', '2028.6', '2026']

Kmax  = {'u18': 1.6563, 'u20': 2.334, 'u22': 1.5426}
Kvs   = {u: np.append(np.arange(0.5, Kmax[u], 0.2), Kmax[u]) for u in devices}

# results[device][attenuator][K] = (power_ll, power_bb)

# --------------------------------------------------------------------------
# load or calculate
# --------------------------------------------------------------------------
if os.path.exists(RESULTS_FILE):
    print(f"Loading results from {RESULTS_FILE}")
    with open(RESULTS_FILE) as f:
        flat = json.load(f)
    results = {}
    for key, val in flat.items():
        device, att, K_str = key.split("|")
        K = float(K_str)
        if device not in results:
            results[device] = {a: {} for a in attenuators}
        results[device][att][round(K, 1)] = tuple(val)
else:
    results = {}
    for device in devices:
        results[device] = {att: {} for att in attenuators}
        for K in Kvs[device]:
            for att in attenuators:
                pBB, pLL = calculate_power(u=device, K=K, attenuators_label=att)
                results[device][att][round(K, 1)] = (pLL, pBB)

    # dump to file
    flat = {}
    for device in results:
        for att in results[device]:
            for K, val in results[device][att].items():
                flat[f"{device}|{att}|{K}"] = list(val)
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
        ll_vals = [results[device][att][round(K, 1)][0] for K in Kvs[device]]
        bb_vals = [results[device][att][round(K, 1)][1] for K in Kvs[device]]
        a[j, i].plot(Kvs[device], ll_vals, 'o-',  label='ll light')
        a[j, i].plot(Kvs[device], bb_vals, 's--', label='bb light')

        a[j, i].set(title=f'{device}  {att}', ylim=(0, None), xlim=(0, None))
        # Add customized grid
        a[j, i].grid(True, linestyle='--', alpha=0.4, color='gray', linewidth=0.5)
        if i == 0:
            a[j, i].legend(loc='upper left', fontsize=6)

for j in range(len(devices)):
    a[j, 0].set(ylabel='Power / W')
for i in range(len(attenuators)):
    a[-1, i].set(xlabel='K')

plt.suptitle('WattDog (light, BM) — power on BB and after LL filter')

plt.tight_layout()
outpng = "plot_wattdog_light_results_bm.png"
plt.savefig(outpng, dpi=150, bbox_inches='tight')
plt.show()
