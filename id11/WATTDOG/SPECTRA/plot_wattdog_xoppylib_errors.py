import numpy as np
import pylab as plt
import json

devices     = ['u18', 'u20', 'u22']
attenuators = ['2028.0', '2028.1', '2028.2', '2028.3', '2028.4', '2028.5', '2028.6', '2026']
Kmax = {'u18': 1.6563, 'u20': 2.334, 'u22': 1.5426}
Kvs  = {u: np.append(np.arange(0.5, Kmax[u], 0.2), Kmax[u]) for u in devices}

# approximation methods to compare against SRW (method=0)
approx_methods = {1: 'WS', 2: 'URGENT', 3: 'BM'}
json_files     = {1: 'plot_wattdog_xoppylib_results_ws.json',
                  2: 'plot_wattdog_xoppylib_results_urgent.json',
                  3: 'plot_wattdog_xoppylib_results_bm.json'}

# --------------------------------------------------------------------------
# load all three JSON files into a single results dict
# results[device][method][att][K] = (power_ll, power_bb)
# --------------------------------------------------------------------------
results = {}

def load_json(path):
    with open(path) as f:
        return json.load(f)

for approx_method, fname in json_files.items():
    flat = load_json(fname)
    for key, val in flat.items():
        device, method_str, att, K_str = key.split("|")
        method = int(method_str)
        K = float(K_str)
        results.setdefault(device, {}).setdefault(method, {}).setdefault(att, {})[round(K, 1)] = tuple(val)

# --------------------------------------------------------------------------
# plot: % error = (approx - SRW) / SRW * 100  for ll and bb
# --------------------------------------------------------------------------
colors  = {1: 'tab:blue', 2: 'tab:orange', 3: 'tab:green'}
markers = {1: 'o',        2: 's',          3: '^'}

f, a = plt.subplots(len(devices), len(attenuators),
                    constrained_layout=True, figsize=(18, 8))

for j, device in enumerate(devices):
    Ks = [round(K, 1) for K in Kvs[device]]
    for i, att in enumerate(attenuators):
        srw = results[device][0][att]
        for method, lbl in approx_methods.items():
            if method not in results[device]:
                continue
            approx = results[device][method][att]
            ll_err, bb_err = [], []
            for K in Ks:
                srw_ll, srw_bb = srw[K]
                ap_ll,  ap_bb  = approx[K]
                ll_err.append((ap_ll - srw_ll) / srw_ll * 100 if srw_ll else np.nan)
                bb_err.append((ap_bb - srw_bb) / srw_bb * 100 if srw_bb else np.nan)
            c, m = colors[method], markers[method]
            a[j, i].plot(Kvs[device], ll_err, marker=m, ls='-',  color=c, label=f'll {lbl}')
            a[j, i].plot(Kvs[device], bb_err, marker=m, ls='--', color=c, label=f'bb {lbl}')

        a[j, i].axhline(0, color='k', lw=0.7, ls=':')
        a[j, i].set(title=f'{device}  {att}', xlim=(0, None))
        # Add customized grid
        a[j, i].grid(True, linestyle='--', alpha=0.4, color='gray', linewidth=0.5)
        if i == 0:
            a[j, i].legend(loc='best', fontsize=5)

for j in range(len(devices)):
    a[j, 0].set(ylabel='Error vs SRW / %')
for i in range(len(attenuators)):
    a[-1, i].set(xlabel='K')

plt.suptitle('WattDog — relative error of WS / URGENT / BM vs SRW (ll and bb)')

plt.tight_layout()
outpng = "plot_wattdog_xoppylib_errors.png"
plt.savefig(outpng, dpi=150, bbox_inches='tight')

plt.show()
