import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv

CSV_PATH = "esrf_undulators_kmax.csv"

rows = []
with open(CSV_PATH, newline='') as f:
    reader = csv.DictReader(f)
    for r in reader:
        if r['Kmax'] == '' or r['Kmax'] is None:
            continue
        rows.append(dict(name=r['Name'], period=float(r['Period [mm]']), kmax=float(r['Kmax'])))

periods = [r['period'] for r in rows]
kmaxs = [r['kmax'] for r in rows]
names = [r['name'] for r in rows]
colors = ['#C0504D' if 'CPMU' in name else '#4F81BD' for name in names]

SHORT_PERIOD_MAX_MM = 40

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, xlim, title in zip(
        axes,
        [None, (10, SHORT_PERIOD_MAX_MM)],
        ["All devices", "Short-period devices (period <= %d mm)" % SHORT_PERIOD_MAX_MM]):
    ax.scatter(periods, kmaxs, s=45, color=colors, edgecolor='black', linewidth=0.5, zorder=3)
    ax.scatter([], [], s=45, color='#C0504D', edgecolor='black', linewidth=0.5, label='CPMU')
    ax.scatter([], [], s=45, color='#4F81BD', edgecolor='black', linewidth=0.5, label='other')
    ax.legend(loc='upper left', fontsize=8)
    for name, p, k in zip(names, periods, kmaxs):
        if xlim is None or p <= xlim[1]:
            ax.annotate(name, (p, k), fontsize=6, xytext=(3, 3), textcoords='offset points')
    if xlim is not None:
        ax.set_xlim(*xlim)
        sub_k = [k for p, k in zip(periods, kmaxs) if xlim[0] <= p <= xlim[1]]
        ax.set_ylim(0, max(sub_k) * 1.15)
    ax.set_xlabel("Period [mm]")
    ax.set_ylabel("Kmax (at Gap min)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

fig.suptitle("ESRF undulators: Kmax vs Period (jsrund.csv)")
plt.tight_layout()
plt.savefig("esrf_kmax_vs_period.png", dpi=150)
plt.close()
print("saved esrf_kmax_vs_period.png (%d points)" % len(rows))
