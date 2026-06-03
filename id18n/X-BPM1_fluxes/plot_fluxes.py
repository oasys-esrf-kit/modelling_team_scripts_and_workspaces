from srxraylib.plot.gol import plot
from scipy import constants as con
import matplotlib.pyplot as plt
import numpy as np

# Compare incoming fluxes on XBPM for
# 1) Worst case (full ML beam, unfiltered)
# 2) Normal working case (beam cut by secondary slit, filter xxx)
objects=[in_object_1,in_object_2,in_object_3,in_object_4,in_object_5,in_object_6]
e=np.zeros((7501,len(objects)))
flux=np.zeros((7501,len(objects)))

for i, obj in enumerate(objects):
    print(i)
    ei = obj.get_content('xoppy_data')[:,0]
    sp_powi = obj.get_content('xoppy_data')[:,-1]
    fluxi = sp_powi/con.e/1000
    
    print('filling e')
    e[:,i]=ei
    print('filling flux')
    flux[:,i]=fluxi


# flux = sp_pow/con.e/1000
# Plots
fig, ax = plot(e[:,0],flux[:,0], e[:,1],flux[:,1], e[:,2],flux[:,2],e[:,3],flux[:,3],e[:,4],flux[:,4],e[:,5],flux[:,5],
    title='ID18 nano, Flux on X-BPM 1',
    legend=['Worst Case  7keV', 'Working Case  7 keV',
            'Worst Case 17keV', 'Working Case 17 keV',
            'Worst Case 35keV', 'Working Case 35 keV'],
    xtitle='Energy (eV)', ytitle='flux (ph/s/0.1%bw)',
    color=['tab:blue','tab:orange','tab:blue','tab:orange','tab:blue','tab:orange'],
    linestyle=['-','-','--','--','-.','-.'],
    xrange=[5000,40000],yrange=[0,2.5e15],
    grid=True, ylog=False)
    
fig.set_constrained_layout(True)
leg = ax.get_legend()
leg.set_loc('upper center')
leg.set_bbox_to_anchor([0.5,-0.15])


# inset Axes....
x1, x2, y1, y2 = 5000, 8000, 0, 2.4e15  # subregion of the original image
axins = ax.inset_axes(
    [0.5, 0.2, 0.17, 0.77],
    xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
# Define the colors matching your main plot
colors = ['tab:blue', 'tab:orange', 'tab:blue', 'tab:orange', 'tab:blue', 'tab:orange']
linestyles = ['-', '-', '--', '--', '-.', '-.']

# Loop through the 6 objects to plot them cleanly
for j in range(len(objects)):
    axins.plot(e[:, j], flux[:, j], color=colors[j], linestyle=linestyles[j])

ax.indicate_inset_zoom(axins, edgecolor="black")
ax.indicate_inset_zoom(axins, edgecolor="black")

fig.show()
