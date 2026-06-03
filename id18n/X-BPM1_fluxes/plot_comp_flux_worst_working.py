from srxraylib.plot.gol import plot
from scipy import constants as con
import matplotlib.pyplot as plt
import numpy as np

# Compare incoming fluxes on XBPM for
# 1) Worst case (full ML beam, unfiltered)
# 2) Normal working case (beam cut by secondary slit, filter xxx)

e1 = in_object_1.get_content('xoppy_data')[:,0]
#flux1 = in_object_1.get_content('xoppy_data')[:,1]
sp_pow1 = in_object_1.get_content('xoppy_data')[:,-1]
flux1 = sp_pow1/con.e/1000

e2 = in_object_2.get_content('xoppy_data')[:,0]
sp_pow2 = in_object_2.get_content('xoppy_data')[:,-1]
flux2 = sp_pow2/con.e/1000
# flux = sp_pow/con.e/1000
# Plots
fig, ax = plot(e1,flux1, e2,flux2, 
    title='ID18 nano, E=35keV, Flux on X-BPM 1',
    legend=['Worst Case (full multilayer beam)', 'Working Case (beam cut by secondary slit & filtered with 0.6mm C + 0.9µm Au)'],
    xtitle='Energy (eV)', ytitle='flux (ph/s/0.1%bw)',
    xrange=[5000,40000],
    grid=True, ylog=False)
    
fig.set_constrained_layout(True)
leg = ax.get_legend()
leg.set_loc('upper center')
leg.set_bbox_to_anchor([0.5,-0.15])
fig.show()

# Export data
if np.array_equal(e1,e2):
    out_data=np.zeros((len(e1),3))
    out_data[:,0]=e1
    out_data[:,1]=flux1
    out_data[:,2]=flux2
    np.savetxt("data_out.csv", out_data, delimiter=",")