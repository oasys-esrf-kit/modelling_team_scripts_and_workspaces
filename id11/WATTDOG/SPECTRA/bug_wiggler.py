
#
# script to make the calculations (created by XOPPY:wiggler)
#
from xoppylib.sources.xoppy_bm_wiggler import xoppy_calc_wiggler_on_aperture
energy, flux, spectral_power, cumulated_power, traj, traj_info =  xoppy_calc_wiggler_on_aperture(
    FIELD=0,
    NPERIODS=111,
    ULAMBDA=0.018,
    K=1.6563,
    ENERGY=6.0,
    PHOT_ENERGY_MIN=1000.0,
    PHOT_ENERGY_MAX=200000.0,
    NPOINTS=500,
    NTRAJPOINTS=51,
    CURRENT=200.0,
    FILE="?",
    SLIT_FLAG=1,
    SLIT_D=23.0,
    SLIT_NY=30,
    SLIT_WIDTH_H_MM=1.8,
    SLIT_HEIGHT_V_MM=1.0,
    SLIT_CENTER_H_MM=0.0,
    SLIT_CENTER_V_MM=0.0,
    SHIFT_X_FLAG=0,
    SHIFT_X_VALUE=0.0,
    SHIFT_BETAX_FLAG=0,
    SHIFT_BETAX_VALUE=0.0,
    TRAJ_RESAMPLING_FACTOR=100.0,
    SLIT_POINTS_FACTOR =1,
    outFile="",
    outFileTraj="",
    )

#
# example plot
#
from srxraylib.plot.gol import plot
plot(energy, flux,
    xtitle="Photon energy [eV]", ytitle="Flux [photons/s/0.1%bw]", title="Wiggler Flux",
    xlog=True, ylog=True, grid=True, show=False)
plot(energy, spectral_power,
    xtitle="Photon energy [eV]", ytitle="Power [W/eV]", title="Wiggler Spectral Power",
    xlog=True, ylog=True, grid=True, show=False)
plot(energy, cumulated_power,
    xtitle="Photon energy [eV]", ytitle="Cumulated Power [W]", title="Wiggler Cumulated Power",
    xlog=False, ylog=False, grid=True, show=True)
#
# end script
#
