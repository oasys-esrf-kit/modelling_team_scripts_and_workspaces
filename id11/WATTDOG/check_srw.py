# -----------------------------------------------------------
# SRW script to compute undulator power and power density
# -----------------------------------------------------------

import numpy as np
try:
    from srwpy.srwlib import *
    from srwpy.uti_plot import *
except:
    from wofrysrw.util.srw import *
import numpy
# -----------------------------------------------------------
# SRW Harmonic 3 Emission Calculation
# -----------------------------------------------------------


# -----------------------------------------------------------
# Electron beam
# -----------------------------------------------------------

eBeam = SRWLPartBeam()
eBeam.Iavg = 0.2
eBeam.partStatMom1.gamma = 6.0 / 0.00051099895

# -----------------------------------------------------------
# Undulator
# -----------------------------------------------------------

K = 1.358
lam_u = 0.018
n = 3

gamma = eBeam.partStatMom1.gamma
hc = 1.239841984e-06  # eV*m

E1 = 2 * gamma**2 * hc / (lam_u * (1 + K**2/2))
En = n * E1

print("Resonant energy harmonic 3 [eV] =", En)

und = SRWLMagFldU([SRWLMagFldH(1, 'v', K, 0)], _per=lam_u, _nPer=111)

# -----------------------------------------------------------
# Energy Scan Around Harmonic 3
# -----------------------------------------------------------

dE = En * 0.01  # ±1% scan

mesh = SRWLRadMesh(
    En - dE, En + dE, 200,  # energy scan
    -0.0005, 0.0005, 51,    # horizontal angle [rad]
    -0.0005, 0.0005, 51,    # vertical angle [rad]
    30.0
)

stk = SRWLStokes()
stk.allocate(mesh.ne, mesh.nx, mesh.ny)
stk.mesh = mesh

# -----------------------------------------------------------
# Calculate Undulator Radiation (flux)
# -----------------------------------------------------------

arPrec = [1,  # method (undulator)
          1,  # precision
          3,  # harmonic start
          3,  # harmonic end
          1,  # longitudinal precision
          1,  # azimuthal precision
          1]  # flux

srwl.CalcStokesUR(stk, eBeam, und, arPrec)

# -----------------------------------------------------------
# Extract maximum flux
# -----------------------------------------------------------

flux_array = np.array(stk.arS).reshape((mesh.ne, mesh.ny, mesh.nx))

max_flux = flux_array.max()

print("Max Flux harmonic 3 [ph/s/0.1%bw] =", max_flux)