import os
import time
import numpy as np
from scipy.signal import savgol_filter
import xrt

from xrt.backends.raycing import BeamLine
from xrt.backends.raycing.sources import Undulator
from xrt.backends.raycing.oes import Plate
from xrt.backends.raycing.materials import Material
from xrt.backends.raycing.apertures import RectangularAperture
from xrt.backends.raycing.oes import DoubleParaboloidLens
from xrt.backends.raycing.screens import Screen
from xrt.plotter import XYCAxis, XYCPlot
from xrt.runner import run_ray_tracing

from id09_xrt import CSS_TO_U17, CSS_TO_M2V, CSS_TO_M2H, CSS_TO_W0
from id09_xrt import CSS_TO_PS, CSS_TO_W1, CSS_TO_W2, CSS_TO_S1
from id09_xrt import CSS_TO_M1, CSS_TO_S2, CSS_TO_BV1, CSS_TO_EH1W2
from id09_xrt import CSS_TO_EH2W1, CSS_TO_XSHUT, CSS_TO_HSC, CSS_TO_CRL
from id09_xrt import CSS_TO_SS2, CSS_TO_XEYE, CSS_TO_SAMPLE
from id09_xrt import M1_TO_S2, M1_TO_BV1, M1_TO_EH1W2
from id09_xrt import M1_TO_EH2W1, M1_TO_XSHUT, M1_TO_HSC, M1_TO_CRL
from id09_xrt import M1_TO_SS2, M1_TO_XEYE, M1_TO_SAMPLE, CRL_TO_SAMPLE
from id09_xrt import ToroidMirrorDistorted


force_beam = True
nrays = 20e3
repeats = 20

energy = 18.07e3  # 15.65e3/15.703e3, 15e3/15.05e3, 18e3/18.07e3

dE = 10*np.round(energy/1e3*4/3)
# eMin = energy - 2.5*dE
eMin = 17000
# eMax = energy + dE
eMax = 18500
phg = 0.5
pvg = 0.5
s1hg = 1.5
s1vg = 1.5
m1Rm = 4.042e6  # 2nd-source
m1incl = 3.37e-3  # 2nd-souce
s2hg = 2
s2vg = 2
lens_thickness = 0.025  # mm
lens_piling_thickness = 1.2  # mm

lens_radius = 0.1
lens_number = 17*2  # 13*2, 17*2

ss2hg = 0.5
ss2vg = 0.5

output_folder = f"res_{energy/1e3:.2f}keV_CRL_N{lens_number}"
output_folder += f"_R{lens_radius*1e3:.0f}"
os.makedirs(output_folder, exist_ok=True)


def build_beamline(force_beam=False):

    global REPETITION
    REPETITION = 0

    bl = BeamLine()
    bl.name = 'id09'
    bl.energy = energy
    bl.erange = eMin, eMax
    bl.force_beam = force_beam

    u17 = Undulator(bl, name='u17', center=[0, CSS_TO_U17, 0], period=17,
                    n=117, eE=6, eI=0.2, eEpsilonX=0.130, eEpsilonZ=0.010,
                    eEspread=9.4e-4, eSigmaX=30, eSigmaZ=5.2, distE='eV',
                    targetE=[energy, 1], eMin=eMin, eMax=eMax, nrays=nrays)

    u17_screen = Screen(bl, 'U17 screen', center=(0, CSS_TO_U17, 0))

    m2v = RectangularAperture(bl, name='m2v', center=(0, CSS_TO_M2V, 0),
                              kind=('left', 'right', 'bottom', 'top'),
                              opening=[-25/2, 25/2, -4/2, 4/2])

    m2h = RectangularAperture(bl, name='m2h', center=(0, CSS_TO_M2H, 0),
                              kind=('left', 'right', 'bottom', 'top'),
                              opening=[-4/2, 4/2, -25/2, 25/2])

    w0 = Plate(bl, name='fe_window', center=(0, CSS_TO_W0, 0), pitch=np.pi/2,
               material=Material('C', rho=3.52, kind='plate'), t=0.3,
               surface=r'diamond 300$\mu$m')

    w0_screen = Screen(bl, 'FE window screen', center=(0, CSS_TO_W0, 0))

    ps = RectangularAperture(bl, name='ps', center=(0, CSS_TO_PS, 0),
                             kind=('left', 'right', 'bottom', 'top'),
                             opening=[-phg/2, phg/2, -pvg/2, pvg/2])

    ps_screen = Screen(bl, 'PS screen', center=(0, CSS_TO_PS, 0))

    w1 = Plate(bl, name='w1', center=(0, CSS_TO_W1, 0), pitch=np.pi/2,
               material=Material('Be', rho=1.848, kind='plate'), t=0.35,
               surface=r'beryllium 350$\mu$m')

    w2 = Plate(bl, name='w2', center=(0, CSS_TO_W2, 0), pitch=np.pi/2,
               material=Material('Be', rho=1.848, kind='plate'), t=0.35,
               surface=r'beryllium 350$\mu$m')

    s1 = RectangularAperture(bl, name='s1', center=(0, CSS_TO_S1, 0),
                             kind=('left', 'right', 'bottom', 'top'),
                             opening=[-s1hg/2, s1hg/2, -s1vg/2, s1vg/2])

    s1_screen = Screen(bl, name='S1 screen', center=(0, CSS_TO_S1, 0))

    m1 = ToroidMirrorDistorted(distorsion_factor=1, bl=bl,
                               name='toroidal mirror',
                               center=(0, CSS_TO_M1, 0),
                               material=Material('Pd', rho=11.8,
                                                 kind='mirror'),
                               R=m1Rm, r=46, pitch=m1incl, yaw=0,
                               limPhysX=[-5, 5], limPhysY=[-300, 300])

    m1_screen = Screen(bl, name='M1 screen', center=(0, CSS_TO_M1, 0))

    s2 = RectangularAperture(bl, name='s2',
                             center=(0, CSS_TO_S2, 2*m1.pitch*M1_TO_S2),
                             kind=('left', 'right', 'bottom', 'top'),
                             opening=[-s2hg/2, s2hg/2, -s2vg/2, s2vg/2])

    s2_screen = Screen(bl, name='S1 screen',
                       center=(0, CSS_TO_S2, 2*m1.pitch*M1_TO_S2))

    bv1_screen = Screen(bl, name='BV1 screen',
                        center=(0, CSS_TO_BV1, 2*m1.pitch*M1_TO_BV1))

    eh1w2_screen = Screen(bl, 'EH2 exit wall',
                          center=(0, CSS_TO_EH1W2, 2*m1.pitch*M1_TO_EH1W2))

    eh2w1_screen = Screen(bl, 'EH1 entry wall',
                          center=(0, CSS_TO_EH2W1, 2*m1.pitch*M1_TO_EH2W1))

    xshut_screen = Screen(bl, name='XSHUT screen',
                          center=(0, CSS_TO_XSHUT, 2*m1.pitch*M1_TO_XSHUT))

    hsc_screen1 = Screen(bl, name='highspeed chopper screen 1',
                         center=(0, CSS_TO_HSC, 2*m1.pitch*M1_TO_HSC))

    hsc = RectangularAperture(bl, name='highspeed chopper',
                              center=(0, CSS_TO_HSC, 2*m1.pitch*M1_TO_HSC),
                              kind=('left', 'right', 'bottom', 'top'),
                              opening=[-1, 1, -0.130/2, +0.130/2])

    hsc_screen2 = Screen(bl, name='highspeed chopper screen 2',
                         center=(0, CSS_TO_HSC, 2*m1.pitch*M1_TO_HSC))

    crl = DoubleParaboloidLens(bl, name='CRL',
                               center=[0, CSS_TO_CRL, 2*m1.pitch*M1_TO_CRL],
                               pitch=np.pi/2+2*m1.pitch,
                               material=Material('C', rho=3.52, kind='lens'),
                               t=lens_thickness,
                               focus=lens_radius/2,  # lens focus f=R/2
                               nCRL=lens_number,
                               zmax=(lens_piling_thickness-lens_thickness)/2
                               )

    crl_screen = Screen(bl, name='CRL screen',
                        center=(0, CSS_TO_CRL, 2*m1.pitch*M1_TO_CRL))

    ss2 = RectangularAperture(bl, name='ss2',
                              center=(0, CSS_TO_SS2, 2*m1.pitch*M1_TO_SS2),
                              kind=('left', 'right', 'bottom', 'top'),
                              opening=[-ss2hg/2, ss2hg/2, -ss2vg/2, ss2vg/2])

    ss2_screen = Screen(bl, name='SS2 screen',
                        center=(0, CSS_TO_SS2, 2*m1.pitch*M1_TO_SS2))

    xeye_screen = Screen(bl, name='XEYE screen',
                         center=(0, CSS_TO_XEYE, 2*m1.pitch*M1_TO_XEYE))

    sample_screen = Screen(bl, name='sample screen',
                           center=[0, CSS_TO_SAMPLE,
                                   2*m1.pitch*(M1_TO_SAMPLE)])

    bl.source = u17
    bl.u17_screen = u17_screen
    bl.m2h = m2h
    bl.m2v = m2v
    bl.w0 = w0
    bl.w0_screen = w0_screen
    bl.ps = ps
    bl.ps_screen = ps_screen
    bl.w1 = w1
    bl.w2 = w2
    bl.s1 = s1
    bl.s1_screen = s1_screen
    bl.m1 = m1
    bl.m1_screen = m1_screen
    bl.s2 = s2
    bl.s2_screen = s2_screen
    bl.bv1_screen = bv1_screen
    bl.eh1w2_screen = eh1w2_screen
    bl.eh2w1_screen = eh2w1_screen
    bl.xshut_screen = xshut_screen
    bl.hsc_screen1 = hsc_screen1
    bl.hsc = hsc
    bl.hsc_screen2 = hsc_screen2
    bl.crl = crl
    bl.crl_screen = crl_screen
    bl.ss2 = ss2
    bl.ss2_screen = ss2_screen
    bl.xeye_screen = xeye_screen
    bl.sample_screen = sample_screen

    return bl


def run_process(bl):

    global REPETITION

    print(f"\nREPETITION = {REPETITION}")
    t0 = time.time()
    folder = f"./beam_{energy*1e-3:.2f}keV_crl/"
    os.makedirs(folder, exist_ok=True)
    fname = folder + f"beam{REPETITION:02d}.npy"
    if os.path.isfile(fname) and not bl.force_beam:
        print(f"Reading beam from file: {fname}")
        beam = np.load(fname, allow_pickle=True).item()
    else:
        beam = bl.source.shine()
        np.save(fname, beam)
    dt = time.time() - t0
    print(f"Time needed to build rays {dt:.3f} sec")

    _ = bl.m2v.propagate(beam)
    _ = bl.m2h.propagate(beam)
    u17_screen = bl.u17_screen.expose(beam)
    beam, _, _ = bl.w0.double_refract(beam)
    w0_screen = bl.w0_screen.expose(beam)
    ps_screen = bl.ps_screen.expose(beam)
    _ = bl.ps.propagate(beam)
    s1_screen = bl.s1_screen.expose(beam)
    m1_screen = bl.m1_screen.expose(beam)
    beam, _ = bl.m1.reflect(beam)
    s2_screen = bl.s2_screen.expose(beam)
    bv1_screen = bl.bv1_screen.expose(beam)
    eh1w2_screen = bl.eh1w2_screen.expose(beam)
    eh2w1_screen = bl.eh2w1_screen.expose(beam)
    xshut_screen = bl.xshut_screen.expose(beam)
    hsc_screen1 = bl.hsc_screen1.expose(beam)
    _ = bl.hsc.propagate(beam)
    hsc_screen2 = bl.hsc_screen2.expose(beam)
    crl_screen = bl.crl_screen.expose(beam)
    beam, _, _ = bl.crl.multiple_refract(beam)  # 2nd source
    _ = bl.ss2.propagate(beam)
    xeye_screen = bl.xeye_screen.expose(beam)
    sample_screen = bl.sample_screen.expose(beam)

    REPETITION += 1

    out = dict(u17_screen=u17_screen,
               w0_screen=w0_screen,
               ps_screen=ps_screen,
               s1_screen=s1_screen,
               m1_screen=m1_screen,
               s2_screen=s2_screen,
               bv1_screen=bv1_screen,
               eh1w2_screen=eh1w2_screen,
               eh2w1_screen=eh2w1_screen,
               xshut_screen=xshut_screen,
               hsc_screen1=hsc_screen1,
               hsc_screen2=hsc_screen2,
               crl_screen=crl_screen,
               xeye_screen=xeye_screen,
               sample_screen=sample_screen)

    return out


xrt.backends.raycing.run.run_process = run_process


def make_plot(bl, screen, size=100, bins=1024, cbins=256):

    global output_folder

    xaxis = XYCAxis(label='x', limits=[-size/2, size/2], unit='um', bins=bins,
                    ppb=int(1024/bins))
    yaxis = XYCAxis(label='z', limits=[-size/2, size/2], unit='um', bins=bins,
                    ppb=int(1024/bins))
    caxis = XYCAxis(label='energy', limits=[bl.source.eMin, bl.source.eMax],
                    unit='eV', bins=cbins, fwhmFormatStr="%.2f",
                    ppb=int(512/cbins))

    fname = output_folder + f"/screen_{screen}.png"

    plot = XYCPlot(beam=screen, xaxis=xaxis, yaxis=yaxis, caxis=caxis,
                   saveName=fname)

    return plot


def rgb2gray(rgb):
    r, g, b = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
    return gray


def crl_info(bl):

    global energy

    N = bl.crl.nCRL
    R = bl.crl.focus*2
    delta = (1 - bl.crl.material.get_refractive_index(energy).real)
    d = bl.crl.t  # wall thickness
    zmax = bl.crl.zmax  # d_depth
    p = 2*zmax + d  # lens thickness
    A = 2*np.sqrt(2*R*zmax)  # lens aperture (Thomas Roth)

    mu = 100*bl.crl.material.get_absorption_coefficient(energy)  # 1/m
    # A = 2*np.sqrt(2*R/1e3/mu)*1e3  # double paraboloid aperture (mm)
    # p = A**2/(4*R) + d  # double paraboloid thickness (mm)

    L = N*p  # crl thickness

    out_str = [
        f"crl position = {bl.crl.center}",
        f"crl-to-sample distance = {CRL_TO_SAMPLE*1e-3:g} m",
        f"energy (fundamental) = {energy*1e-3:.3f} keV",
        f"material mu = {mu:g} 1/m",
        f"material delta = {delta:.4e}",
        f"number of lenses = {N}",
        f"lens ideal radius of curvature = {2*CRL_TO_SAMPLE*delta*N:.3f} mm",
        f"lens chosen radius of curvature = {R:.3f} mm",
        f"crl focal length (with chosen R) = {R/(2*delta*N)*1e-3:.3f} m",
        f"lens wall thickness = {d*1e3:g} um",
        f"lens zmax (d_depth) = {zmax*1e3:g} um",
        f"lens piling thickness = {p:g} mm",
        f"lens aperture (D_phys) = {A*1e3:.1f} um",
        f"crl total thickness = {L:g} mm",
    ]
    print()
    out_str = "\n".join(out_str) + "\n"
    print(out_str)
    return out_str


def do_after_script(plots, d, bl):

    global output_folder

    out_str1 = crl_info(bl)
    d = np.array(d)/1e3
    fwhm_h = np.array([plot.dx for plot in plots])
    fwhm_v = np.array([plot.dy for plot in plots])
    ecen = np.array([plot.cE/1e3 for plot in plots])
    bw = np.array([plot.dE/energy*1e2 for plot in plots])
    flux = np.array([plot.flux for plot in plots])
    profile_h = np.array([[p.xaxis.binCenters, p.xaxis.total1D]
                          for p in plots])
    profile_v = np.array([[p.yaxis.binCenters, p.yaxis.total1D]
                          for p in plots])
    profile_e = np.array([[p.caxis.binCenters, p.caxis.total1D]
                          for p in plots])
    # img_h = [rgb2gray(p.ax1dHistX.images[0].get_array()[::-1, :, :])
    #          for p in plots]
    # img_v = [rgb2gray(p.ax1dHistY.images[0].get_array()) for p in plots]
    # img_e = [rgb2gray(p.ax1dHistE.images[0].get_array()) for p in plots]

    epeak = []
    for k, plot in enumerate(plots):
        fname = output_folder + f"/res_{plot.beam}_profile_h.txt"
        np.savetxt(fname, profile_h[k].T)
        fname = output_folder + f"/res_{plot.beam}_profile_v.txt"
        np.savetxt(fname, profile_v[k].T)
        fname = output_folder + f"/res_{plot.beam}_profile_e.txt"
        np.savetxt(fname, profile_e[k].T)
        e, y = np.loadtxt(fname, unpack=True)
        epeak.append(e[np.argmax(savgol_filter(y, 51, 3))])
    epeak = np.array(epeak)/1e3

    out_str2 = [
        f"epeak = {epeak}",
        f"ecen = {ecen}",
        f"bw = {bw}",
        f"fwhm_h = {fwhm_h}",
        f"fwhm_v = {fwhm_v}",
        f"flux = {flux}",
        f"CRL transmission = {flux[-2]/flux[-1]*100:.1f} per cent"
    ]

    out_str2 = '\n'.join(out_str2)
    print(out_str2)
    out_str1 += '\n' + out_str2 + '\n'
    with open(output_folder + "/results.txt", 'w') as f:
        f.write(out_str1)


def main():
    bl = build_beamline(force_beam=force_beam)
    screens = ['hsc_screen1', 'crl_screen', 'sample_screen']
    sizes = [500, 750, 100]
    d = [CSS_TO_HSC, CSS_TO_CRL, CSS_TO_SAMPLE]
    n_screens_to_keep = 0
    screens = screens[-n_screens_to_keep:]
    sizes = sizes[-n_screens_to_keep:]
    d = d[-n_screens_to_keep:]
    plots = [make_plot(bl, screen, size) for screen, size
             in zip(screens[::-1], sizes[::-1])]

    run_ray_tracing(plots=plots, repeats=repeats, pickleEvery=1,
                    updateEvery=1, beamLine=bl, afterScript=do_after_script,
                    afterScriptArgs=[plots, d, bl])


if __name__ == '__main__':
    main()
