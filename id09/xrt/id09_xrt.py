import os
import time
import scipy
import numpy as np
import matplotlib.pyplot as plt
import xrt

from matplotlib.ticker import FormatStrFormatter
from xrt.backends.raycing import BeamLine
from xrt.backends.raycing.sources import Undulator
from xrt.backends.raycing.oes import Plate
from xrt.backends.raycing.materials import Material, CrystalSi, Multilayer
from xrt.backends.raycing.apertures import RectangularAperture
from xrt.backends.raycing.apertures import PolygonalAperture
from xrt.backends.raycing.oes import DCM
from xrt.backends.raycing.oes import ToroidMirror
from xrt.backends.raycing.oes import DoubleParaboloidLens
from xrt.backends.raycing.screens import Screen
# from xrt.plotter import XYCAxis, XYCPlot
# from xrt.runner import run_ray_tracing


Z0 = 1400  # ASSUMED BEAM HEIGHT WITH RESPECT TO FLOOR

CSS_TO_U23 = -1250
CSS_TO_U17 = 1250
CSS_TO_M2V = 22225
CSS_TO_M2H = 23211
CSS_TO_W0 = 23428

CSS_TO_OHW1 = 26039
CSS_TO_PS = 27100
CSS_TO_HPDA = 27800
CSS_TO_W1 = 28500
CSS_TO_PIC = 28800
CSS_TO_HLC = 30000
CSS_TO_W2 = 30430
CSS_TO_MONO1 = 31430
CSS_TO_S1 = 37500
CSS_TO_OHW2 = 38950

CSS_TO_EH1W1 = 39100
CSS_TO_M1 = 44540
CSS_TO_S2 = 45800
CSS_TO_BV1 = 46300
CSS_TO_BV2 = 46540
CSS_TO_EH1W2 = 50000
CSS_TO_SRC2 = CSS_TO_M1 + 8050
CSS_TO_EH2W1 = 53505  # (-2.784 m from sample)
CSS_TO_XSHUT = 54195 # (-2.094 m from sample) 
CSS_TO_HSC = 54535  # (-1.754 m from sample)
CSS_TO_ATT = 54875  # (-1.414 m from sample)
CSS_TO_SS1 = 54950  # (-1.339 m from sample)
CSS_TO_CRL = 55439  # (-85 cm from sample)
CSS_TO_SS2 = 55641  # (-64.8 cm from sample)
CSS_TO_W3 = 55691  # (-59.8 cm from sample)
CSS_TO_XEYE = 55825  # (-46.4 cm from sample
CSS_TO_SAMPLE = 56289
CSS_TO_EH2W2 = 62500

M1_TO_S2 = CSS_TO_S2 - CSS_TO_M1
M1_TO_BV1 = CSS_TO_BV1 - CSS_TO_M1
M1_TO_BV2 = CSS_TO_BV2 - CSS_TO_M1
M1_TO_EH1W2 = CSS_TO_EH1W2 - CSS_TO_M1
M1_TO_SRC2 = CSS_TO_SRC2 - CSS_TO_M1
M1_TO_EH2W1 = CSS_TO_EH2W1 - CSS_TO_M1
M1_TO_XSHUT = CSS_TO_XSHUT - CSS_TO_M1
M1_TO_HSC = CSS_TO_HSC - CSS_TO_M1
M1_TO_CRL = CSS_TO_CRL - CSS_TO_M1
M1_TO_SS2 = CSS_TO_SS2 - CSS_TO_M1
M1_TO_XEYE = CSS_TO_XEYE - CSS_TO_M1
M1_TO_SAMPLE = CSS_TO_SAMPLE - CSS_TO_M1
HSC_TO_SAMPLE = CSS_TO_SAMPLE - CSS_TO_HSC
CRL_TO_SAMPLE = CSS_TO_SAMPLE - CSS_TO_CRL

U17_TO_M2V = CSS_TO_M2V - CSS_TO_U17
U17_TO_M2H = CSS_TO_M2H - CSS_TO_U17
U17_TO_W0 = CSS_TO_W0 - CSS_TO_U17
U17_TO_PS = CSS_TO_PS - CSS_TO_U17
U17_TO_W1 = CSS_TO_W1 - CSS_TO_U17
U17_TO_W2 = CSS_TO_W2 - CSS_TO_U17
U17_TO_S1 = CSS_TO_S1 - CSS_TO_U17
U17_TO_M1 = CSS_TO_M1 - CSS_TO_U17
U17_TO_S2 = CSS_TO_S2 - CSS_TO_U17
U17_TO_BV1 = CSS_TO_BV1 - CSS_TO_U17
U17_TO_BV2 = CSS_TO_BV2 - CSS_TO_U17
U17_TO_XSHUT = CSS_TO_XSHUT - CSS_TO_U17
U17_TO_HSC = CSS_TO_HSC - CSS_TO_U17
U17_TO_CRL = CSS_TO_CRL - CSS_TO_U17
U17_TO_SS2 = CSS_TO_SS2 - CSS_TO_U17
U17_TO_XEYE = CSS_TO_XEYE - CSS_TO_U17
U17_TO_SAMPLE = CSS_TO_SAMPLE - CSS_TO_U17

U23_TO_PS = CSS_TO_PS - CSS_TO_U23
U23_TO_M1 = CSS_TO_M1 - CSS_TO_U23
U23_TO_BV1 = CSS_TO_BV1 - CSS_TO_U23
U23_TO_BV2 = CSS_TO_BV2 - CSS_TO_U23
U23_TO_XSHUT = CSS_TO_XSHUT - CSS_TO_U23
U23_TO_HSC = CSS_TO_HSC - CSS_TO_U23
U23_TO_CRL = CSS_TO_CRL - CSS_TO_U23
U23_TO_XEYE = CSS_TO_XEYE - CSS_TO_U23
U23_TO_SAMPLE = CSS_TO_SAMPLE - CSS_TO_U23


def get_theta_psi(slit_width, slit_height, slit_distance, slit_npt=51):
    theta = np.arctan(0.5*slit_width/slit_distance)
    theta = np.linspace(-theta, theta, slit_npt)
    psi = np.arctan(0.5*slit_height/slit_distance)
    psi = np.linspace(-psi, psi, slit_npt)
    return theta, psi


def get_power_through_slit(und, K=0.8525, e_min=12e3, e_max=16e3, e_npt=201,
                           slit_width=0.5, slit_height=0.5,
                           slit_distance=U17_TO_PS, slit_npt=21,
                           harmonic=None, shine=False, return_res=False,
                           return_spectrum=True):

    EV2ERG = 1.602176565e-12

    und.K = K
    und.eMin = e_min
    und.eMax = e_max
    if shine:
        und.shine()

    energy = np.linspace(e_min, e_max, e_npt)

    theta, psi = get_theta_psi(slit_width, slit_height, slit_distance,
                               slit_npt)

    dE = energy[1] - energy[0]
    dtheta = theta[1] - theta[0]
    dpsi = psi[1] - psi[0]

    print("\nParameters")
    print("----------")
    print(f"Undulator: {und.name}")
    print(f"K: {und.K}")
    print(f"distE: {und.distE}")
    print(f"e_min: {und.eMin} eV")
    print(f"e_max: {und.eMax} eV")
    print(f"slit_width: {slit_width} mm")
    print(f"slit_height: {slit_height} mm")
    print(f"slit_distance: {slit_distance*1e-3} m")
    print(f"theta range: ({theta[0]*1e6:.2f}, {theta[-1]*1e6:.2f}) urad")
    print(f"psi range: ({psi[0]*1e6:.2f}, {psi[-1]*1e6:.2f}) urad")
    print(f"dE: {dE} eV")
    print(f"dtheta: {dtheta*1e6:.2f} urad")
    print(f"dpsi: {dpsi*1e6:.2f} urad")

    t0 = time.time()
    print("\nCalculating intensities on mesh...")
    I0 = und.intensities_on_mesh(energy, theta, psi, harmonic=harmonic)[0]
    print(f"Calculation took {time.time()-t0:g} sec")

    if und.distE == 'BW':
        I0 *= 1e3
    else:  # 'eV'
        I0 *= energy[:, np.newaxis, np.newaxis]

    power = I0.sum() * dtheta * dpsi * dE * EV2ERG * 1e-7

    print(f"\npower: {power} W")

    if return_res:
        return energy, theta, psi, I0

    if return_spectrum:
        flux = I0.sum(axis=-1).sum(axis=-1) * dtheta * dpsi
        flux *= 1e-3  # ph/s/0.1%bw
        spectral_power = flux * EV2ERG * 1e-7  # W/eV
        return energy, flux, spectral_power

    return power


def get_u17_theta_psi(slit_width=0.5, slit_height=0.5,
                      slit_position=U17_TO_PS):
    theta, psi = get_theta_psi(slit_width, slit_height, slit_position)
    return theta, psi


def get_u23_theta_psi(slit_width=0.5, slit_height=0.5,
                      slit_position=U23_TO_PS):
    theta, psi = get_theta_psi(slit_width, slit_height, slit_position)
    return theta, psi


def get_u17_K(gap, Kmax=0.8525, gap_min=6, lu=17, b1=np.pi, b2=0):
    K = Kmax*np.exp(-b1*(gap-gap_min)/lu)*np.exp(-b2*((gap-gap_min)/lu)**2)
    return K


def get_u23_K(gap, Kmax=0.751, gap_min=6, lu=17, b1=np.pi, b2=0):
    K = Kmax*np.exp(-b1*(gap-gap_min)/lu)*np.exp(-b2*((gap-gap_min)/lu)**2)
    return K


def get_u17_flux_through_aperture(slit_width=0.5, slit_height=0.5,
                                  slit_position=CSS_TO_PS, emin=5e3, emax=30e3,
                                  en=1000, gap=6):

    energy = np.linspace(emin, emax, en)

    theta = np.arctan(slit_width/2/slit_position)
    psi = np.arctan(slit_height/2/slit_position)
    theta = np.linspace(-theta, theta, 51)
    psi = np.linspace(-psi, psi, 51)

    bl = build_beamline()
    u17 = bl.u17
    u17.K = get_u17_K(gap=gap)
    u17.eMin = emin
    u17.eMax = emax

    I0, l1, l2, l3 = u17.intensities_on_mesh(energy, theta, psi)
    dtheta, dpsi = theta[1] - theta[0], psi[1] - psi[0]
    flux = I0.sum(axis=(1, 2)) * dtheta * dpsi

    return energy*1e-3, flux


def get_u17_power_through_aperture(slit_width=0.5, slit_height=0.5,
                                   slit_position=CSS_TO_PS, emin=5e3,
                                   emax=30e3, en=1001, gap=6):

    energy = np.linspace(emin, emax, en)

    theta = np.arctan(slit_width/2/slit_position)
    theta = np.linspace(-theta, theta, 21)
    psi = np.arctan(slit_height/2/slit_position)
    psi = np.linspace(-psi, psi, 21)

    gap = np.linspace(6, 28, 20)
    K = get_u17_K(gap)
    harmonic = range(1, 9)

    u17 = Undulator(name='u17', center=[0, CSS_TO_U17, 0], period=17, n=117,
                    eE=6, eI=0.2, eEpsilonX=0.130, eEpsilonZ=0.010,
                    eEspread=9.4e-4, eSigmaX=30, eSigmaZ=5.2, distE='BW',
                    K=0.8525)

    power = u17.power_vs_K(energy, theta, psi, harmonic, K)

    plt.plot(gap, power, 'o-')
    ax = plt.gca()
    ax.set_xlabel(u'gap (mm)')
    ax.set_ylabel(u'power through 0.5x0.5mm2 primary slits (W)')

    maxPower = power.max()
    for i, (pi, Ki, gi) in enumerate(zip(power, K, gap)):
        ax.text(gi, pi + maxPower*0.03 * (i % 2 * 2 - 1),
                '{0}{1:.1f}'.format('K=' if i == 0 else '', Ki), fontsize=9,
                ha='center', va='center')

    # plt.savefig("id09_power_u17.png")
    plt.show()


def get_u17_tuning_curve(slit_width=0.5, slit_height=0.5,
                         slit_position=CSS_TO_PS, emin=5e3, emax=55e3,
                         en=1000):

    energy = np.linspace(emin, emax, en)
    bl = build_beamline()
    u17 = bl.u17
    u17.eMin = emin
    u17.eMax = emax

    theta = 3*np.arctan(slit_width/2/slit_position)
    psi = 3*np.arctan(slit_height/2/slit_position)
    print(f"Calculating flux through {theta*1e6:.2e} x {psi*1e6:.2e} urad^2")
    theta = np.linspace(-theta, theta, 11)
    psi = np.linspace(-psi, psi, 11)

    xlims, ylims = (10, 40), (1e12, 1e16)
    gaps = np.array([6, 7.5, 9, 10.5, 12, 14])
    Ks = get_u17_K(gaps)
    harmonics = [1, 2]
    colors = ['b', 'g']

    tunesE, tunesF = u17.tuning_curves(energy, theta, psi, harmonics, Ks)

    for tuneE, tuneF, harmonic, color in zip(tunesE, tunesF, harmonics,
                                             colors):
        plt.loglog(tuneE, tuneF, 'o-', label='{0}'.format(harmonic),
                   color=color)

    ax = plt.gca()
    ax.set_xlabel(u'Photon energy (keV)')
    ax.set_ylabel(u'Flux through primary slits (ph/s/0.1% bw)')
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    # ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    plt.grid(which='both')
    ax.legend(loc='upper right', title='harmonics')

    for tuneE, tuneF, harmonic in zip(tunesE, tunesF, harmonics):
        for x, y, K, gap in zip(tuneE, tuneF, Ks, gaps):
            if xlims[0] < x < xlims[1] and ylims[0] < y < ylims[1]:
                ax.text(x*1.02, y*1.05, '{0:g} mm'.format(gap), fontsize=9)


def get_mono_reflectivity(emin=5e3, emax=55e3, en=1000):
    energy = np.linspace(emin, emax, en)
    bl = build_beamline()
    mono = bl.mono
    theta = mono.get_Bragg_angle(energy)
    curS, curP = mono.get_amplitude(energy, np.sin(theta))
    plt.plot(energy, abs(curS)**2)
    plt.show()


def get_m1_slope_error(x=[-5, 5], nx=51, distorsion_factor=1,
                       fname="ID9-toreSESO-D2-P2.slp"):
    if fname == "2018.05.03_id09_slope_error.txt":
        y_mm, dz_dy_mrad = np.loadtxt(fname, unpack=True)
        dz_dy_rad = dz_dy_mrad*1e-3
    elif fname == "ID9-toreSESO-D2-P2.slp":
        y_mm, dz_dy_urad = np.loadtxt(fname, skiprows=6, unpack=True)
        poly = np.polyfit(y_mm, dz_dy_urad, 1)
        dz_dy_urad = dz_dy_urad - np.polyval(poly, y_mm)
        dz_dy_rad = dz_dy_urad*1e-6
    elif fname == "dabam-030.dat":
        y_mm, dz_dy_urad = np.loadtxt(fname, unpack=True)
        dz_dy_rad = dz_dy_urad*1e-6
    else:
        raise ValueError(f"'fname' = {fname} not recognized")
    dz_dy_rad *= distorsion_factor
    y_mm -= np.mean(y_mm)
    dz_dy_rad -= np.mean(dz_dy_rad)
    z_mm = scipy.integrate.cumtrapz(dz_dy_rad, x=y_mm, initial=0)
    x_mm = np.linspace(x[0], x[1], nx)
    z_mm = z_mm[:, np.newaxis]
    z_mm = np.repeat(z_mm, nx, axis=1).T
    return x_mm, y_mm, z_mm


class ToroidMirrorDistorted(ToroidMirror):

    def __init__(self, distorsion_factor=1, *args, **kwargs):
        ToroidMirror.__init__(self, *args, **kwargs)
        x = self.limPhysX
        x_dist, y_dist, z_dist = get_m1_slope_error(
            x=x, distorsion_factor=distorsion_factor)
        self.n_x_dist = len(x_dist)
        self.n_y_dist = len(y_dist)
        self.limPhysX = np.min(x_dist), np.max(x_dist)
        self.limPhysY = np.min(y_dist), np.max(y_dist)
        self.get_surface_limits()
        self.x_grad, self.y_grad = np.gradient(z_dist, x_dist, y_dist)
        self.x_grad = np.arctan(self.x_grad)
        self.y_grad = np.arctan(self.y_grad)
        self.z_spline = scipy.ndimage.spline_filter(z_dist)
        self.x_grad_spline = scipy.ndimage.spline_filter(self.x_grad)
        self.y_grad_spline = scipy.ndimage.spline_filter(self.y_grad)

    def local_z_distorted(self, x, y):
        coords = np.array(
            [(x-self.limPhysX[0]) /
             (self.limPhysX[1]-self.limPhysX[0]) * (self.n_x_dist-1),
             (y-self.limPhysY[0]) /
             (self.limPhysY[1]-self.limPhysY[0]) * (self.n_y_dist-1)])
        z = scipy.ndimage.map_coordinates(self.z_spline, coords,
                                          prefilter=True)
        return z

    def local_n_distorted(self, x, y):
        coords = np.array(
            [(x-self.limPhysX[0]) /
             (self.limPhysX[1]-self.limPhysX[0]) * (self.n_x_dist-1),
             (y-self.limPhysY[0]) /
             (self.limPhysY[1]-self.limPhysY[0]) * (self.n_y_dist-1)])
        a = scipy.ndimage.map_coordinates(self.x_grad_spline, coords,
                                          prefilter=True)
        b = scipy.ndimage.map_coordinates(self.y_grad_spline, coords,
                                          prefilter=True)
        return b, -a


def get_mRu_reflectivity(e=15e3, en=400, theta=0.615):
    emin = e-5e3
    emax = e+5e3
    energy = np.linspace(emin, emax, en)
    theta = np.deg2rad(theta)
    bl = build_beamline()
    mlRu = bl.mlRu
    # theta = mlRu.get_Bragg_angle(e)
    rs, rp = mlRu.get_amplitude(energy, np.sin(theta))[0:2]
    plt.plot(energy*1e-3, abs(rs)**2)
    plt.xlabel("Photon energy (keV)")
    plt.ylabel("Reflectivity")
    plt.grid()


def get_mlW_reflectivity(e=25e3, en=400, theta=0.558):
    emin = e-5e3
    emax = e+5e3
    energy = np.linspace(emin, emax, en)
    theta = np.deg2rad(theta)
    bl = build_beamline()
    mlW = bl.mlW
    # theta = mlW.get_Bragg_angle(e)
    rs, rp = mlW.get_amplitude(energy, np.sin(theta))[0:2]
    plt.plot(energy*1e-3, abs(rs)**2)
    plt.xlabel("Photon energy (keV)")
    plt.ylabel("Reflectivity")
    plt.grid()


def build_beamline(nrays=20e3, force_beam=False, phg=0.5, pvg=0.5, s1hg=1.5,
                   s1vg=1.5, s2hg=2, s2vg=2, m1Rm=7.5e6, m1rs=46,
                   m1pitch=2.495e-3, m1distorsion=1, xshuty=-0.15,
                   ss1hg=0.5, ss1vg=0.5, ss2hg=0.3, ss2vg=0.3):

    bl = BeamLine()
    bl.name = 'id09'
    bl.force_beam = force_beam

    u23 = Undulator(bl, name='u23', center=[0, CSS_TO_U23, Z0], period=23,
                    n=86, eE=6, eI=0.2, eEpsilonX=0.130, eEpsilonZ=0.010,
                    eEspread=9.4e-4, eSigmaX=30, eSigmaZ=5.2, distE='eV',
                    K=1.613, nrays=nrays)

    u23_screen = Screen(bl, 'u23 screen', center=[0, CSS_TO_U23, Z0])
    u23_screen.size = 0.2  # mm

    u17 = Undulator(bl, name='u17', center=[0, CSS_TO_U17, Z0], period=17,
                    n=117, eE=6, eI=0.2, eEpsilonX=0.130, eEpsilonZ=0.010,
                    eEspread=9.4e-4, eSigmaX=30, eSigmaZ=5.2, distE='eV',
                    K=0.8525, nrays=nrays)

    u17_screen = Screen(bl, 'u17 screen', center=[0, CSS_TO_U17, Z0])
    u17_screen.size = 0.2  # mm

    m2v = RectangularAperture(bl, name='m2v', center=[0, CSS_TO_M2V, Z0],
                              kind=('left', 'right', 'bottom', 'top'),
                              opening=[-25/2, 25/2, -4/2, 4/2])

    m2v_screen = Screen(bl, 'm2v screen', center=[0, CSS_TO_M2V, Z0])
    m2v_screen.size = 2  # mm

    m2h = RectangularAperture(bl, name='m2h', center=[0, CSS_TO_M2H, Z0],
                              kind=('left', 'right', 'bottom', 'top'),
                              opening=[-4/2, 4/2, -25/2, 25/2])

    m2h_screen = Screen(bl, 'm2h screen', center=[0, CSS_TO_M2H, Z0])
    m2h_screen.size = 2  # mm

    w0 = Plate(bl, name='w0', center=[0, CSS_TO_W0, Z0], pitch=np.pi/2,
               material=Material('C', rho=3.52, kind='plate'), t=0.3)
    w0_screen = Screen(bl, 'w0 screen', center=[0, CSS_TO_W0, Z0])
    w0_screen.size = 2  # mm

    ps = RectangularAperture(bl, name='ps', center=[0, CSS_TO_PS, Z0],
                             kind=('left', 'right', 'bottom', 'top'),
                             opening=[-phg/2, phg/2, -pvg/2, pvg/2])

    ps_screen = Screen(bl, 'ps screen', center=[0, CSS_TO_PS, Z0])
    ps_screen.size = phg + pvg

    w1 = Plate(bl, name='w1', center=[0, CSS_TO_W1, Z0], pitch=np.pi/2,
               material=Material('Be', rho=1.848, kind='plate'), t=0.35)
    w1_screen = Screen(bl, 'w1 screen', center=[0, CSS_TO_W1, Z0])
    w1_screen.size = phg + pvg

    w2 = Plate(bl, name='w2', center=[0, CSS_TO_W2, Z0], pitch=np.pi/2,
               material=Material('Be', rho=1.848, kind='plate'), t=0.35)
    w2_screen = Screen(bl, 'w2 screen', center=[0, CSS_TO_W2, Z0])
    w2_screen.size = phg + pvg

    mono1 = DCM(bl, name='mono1', center=[0, CSS_TO_MONO1, Z0],
                material=CrystalSi(hkl=(1, 1, 1), tK=273.15-186),
                alpha=np.radians(0),
                limPhysX=(-10, 10), limPhysY=(-30, 30),
                cryst2perpTransl=20, cryst2longTransl=65,
                limPhysX2=(-10, 10), limPhysY2=(-90, 90))

    s1 = RectangularAperture(bl, name='s1', center=[0, CSS_TO_S1, 'auto'],
                             kind=('left', 'right', 'bottom', 'top'),
                             opening=[-s1hg/2, s1hg/2, -s1vg/2, s1vg/2])

    s1_screen = Screen(bl, 's1 screen', center=[0, CSS_TO_S1, 'auto'])
    s1_screen.size = 1

    ohw2_screen = Screen(bl, 'ohw2 screen',
                         center=[0, CSS_TO_OHW2, 'auto'])
    ohw2_screen.size = 1  # mm
    
    eh1w1_screen = Screen(bl, 'eh1w1 screen',
                         center=[0, CSS_TO_EH1W1, 'auto'])
    eh1w1_screen.size = 1  # mm
    
    m1 = ToroidMirrorDistorted(distorsion_factor=m1distorsion, bl=bl,
                               name='m1', center=[0, CSS_TO_M1, 'auto'],
                               material=Material('Pd', rho=11.8,
                                                 kind='mirror'),
                               R=m1Rm, r=m1rs, pitch=m1pitch, yaw=0,
                               limPhysX=(-5, 5), limPhysY=(-300, 300))
    
    m1_screen = Screen(bl, 'm1 screen', center=[0, CSS_TO_M1, 'auto'])
    m1_screen.size = 1  # mm

    bv1 = Plate(bl, name='bv1', center=[0, CSS_TO_BV1, 'auto'],
                material=Material('C', rho=3.52, kind='plate'), t=0.3,
                pitch=np.pi/4)
    
    bv1_screen = Screen(bl, 'bv1 screen', center=[0, CSS_TO_BV1, 'auto'])
    bv1_screen.size = 1  # mm

    bv2 = Plate(bl, name='bv2', center=[0, CSS_TO_BV2, 'auto'],
                material=Material(('Y', 'Al', 'O'), (3, 5, 12),
                                  rho=4.57, kind='plate'),
                t=0.1,  # ??
                pitch=np.pi/4)
    
    bv2_screen = Screen(bl, 'bv2 screen', center=[0, CSS_TO_BV2, 'auto'])
    bv2_screen.size = 1  # mm

    s2 = RectangularAperture(bl, name='s2',
                             center=[0, CSS_TO_S2, 'auto'],
                             kind=('left', 'right', 'bottom', 'top'),
                             opening=[-s2hg/2, s2hg/2, -s2vg/2, s2vg/2])

    s2_screen = Screen(bl, 's2 screen', center=[0, CSS_TO_S2, 'auto'])
    s2_screen.size = 1  # mm

    eh1w2_screen = Screen(bl, 'eh1w2 screen',
                          center=[0, CSS_TO_EH1W2, 'auto'])
    eh1w2_screen.size = 1  # mm

    src2_screen = Screen(bl, 'src2 screen',
                         center=[0, CSS_TO_SRC2, 'auto'])
    src2_screen.size = 0.5  # mm

    eh2w1_screen = Screen(bl, 'eh2w1 screen',
                          center=[0, CSS_TO_EH2W1, 'auto'])
    eh2w1_screen.size = 0.5  # mm

    xshut = PolygonalAperture(bl, name='xshut',
                              center=[0, CSS_TO_XSHUT, 'auto'],
                              opening=[(4 + xshuty, 0.15), (xshuty, 1),
                                       (xshuty, -1), (4 + xshuty, -0.15)])
    xshut_screen = Screen(bl, 'xshut screen',
                          center=[0, CSS_TO_XSHUT, 'auto'])
    xshut_screen.size = 0.3  # mm

    hsc = RectangularAperture(bl, name='hsc',
                              center=[0, CSS_TO_HSC, 'auto'],
                              kind=('left', 'right', 'bottom', 'top'),
                              opening=[-1, 1, -0.130/2, +0.130/2])

    hsc_screen = Screen(bl, 'hsc screen', center=[0, CSS_TO_HSC, 'auto'])
    hsc_screen.size = 0.3  # mm

    ss1 = RectangularAperture(bl, name='ss1',
                              center=[0, CSS_TO_SS1, 'auto'],
                              kind=('left', 'right', 'bottom', 'top'),
                              opening=[-ss1hg/2, ss1hg/2, -ss1vg/2, ss1vg/2])

    ss1_screen = Screen(bl, 'ss1 screen', center=[0, CSS_TO_SS1, 'auto'])
    ss1_screen.size = 0.3

    mlRu = Multilayer(tLayer=Material('Ru', rho=9.20),
                      tThickness=20.2,
                      bLayer=Material(['B', 'C'], [4, 1], rho=2.55),
                      bThickness=20.2, nPairs=60,
                      substrate=Material('Si', rho=2.33))

    mlW = Multilayer(tLayer=Material('W', rho=14.85),
                     tThickness=26*0.346,
                     bLayer=Material(['B', 'C'], [4, 1], rho=2.65),
                     bThickness=26*0.654, nPairs=100,
                     substrate=Material('Si', rho=2.33))

    crl = DoubleParaboloidLens(bl, name='crl',
                               center=[0, CSS_TO_CRL, 'auto'],
                               pitch=np.pi/2+2*m1.pitch,
                               material=Material('C', rho=3.52, kind='lens'),
                               t=0.025,  # lens thickness
                               focus=0.025,  # lens focus f=R/2
                               nCRL=17,  # number of lenses
                               zmax=0.5875,  # lens piling thickness = 2*zmax+t
                               )
    
    crl_screen = Screen(bl, 'crl screen', center=[0, CSS_TO_CRL, 'auto'])
    crl_screen.size = 0.3  # mm

    ss2 = RectangularAperture(bl, name='ss2',
                              center=[0, CSS_TO_SS2, 'auto'],
                              kind=('left', 'right', 'bottom', 'top'),
                              opening=[-ss2hg/2, ss2hg/2, -ss2vg/2, ss2vg/2])

    ss2_screen = Screen(bl, 'ss2 screen', center=[0, CSS_TO_SS2, 'auto'])
    ss2_screen.size = 0.3

    w3 = Plate(bl, name='w3', center=[0, CSS_TO_W3, Z0], pitch=np.pi/2,
               material=Material('Be', rho=1.848, kind='plate'), t=0.2)
    
    w3_screen = Screen(bl, 'w3 screen', center=[0, CSS_TO_W3, Z0])
    w3_screen.size = 0.3  # mm

    xeye = Plate(bl, name='xeye', center=[0, CSS_TO_XEYE, 'auto'],
                material=Material('C', rho=3.52, kind='plate'), t=0.3,
                pitch=np.pi/2)
    
    xeye_screen = Screen(bl, 'xeye screen',
                         center=[0, CSS_TO_XEYE, 'auto'])
    xeye_screen.size = 0.3  # mm

    pre_screen = Screen(bl, 'pre screen',
                        center=[0, CSS_TO_SAMPLE-150, 'auto'])
    pre_screen.size = 0.3  # mm

    sample_screen = Screen(bl, 'sample screen',
                           center=[0, CSS_TO_SAMPLE, 'auto'])
    sample_screen.size = 0.3  # mm

    post_screen = Screen(bl, 'post screen',
                         center=[0, CSS_TO_SAMPLE+150, 'auto'])
    post_screen.size = 0.3  # mm
    
    bl.u23 = u23
    bl.u23_screen = u23_screen
    bl.u17 = u17
    bl.u17_screen = u17_screen
    bl.m2h = m2h
    bl.m2h_screen = m2h_screen
    bl.m2v = m2v
    bl.m2v_screen = m2v_screen
    bl.w0 = w0
    bl.w0_screen = w0_screen
    bl.ps = ps
    bl.ps_screen = ps_screen
    bl.w1 = w1
    bl.w1_screen = w1_screen
    bl.w2 = w2
    bl.w2_screen = w2_screen
    bl.mono1 = mono1
    bl.s1 = s1
    bl.s1_screen = s1_screen
    bl.ohw2_screen = ohw2_screen
    bl.eh1w1_screen = eh1w1_screen
    bl.m1 = m1
    bl.m1_screen = m1_screen
    bl.s2 = s2
    bl.s2_screen = s2_screen
    bl.bv1 = bv1
    bl.bv1_screen = bv1_screen
    bl.bv2 = bv2
    bl.bv2_screen = bv2_screen
    bl.eh1w2_screen = eh1w2_screen
    bl.src2_screen = src2_screen
    bl.eh2w1_screen = eh2w1_screen
    bl.xshut = xshut
    bl.xshut_screen = xshut_screen
    bl.hsc = hsc
    bl.hsc_screen = hsc_screen
    bl.ss1 = ss1
    bl.ss1_screen = ss1_screen
    bl.mlRu = mlRu
    bl.mlW = mlW
    bl.crl = crl
    bl.crl_screen = crl_screen
    bl.ss2 = ss2
    bl.ss2_screen = ss2_screen
    bl.w3 = w3
    bl.w0_screen = w0_screen
    bl.xeye = xeye
    bl.xeye_screen = xeye_screen
    bl.pre_screen = pre_screen
    bl.sample_screen = sample_screen
    bl.post_screen = post_screen

    return bl


def run_process(bl):

    global REPETITION

    t0 = time.time()
    folder = f"./beam_{bl.energy*1e-3:.2f}keV_std/"
    os.makedirs(folder, exist_ok=True)
    fname = folder + f"beam{REPETITION:02d}.npy"
    if os.path.isfile(fname) and not bl.force_beam:
        print(f"Reading beam from file: {fname}")
        beam = np.load(fname, allow_pickle=True).item()
    else:
        beam = bl.source.shine()
        dt = time.time() - t0
        print(f"Time needed to build rays {dt:.3f} sec")
        np.save(fname, beam)

    u17_screen = bl.u17_screen.expose(beam)
    _ = bl.m2v.propagate(beam)
    _ = bl.m2h.propagate(beam)
    beam, _, _ = bl.w0.double_refract(beam)
    w0_screen = bl.w0_screen.expose(beam)
    ps_screen = bl.ps_screen.expose(beam)
    _ = bl.ps.propagate(beam)
    s1_screen = bl.s1_screen.expose(beam)
    _ = bl.s1.propagate(beam)
    m1_screen = bl.m1_screen.expose(beam)
    beam, _ = bl.m1.reflect(beam)
    s2_screen = bl.s2_screen.expose(beam)
    _ = bl.s2.propagate(beam)
    bv1_screen = bl.bv1_screen.expose(beam)
    eh1w2_screen = bl.eh1w2_screen.expose(beam)
    eh2w1_screen = bl.eh2w1_screen.expose(beam)
    xshut_screen = bl.xshut_screen.expose(beam)
    hsc_screen = bl.hsc_screen.expose(beam)
    _ = bl.hsc.propagate(beam)
    ss1_screen = bl.ss1_screen.expose(beam)
    _ = bl.ss1.propagate(beam)
    crl_screen = bl.crl_screen.expose(beam)
    beam, _, _ = bl.crl.multiple_refract(beam)  # 2nd source
    ss2_screen = bl.ss2_screen.expose(beam)
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
               hsc_screen=hsc_screen,
               ss1_screen=ss1_screen,
               crl_screen=crl_screen,
               ss2_screen=ss2_screen,
               xeye_screen=xeye_screen,
               sample_screen=sample_screen)

    return out


xrt.backends.raycing.run.run_process = run_process
