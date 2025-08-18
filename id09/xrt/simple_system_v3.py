import os
import time
import numpy as np
from scipy.signal import savgol_filter
import xrt

from xrt.backends.raycing import BeamLine
from xrt.backends.raycing.sources import Undulator
from xrt.backends.raycing.screens import Screen
from xrt.plotter import XYCAxis, XYCPlot
from xrt.runner import run_ray_tracing

#
# !!!! ADD HERE THE FUNCTION oasys_xrt_components_code() !!!!!!
#
def oasys_xrt_components_code():

    oasys_list_of_elements = []

    oasys_list_of_elements.append('''
from xrt.backends.raycing import BeamLine
from xrt.backends.raycing.sources import Undulator

xrt_component = Undulator(
    BeamLine(),
    name="u17",
    center=[0,1250,0],
    period=17,
    n=117,
    eE=6.0,
    eI=0.2,
    eEpsilonX=0.130,
    eEpsilonZ=0.010,
    eEspread=9.4e-4,
    eSigmaX=30.0,
    eSigmaZ=5.2,
    distE="eV",
    targetE=[18070.0, 1],
    eMin=17000.0,
    eMax=18500.0,
    nrays=20000,
    )
''')

    oasys_list_of_elements.append('''

from xrt.backends.raycing import BeamLine
from xrt.backends.raycing.screens import Screen

xrt_component = Screen(
    BeamLine(),
    name="sample_screen",
    center=[0, 56289, 0],
    )
''')

    return oasys_list_of_elements


def oasys_xrt_components_objects(oasys_list_of_elements):
    # run components to get xrt objects
    oasys_list_of_elements_objects = []
    for i in range(len(oasys_list_of_elements)):
        code = oasys_list_of_elements[i]
        namespace = dict()
        exec(code, namespace)
        element = namespace["xrt_component"]
        oasys_list_of_elements_objects.append(element)
    return  oasys_list_of_elements_objects

def build_beamline(name=""):
    list_of_components_code = oasys_xrt_components_code()
    list_of_components_objects = oasys_xrt_components_objects(list_of_components_code)

    # add to bl
    bl = BeamLine()
    bl.name = name

    for element in list_of_components_objects:
        setattr(bl, element.name, element)

    bl.list_of_elements = list_of_components_code
    bl.list_of_elements_objects  = list_of_components_objects

    return bl

def run_process(bl):
    global REPETITION, dump_beams_folder, dump_beams_flag
    print("REPETITION = %d" % REPETITION)
    t0 = time.time()

    beam_out = dict()
    beam_out_list = []

    for i, element in enumerate(bl.list_of_elements_objects):
        if isinstance(element, Undulator):
            out_i = element.shine()
            beam_out[element.name] = out_i
            beam_out_list.append(out_i)
        elif isinstance(element, Screen):
            out_i = bl.sample_screen.expose(beam_out_list[i-1])
            beam_out[element.name] = out_i
            beam_out_list.append(out_i)

    if dump_beams_flag:
        for element in bl.list_of_elements_objects:
            fname = "%s%s_%02d.npy" % (dump_beams_folder, element.name, REPETITION)
            np.save(fname, beam_out_list[i])

    dt = time.time() - t0
    print("Time needed to create source and trace system %.3f sec" % dt)

    REPETITION += 1

    return beam_out


def make_plot(bl, screen, size=100, bins=1024, cbins=256):
    global dump_scores_folder

    xaxis = XYCAxis(label='x', limits=[-size/2, size/2], unit='um', bins=bins,
                    ppb=int(1024/bins))
    yaxis = XYCAxis(label='z', limits=[-size/2, size/2], unit='um', bins=bins,
                    ppb=int(1024/bins))
    caxis = XYCAxis(label='energy', limits=[bl.u17.eMin, bl.u17.eMax],
                    unit='eV', bins=cbins, fwhmFormatStr="%.2f",
                    ppb=int(512/cbins))

    fname = dump_scores_folder + "/screen_%s.png" % screen

    plot = XYCPlot(beam=screen, xaxis=xaxis, yaxis=yaxis, caxis=caxis, saveName=fname)

    return plot

def do_after_script(plots):
    global dump_scores_folder

    out_str1 = ""
    fwhm_h = np.array([plot.dx for plot in plots])
    fwhm_v = np.array([plot.dy for plot in plots])
    ecen = np.array([plot.cE/1e3 for plot in plots])
    flux = np.array([plot.flux for plot in plots])
    profile_h = np.array([[p.xaxis.binCenters, p.xaxis.total1D]
                          for p in plots])
    profile_v = np.array([[p.yaxis.binCenters, p.yaxis.total1D]
                          for p in plots])
    profile_e = np.array([[p.caxis.binCenters, p.caxis.total1D]
                          for p in plots])

    epeak = []
    for k, plot in enumerate(plots):
        fname = dump_scores_folder + "/res_%s_profile_h.txt" % plot.beam
        np.savetxt(fname, profile_h[k].T)
        fname = dump_scores_folder + "/res_%s_profile_v.txt" % plot.beam
        np.savetxt(fname, profile_v[k].T)
        fname = dump_scores_folder + "/res_%s_profile_e.txt" % plot.beam
        np.savetxt(fname, profile_e[k].T)
        e, y = np.loadtxt(fname, unpack=True)
        epeak.append(e[np.argmax(savgol_filter(y, 51, 3))])
    epeak = np.array(epeak)/1e3

    out_str2 = [
        "epeak  =  " + repr(epeak.tolist()),
        "ecen   =  " + repr(ecen.tolist()),
        "fwhm_h =  " + repr(fwhm_h.tolist()),
        "fwhm_v =  " + repr(fwhm_v.tolist()),
        "flux   =  " + repr(flux.tolist()),
    ]

    out_str2 = chr(10).join(out_str2)
    print(out_str2)
    out_str1 += chr(10) + out_str2 + chr(10)
    with open(dump_scores_folder + "/results.txt", 'w') as f:
        f.write(out_str1)

#
# main
#
global REPETITION
REPETITION = 0

#
# xrt beamline
#
bl = build_beamline(name="id09")

#
# prepare output folders (in the directory with bl.name)
#
global dump_scores_folder
dump_scores_flag = 1
dump_scores_folder = "%s/scores/" % bl.name
os.makedirs(dump_scores_folder, exist_ok=True)

global dump_beams_flag, dump_beams_folder
dump_beams_flag = 1
dump_beams_folder = "%s/beams/" % bl.name
if dump_beams_flag: os.makedirs(dump_beams_folder, exist_ok=True)


#
# declare run_process() that makes the tracing
#
xrt.backends.raycing.run.run_process = run_process

#
# define plots
#
screens_to_plot = ['sample_screen']
sizes_in_um     = [10000]

plots = []
for i in range(len(screens_to_plot)):
    plots.append(make_plot(bl, screens_to_plot[i], sizes_in_um[i]))

#
# run
#
run_ray_tracing(plots=plots,
                repeats=2,
                pickleEvery=1,
                updateEvery=1,
                beamLine=bl,
                afterScript=do_after_script,
                afterScriptArgs=[plots])


