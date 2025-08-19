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

# def get_oasys_list_of_elements():

def get_oasys_list_of_elements():

    oasys_list_of_elements = []

    oasys_list_of_elements.append(
        Undulator(BeamLine(),
                  name='u17', center=[0, 1250, 0], period=17,
                    n=117, eE=6, eI=0.2, eEpsilonX=0.130, eEpsilonZ=0.010,
                    eEspread=9.4e-4, eSigmaX=30, eSigmaZ=5.2, distE='eV',
                    targetE=[18.07e3, 1], eMin=17000, eMax=18500, nrays=20e3)
    )


    oasys_list_of_elements.append(
        Screen(BeamLine(), name='sample_screen', center=[0, 56289, 0])
    )

    return oasys_list_of_elements


def build_beamline(name=""):

    # bl = BeamLine()
    # bl.name = 'id09'

    #
    # define a sorted list of elements
    #
    oasys_list_of_elements = get_oasys_list_of_elements()

    #
    # add to bl
    #
    bl = BeamLine()
    bl.name = name

    for i, element in enumerate(oasys_list_of_elements):
        setattr(bl, element.name, element)

    bl.oasys_list_of_elements  = oasys_list_of_elements

    return bl


def run_process(bl):

    global REPETITION, output_folder_beams, output_beams_dump

    print(f"\nREPETITION = {REPETITION}")
    t0 = time.time()

    beam_out = {}
    beam_out_list = []

    oasys_list_of_elements = bl.oasys_list_of_elements
    for i in range(len(oasys_list_of_elements)):
        element = oasys_list_of_elements[i]
        element_name = element.name
        if isinstance(element, Undulator):
            out_i = element.shine()
            beam_out[element_name] = out_i
            beam_out_list.append(out_i)
        elif isinstance(element, Screen):
            out_i = bl.sample_screen.expose(beam_out_list[i-1])
            beam_out[element_name] = out_i
            beam_out_list.append(out_i)

    if output_beams_dump:
        for i in range(len(oasys_list_of_elements)):
            element = oasys_list_of_elements[i]
            element_name = element.name
            fname = output_folder_beams + f"{element_name}_{REPETITION:02d}.npy"
            np.save(fname, beam_out_list[i])

    dt = time.time() - t0
    print(f"Time needed to create source and trace system {dt:.3f} sec")

    REPETITION += 1

    return beam_out


def make_plot(bl, screen, size=100, bins=1024, cbins=256):

    global output_folder_scores

    xaxis = XYCAxis(label='x', limits=[-size/2, size/2], unit='um', bins=bins,
                    ppb=int(1024/bins))
    yaxis = XYCAxis(label='z', limits=[-size/2, size/2], unit='um', bins=bins,
                    ppb=int(1024/bins))
    caxis = XYCAxis(label='energy', limits=[bl.u17.eMin, bl.u17.eMax],
                    unit='eV', bins=cbins, fwhmFormatStr="%.2f",
                    ppb=int(512/cbins))

    fname = output_folder_scores + f"/screen_{screen}.png"

    plot = XYCPlot(beam=screen, xaxis=xaxis, yaxis=yaxis, caxis=caxis,
                   saveName=fname)

    return plot



def do_after_script(plots):

    global output_folder_scores

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
        fname = output_folder_scores + f"/res_{plot.beam}_profile_h.txt"
        np.savetxt(fname, profile_h[k].T)
        fname = output_folder_scores + f"/res_{plot.beam}_profile_v.txt"
        np.savetxt(fname, profile_v[k].T)
        fname = output_folder_scores + f"/res_{plot.beam}_profile_e.txt"
        np.savetxt(fname, profile_e[k].T)
        e, y = np.loadtxt(fname, unpack=True)
        epeak.append(e[np.argmax(savgol_filter(y, 51, 3))])
    epeak = np.array(epeak)/1e3

    out_str2 = [
        f"epeak = {epeak}",
        f"ecen = {ecen}",
        f"fwhm_h = {fwhm_h}",
        f"fwhm_v = {fwhm_v}",
        f"flux = {flux}",
    ]

    out_str2 = '\n'.join(out_str2)
    print(out_str2)
    out_str1 += '\n' + out_str2 + '\n'
    with open(output_folder_scores + "/results.txt", 'w') as f:
        f.write(out_str1)


def main():
    global REPETITION
    REPETITION = 0

    global output_folder_scores
    output_folder_scores = f"tmp/scores/"
    os.makedirs(output_folder_scores, exist_ok=True)

    global output_beams_dump, output_folder_beams
    output_beams_dump = False
    output_folder_beams = f"tmp/beams/"
    if output_beams_dump: os.makedirs(output_folder_beams, exist_ok=True)

    #
    # beamline elements
    #
    bl = build_beamline(name="id09")

    #
    # run_process() makes the tracing
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

if __name__ == '__main__':
    main()
