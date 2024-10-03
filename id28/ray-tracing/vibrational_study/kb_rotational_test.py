import Shadow
import numpy
from srxraylib.util.h5_simple_writer import H5SimpleWriter
from srxraylib.plot.gol import plot as manolo_plot


def run_beamline_shadow(beam, mono_pitch=0.0, mono_roll=0.0, vfm_pitch=0.0, hfm_pitch=0.0):
    
    iwrite = 0
    
    oe0 = Shadow.OE()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()
    
    #
    # Define variables. See meaning of variables in: 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #
    #mirror used as BS-mono
    oe0.DUMMY = 100.0
    oe0.FWRITE = 1
    oe0.F_MOVE = 1
    oe0.T_IMAGE = 0.0
    oe0.T_INCIDENCE = 0.2742363247
    oe0.T_REFLECTION = 0.2742363247
    oe0.T_SOURCE = 24.0
    oe0.X_ROT = mono_pitch
    oe0.Y_ROT = mono_roll
    
    oe1.DUMMY = 100.0
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.N_SCREEN = 1
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 26.0
    
    oe2.DUMMY = 100.0
    oe2.FCYL = 1
    oe2.FHIT_C = 1
    oe2.FILE_REFL = b'/home/esrf/reyesher/OASYS/ID28/ray-tracing/Pt_mirror'
    oe2.FILE_RIP = b'/home/esrf/reyesher/OASYS/ID28/ray-tracing/M_V_present_shadow.dat'
    oe2.FMIRR = 2
    oe2.FWRITE = 1
    oe2.F_DEFAULT = 0
    oe2.F_G_S = 2
    oe2.F_MOVE = 1
    oe2.F_REFLEC = 1
    oe2.F_RIPPLE = 1
    oe2.RLEN1 = 0.12
    oe2.RLEN2 = 0.12
    oe2.RWIDX1 = 0.02
    oe2.RWIDX2 = 0.02
    oe2.SIMAG = 2.6
    oe2.SSOUR = 99.5
    oe2.THETA = 89.8453013953
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 89.8453013953
    oe2.T_REFLECTION = 89.8453013953
    oe2.T_SOURCE = 0.5
    oe2.X_ROT = vfm_pitch
    
    oe3.ALPHA = 90.0
    oe3.DUMMY = 100.0
    oe3.FCYL = 1
    oe3.FHIT_C = 1
    oe3.FILE_RIP = b'/home/esrf/reyesher/OASYS/ID28/ray-tracing/M_H_present_shadow.dat'
    oe3.FMIRR = 2
    oe3.FWRITE = 1
    oe3.F_DEFAULT = 0
    oe3.F_G_S = 2
    oe3.F_MOVE = 1
    oe3.F_RIPPLE = 1
    oe3.RLEN1 = 0.12
    oe3.RLEN2 = 0.12
    oe3.RWIDX1 = 0.02
    oe3.RWIDX2 = 0.02
    oe3.SIMAG = 1.511
    oe3.SSOUR = 100.5
    oe3.THETA = 89.33995262
    oe3.T_IMAGE = 0.0
    oe3.T_INCIDENCE = 89.33995262
    oe3.T_REFLECTION = 89.33995262
    oe3.T_SOURCE = 1.0
    oe3.X_ROT = hfm_pitch
    
    oe4.ALPHA = 90.0
    oe4.DUMMY = 100.0
    oe4.FWRITE = 3
    oe4.F_REFRAC = 2
    oe4.T_IMAGE = 0.0
    oe4.T_INCIDENCE = 0.0
    oe4.T_REFLECTION = 180.0
    oe4.T_SOURCE = 1.5    
    
    #Run SHADOW to create the source
    
    if iwrite:
        oe0.write("start.00")
    
    beam.traceOE(oe0,0)
    
    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")    
    
    #
    #run optical element 1
    #
    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")
    
    beam.traceOE(oe1,1)
    
    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")    
    
    #
    #run optical element 2
    #
    print("    Running optical element: %d"%(2))
    if iwrite:
        oe2.write("start.02")
    
    beam.traceOE(oe2,2)
    
    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")    
    
    #
    #run optical element 3
    #
    print("    Running optical element: %d"%(3))
    if iwrite:
        oe3.write("start.03")
    
    beam.traceOE(oe3,3)
    
    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")    
    
    #
    #run optical element 4
    #
    print("    Running optical element: %d"%(4))
    if iwrite:
        oe4.write("start.04")
    
    beam.traceOE(oe4,4)
    
    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")
        
    tkt = beam.histo2(1, 3, ref=23, xrange=[-5e-6, 5e-6], yrange=[-5e-6, 5e-6], nbins=201, nolost=1)

    x = tkt['bin_h_center']
    y = tkt['bin_v_center']
    histogram = tkt['histogram'] 
        
    return x, y, histogram

def run_loop(nruns, op_elem = 'vfm', sigma_vfm_pitch_urad = 0.0, sigma_hfm_pitch_urad = 0.0,
             sigma_mono_pitch_urad = 0.0, sigma_mono_roll_urad = 0.0, save_file=True):

    beam_source = Shadow.Beam()
    beam_source.load('/home/esrf/reyesher/OASYS/ID28/ray-tracing/vibrational_study/id28_3xU32_21_7keV_SS_49m.dat')
    
    sigma_vfm_pitch_deg = numpy.rad2deg(sigma_vfm_pitch_urad * 1e-6)
    sigma_hfm_pitch_deg = numpy.rad2deg(sigma_hfm_pitch_urad * 1e-6)
    sigma_mono_pitch_deg = numpy.rad2deg(sigma_mono_pitch_urad * 1e-6)
    sigma_mono_roll_deg = numpy.rad2deg(sigma_mono_roll_urad * 1e-6)

    vfm_pitch = numpy.random.normal(loc=0.0, scale=sigma_vfm_pitch_deg, size=nruns)
    hfm_pitch = numpy.random.normal(loc=0.0, scale=sigma_hfm_pitch_deg, size=nruns)
    mono_pitch = numpy.random.normal(loc=0.0, scale=sigma_mono_pitch_deg, size=nruns)
    mono_roll = numpy.random.normal(loc=0.0, scale=sigma_mono_roll_deg, size=nruns)
    
    for i in range(nruns):
        if op_elem == 'vfm':
            x, y, histogram = run_beamline_shadow(beam_source.duplicate(), vfm_pitch=vfm_pitch[i])
        elif op_elem == 'hfm':
            x, y, histogram = run_beamline_shadow(beam_source.duplicate(), hfm_pitch=hfm_pitch[i])
        elif op_elem == 'mono_pitch':
            x, y, histogram = run_beamline_shadow(beam_source.duplicate(), mono_pitch=mono_pitch[i])
        elif op_elem == 'mono_roll':
            x, y, histogram = run_beamline_shadow(beam_source.duplicate(), mono_roll=mono_roll[i])
        elif op_elem == 'mono_both':
            x, y, histogram = run_beamline_shadow(beam_source.duplicate(), mono_pitch=mono_pitch[i], mono_roll=mono_roll[i])
        elif op_elem == 'all':
            x, y, histogram = run_beamline_shadow(beam_source.duplicate(),
                                                 vfm_pitch=vfm_pitch[i],
                                                 hfm_pitch=hfm_pitch[i],
                                                 mono_pitch=mono_pitch[i],
                                                 mono_roll=mono_roll[i])
        else:
            raise RuntimeError("ERROR: Unexpected optical element")
        
        if i == 0:
            result_profile = numpy.zeros_like(histogram)
            result_profile += histogram
        else:
            result_profile += histogram

    if save_file:

        h5w = H5SimpleWriter.initialize_file(f"new_profile_opt_elem_vib_{op_elem}.h5", creator="h5_basic_writer.py")

        h5w.add_image(result_profile, x * 1e6, y * 1e6,
                      image_name="s4_weight_profile",
                      title_x="h [um]",title_y="v [um]")
        print(f"new_profile_opt_elem_vib_{op_elem}.h5 has been saved to disk")

def check_distribution(nruns, sigma_rot_urad = 0.046):

    sigma_rot_deg = numpy.rad2deg(sigma_rot_urad * 1e-6)

    rotations_deg = numpy.random.normal(loc=0.0, scale=1.0*sigma_rot_deg, size=nruns)

    t = numpy.linspace(-10, 10, nruns)

    max_rot_deg = max(rotations_deg)
    min_rot_deg = min(rotations_deg)

    max_rot_rad = numpy.deg2rad(max_rot_deg)
    min_rot_rad = numpy.deg2rad(min_rot_deg)

    max_3_sigma = 3*sigma_rot_urad * numpy.ones_like(rotations_deg)
    min_3_sigma = -3*sigma_rot_urad * numpy.ones_like(rotations_deg)

    print(f'Max rot angle {max_rot_rad} rad')
    print(f'Max rot angle {max_rot_deg} deg')

    print(f'Max rot over sigma = {max_rot_rad * 1e6 / sigma_rot_urad}')

    print(f'Min rot over sigma = {min_rot_rad * 1e6 / sigma_rot_urad}')


    manolo_plot(t, numpy.deg2rad(rotations_deg*1e6), t, max_3_sigma, t, min_3_sigma)  



if __name__  == "__main__":
    pass
    #runs = 200
    #op_elem_runs = ['hfm', 'mono_pitch', 'mono_roll', 'mono_both', 'all']
#
    #for op_elem in op_elem_runs:
    #    
    #    run_loop(runs, op_elem = op_elem, sigma_vfm_pitch_urad = 0.054,
    #             sigma_hfm_pitch_urad = 0.195, sigma_mono_pitch_urad = 0.046,
    #             sigma_mono_roll_urad = 0.120, save_file=True)
    
    