import Shadow
import numpy
import pandas as pd
import matplotlib.pyplot as plt

def run_shadow_save_beams(n_rays=100000):
    """ Function to get all for each optical element: BEAM, P and Q """
    BEAM = []
    P = []
    Q = []

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0
    
    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()
    oe5 = Shadow.OE()
    oe6 = Shadow.OE()
    oe7 = Shadow.IdealLensOE()
    oe8 = Shadow.OE()
    
    #
    # Define variables. See meaning of variables in: 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #
    
    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = n_rays
    oe0.PH1 = 15000.001
    oe0.PH2 = 15000.001
    oe0.SIGDIX = 6.225556166616453e-06
    oe0.SIGDIZ = 4.644013456052442e-06
    oe0.SIGMAX = 3.0313512755844026e-05
    oe0.SIGMAZ = 4.5916036770194325e-06
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0
    
    oe1.DUMMY = 100.0
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.N_SCREEN = 1
    oe1.RX_SLIT = numpy.array([0.0005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.RZ_SLIT = numpy.array([0.0005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 27.066
    
    oe2.DUMMY = 100.0
    oe2.FHIT_C = 1
    oe2.FILE_REFL = b'Pd_reflec.dat'
    oe2.FILE_RIP = b'mirror_shadow.dat'
    oe2.FMIRR = 3
    oe2.FWRITE = 1
    oe2.F_EXT = 1
    oe2.F_G_S = 2
    oe2.F_REFLEC = 1
    oe2.F_RIPPLE = 1
    oe2.RLEN1 = 0.3
    oe2.RLEN2 = 0.3
    oe2.RWIDX1 = 0.065
    oe2.RWIDX2 = 0.065
    oe2.R_MAJ = 4042.4
    oe2.R_MIN = 0.046
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 89.806913223
    oe2.T_REFLECTION = 89.806913223
    oe2.T_SOURCE = 17.474
    
    oe3.DUMMY = 100.0
    oe3.FWRITE = 3
    oe3.F_REFRAC = 2
    oe3.F_SCREEN = 1
    oe3.N_SCREEN = 1
    oe3.T_IMAGE = 0.0
    oe3.T_INCIDENCE = 0.0
    oe3.T_REFLECTION = 180.0
    oe3.T_SOURCE = 4.0
    
    oe4.DUMMY = 100.0
    oe4.FWRITE = 3
    oe4.F_REFRAC = 2
    oe4.F_SCREEN = 1
    oe4.N_SCREEN = 1
    oe4.T_IMAGE = 0.0
    oe4.T_INCIDENCE = 0.0
    oe4.T_REFLECTION = 180.0
    oe4.T_SOURCE = 4.05
    
    oe5.DUMMY = 100.0
    oe5.FWRITE = 3
    oe5.F_REFRAC = 2
    oe5.F_SCREEN = 1
    oe5.N_SCREEN = 1
    oe5.T_IMAGE = 0.0
    oe5.T_INCIDENCE = 0.0
    oe5.T_REFLECTION = 180.0
    oe5.T_SOURCE = 1.5
    
    oe6.DUMMY = 100.0
    oe6.FWRITE = 3
    oe6.F_REFRAC = 2
    oe6.F_SCREEN = 1
    oe6.N_SCREEN = 1
    oe6.T_IMAGE = 0.0
    oe6.T_INCIDENCE = 0.0
    oe6.T_REFLECTION = 180.0
    oe6.T_SOURCE = 0.43
    
    oe7.T_SOURCE = 0.97
    oe7.T_IMAGE = 0.0
    oe7.focal_x = 0.627
    oe7.focal_z = 0.627
    
    oe8.DUMMY = 100.0
    oe8.FWRITE = 3
    oe8.F_REFRAC = 2
    oe8.T_IMAGE = 0.0
    oe8.T_INCIDENCE = 180.0
    oe8.T_REFLECTION = 0.0
    oe8.T_SOURCE = 0.8
    
    
    
    #Run SHADOW to create the source
    
    if iwrite:
        oe0.write("start.00")
    
    beam.genSource(oe0)
    
    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")
    
    BEAM.append(beam.duplicate())
    # Please notice that for the source
    # we are not saving P or Q, only the shadow beam 
    
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

    BEAM.append(beam.duplicate())
    P.append(oe1.T_SOURCE)
    Q.append(oe1.T_IMAGE)
    
    
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

    BEAM.append(beam.duplicate())
    P.append(oe2.T_SOURCE)
    Q.append(oe2.T_IMAGE)    
    
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
    
    BEAM.append(beam.duplicate())
    P.append(oe3.T_SOURCE)
    Q.append(oe3.T_IMAGE)  
    
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
    
    BEAM.append(beam.duplicate())
    P.append(oe4.T_SOURCE)
    Q.append(oe4.T_IMAGE)  
    #
    #run optical element 5
    #
    print("    Running optical element: %d"%(5))
    if iwrite:
        oe5.write("start.05")
    
    beam.traceOE(oe5,5)
    
    if iwrite:
        oe5.write("end.05")
        beam.write("star.05")

    BEAM.append(beam.duplicate())
    P.append(oe5.T_SOURCE)
    Q.append(oe5.T_IMAGE)  
    
    #
    #run optical element 6
    #
    print("    Running optical element: %d"%(6))
    if iwrite:
        oe6.write("start.06")
    
    beam.traceOE(oe6,6)
    
    if iwrite:
        oe6.write("end.06")
        beam.write("star.06")

    BEAM.append(beam.duplicate())
    P.append(oe6.T_SOURCE)
    Q.append(oe6.T_IMAGE)  
    
    #
    #run optical element 7
    #
    print("    Running optical element: %d"%(7))
    if iwrite:
        oe7.write("start.07")
    
    beam.traceIdealLensOE(oe7,7)
    
    if iwrite:
        oe7.write("end.07")
        beam.write("star.07")
    
    BEAM.append(beam.duplicate())
    P.append(oe7.T_SOURCE)
    Q.append(oe7.T_IMAGE)  
    
    
    #
    #run optical element 8
    #
    print("    Running optical element: %d"%(8))
    if iwrite:
        oe8.write("start.08")
    
    beam.traceOE(oe8,8)
    
    if iwrite:
        oe8.write("end.08")
        beam.write("star.08")
    
    BEAM.append(beam.duplicate())
    P.append(oe8.T_SOURCE)
    Q.append(oe8.T_IMAGE)  
    
    
    Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    
    # BEAM has an extra element than P and Q, since source has neither P or Q
    return BEAM, P, Q

def cal_dist_beamsizes(beam_list, p_list, num_dist=10, save_file=True):

    """ Function to get distances and beam sizes from BEAMs, Ps
    assuming that all Qs = 0 """   
    
    df_results = pd.DataFrame(columns=['Distance(m)', 'FWHM_H(um)', 'FWHM_V(um)', 'STD_H(um)', 'STD_V(um)'])
   
    for i in range(len(beam_list) - 1):
        if i == 0:
            d_tot = 0
        else:
            d_tot += p_list[i-1]        

        for step_p in numpy.linspace(0, p_list[i], num_dist):          

            b = beam_list[i].duplicate()
            b.retrace(step_p)
            # get the FWHM
            tkt = b.histo2(1, 3, ref=23, nbins=201, nolost=1)
            fwhm_h = tkt['fwhm_h'] * 1e6 #microns
            fwhm_v = tkt['fwhm_v'] * 1e6 #microns
            #STD
            std_h = b.get_standard_deviation(1, nolost=1, ref=23) * 1e6 #microns
            std_v = b.get_standard_deviation(3, nolost=1, ref=23) * 1e6 #microns

            df_results = pd.concat([df_results, pd.DataFrame([[d_tot + step_p, fwhm_h, fwhm_v, std_h, std_v]],
                                   columns=df_results.columns)], ignore_index=True)
            
    if save_file:
        df_results.to_csv('Beamsize_along_id09.csv', index=False)
        print('Beamsize_along_id09.csv has been saved to disk')

    return df_results

def plot_results(csv_file, oe_p=None, oe_names=None, f_size=12, x_margin=0.5):
    
    if oe_p == None:
        oe_p = [27.066, 17.474, 4.0, 4.05, 1.5, 0.43, 0.97, 0.8]
    if oe_names == None:
       oe_list = ['Source', 'PS', 'TM', 'S-OH2', '2nd source', 'S-ExpH', 'Chopper', '2nd focusing', 'Sample']

    oe_dist_to_source = [0] # first the source
    oe_to_src = 0 

    for oe_p in oe_p:
        oe_to_src += oe_p 
        oe_dist_to_source.append(oe_to_src)
    
    df = pd.read_csv(csv_file)

    fig, (ax1, ax2) = plt.subplots(2)    
    fig.suptitle('Beamsize along ID09: Secondary source project')
    #plot FWHM
    ax1.plot(df['Distance(m)'], df['FWHM_H(um)'], label='Horizontal')
    ax1.plot(df['Distance(m)'], df['FWHM_V(um)'], label='Vertical')
    ax1.legend(fontsize=f_size)
    
    ax1.set_xlabel("Distance to source (m)", fontsize= f_size)      
    ax1.set_ylabel("FWHM ($\mu$m)", fontsize= f_size)

    ax1.set_xlim(xmin=df['Distance(m)'].iloc[0] - x_margin, xmax=df['Distance(m)'].iloc[-1] + x_margin)
    ax1.tick_params(axis='both', which='major', labelsize=f_size)

    ax1.grid(axis='y')

    # plot for STD
    ax2.plot(df['Distance(m)'], 2.33 * df['STD_H(um)'])
    ax2.plot(df['Distance(m)'], 2.33 * df['STD_V(um)'])

    ax2.set_xlabel("Main optical elements", fontsize= f_size)    
    ax2.set_ylabel("2.35 x STD ($\mu$m)", fontsize= f_size)

    ax2.set_xticks(oe_dist_to_source, oe_list, rotation='vertical')
    ax2.tick_params(axis='both', which='major', labelsize=f_size)

    ax2.set_xlim(xmin=df['Distance(m)'].iloc[0] - x_margin, xmax=df['Distance(m)'].iloc[-1] + x_margin)
    ax2.grid(axis='both')     
    
    plt.show()   


if __name__ == "__main__":
    #pass
    
    #b, p, q = run_shadow_save_beams(n_rays=500000)
    #df = cal_dist_beamsizes(b, p, num_dist=10, save_file=True)

    plot_results('Beamsize_along_id09.csv')    
    