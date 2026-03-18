
def get_reflectivities():
    import numpy
    from crystalpy.util.calc_xcrystal import calc_xcrystal_angular_scan
    from dabax.dabax_xraylib import DabaxXraylib

    bunch_out_dict, diffraction_setup, deviations = calc_xcrystal_angular_scan(
        crystal_name='Si',
        thickness=1e-05,  # in m
        miller_h=1,
        miller_k=1,
        miller_l=1,
        asymmetry_angle=0.0,
        energy=2838.0,  # in eV
        angle_deviation_min=-0.00727221,  # in rad
        angle_deviation_max=0.00727221,  # in rad
        angle_deviation_points=1000, # 10000,
        angle_center_flag=1,
        flag_calculate_stokes=2,  # 0=No, 1=Yes, 2=Yes normalized
        jones_in=[1, 0],  # for flag_calculate_stokes > 0
        chi_deg=45,  # for flag_calculate_stokes > 0
        calculation_method=0,  # 0=Zachariasen, 1=Guigay
        is_thick=0,
        use_transfer_matrix=0,
        geometry_type_index=2,
        calculation_strategy_flag=0,  # 0=mpmath 1=numpy 2=numpy-truncated
        material_constants_library_flag=1,  # 0=xraylib, 1=dabax
        dabax=DabaxXraylib(file_f1f2="f1f2_Windt.dat"),
        # for material_constants_library_flag=1, the pointer to DabaxXraylib()
    )

    tmp = numpy.zeros((bunch_out_dict["energies"].size, 7 + 4))
    tmp[:, 0] = deviations / 4.84813681109536e-06
    tmp[:, 1] = 2838.0
    tmp[:, 2] = bunch_out_dict["phaseP"]
    tmp[:, 3] = bunch_out_dict["phaseS"]
    # tmp[:, 4] = circular polarization
    tmp[:, 5] = bunch_out_dict["intensityP"]
    tmp[:, 6] = bunch_out_dict["intensityS"]

    # from crystalpy.util.calc_xcrystal import apply_convolution_on_bunch_out_dict
    # bunch_out_dict = apply_convolution_on_bunch_out_dict(bunch_out_dict, sigma=3.5e-05)
    tmp[:, 7] = bunch_out_dict['P0']
    tmp[:, 8] = bunch_out_dict['P1']
    tmp[:, 9] = bunch_out_dict['P2']
    tmp[:, 10] = bunch_out_dict['P3']

    # from srxraylib.plot.gol import plot
    # plot(tmp[:, 0], tmp[:, 6], tmp[:, 0], tmp[:, 5], xtitle="angle", legend=["S-pol", "P-pol"])

    r_s = np.sqrt(bunch_out_dict["intensityS"]) * np.exp(1j * bunch_out_dict["phaseS"])
    r_p = np.sqrt(bunch_out_dict["intensityP"]) * np.exp(1j * bunch_out_dict["phaseP"])
    return tmp, r_s, r_p

def get_experiment(fname="2838_10microns_Si111_111SB_Stokes.dat"):
    #
    # experimental
    #
    a = np.loadtxt(fname, skiprows=5)
    angle_arcsec = a[:, 0]
    I0 = a[:, 1]
    P1 = a[:, 3]
    P2 = a[:, 5]
    P3 = np.sqrt(1 - P1 ** 2 - P2 ** 2)
    return angle_arcsec, I0, P1, P2, P3

if __name__ == "__main__":
    import numpy
    import numpy as np
    from dabax.dabax_xraylib import DabaxXraylib
    from srxraylib.plot.gol import plot, plot_show


    angle_arcsec, I0, P1, P2, P3 = get_experiment(fname="2838_10microns_Si111_111SB_Stokes.dat")
    data, r_s, r_p = get_reflectivities()


    #
    # from paper
    #
    from check_paper import _stokes_Sout_sp

    chi = np.pi/4
    S_in = np.array([1., 1., 0., 0.])

    Pc = np.zeros(r_s.size, dtype=float)
    S_out_1 = np.zeros_like(Pc)
    S_out_2 = np.zeros_like(Pc)
    S_out_3 = np.zeros_like(Pc)

    # get equation
    import sympy as sp

    S0_s = sp.Symbol("S0", real=True)
    S1_s = sp.Symbol("S1", real=True)
    S2_s = sp.Symbol("S2", real=True)
    S3_s = sp.Symbol("S3", real=True)

    S_in_s = sp.Matrix([S0_s,
                        S1_s,
                        S2_s,
                        S3_s,
                        ])
    r_s_s = sp.Symbol("r_s", complex=True)
    r_p_s = sp.Symbol("r_p", complex=True)
    chi_s = sp.Symbol("chi", real=True)

    S_out_s = _stokes_Sout_sp(S_in_s, r_s_s, r_p_s, chi_s)

    for i in range(Pc.size):

        S_out = S_out_s.subs({S0_s: S_in[0], S1_s: S_in[1], S2_s: S_in[2], S3_s: S_in[3],
                              r_s_s: r_s[i], r_p_s: r_p[i], chi_s: chi }).evalf(15)

        S_out_1[i] = sp.re(S_out[1]) / sp.re(S_out[0])
        S_out_2[i] = sp.re(S_out[2]) / sp.re(S_out[0])
        S_out_3[i] = sp.re(S_out[3]) / sp.re(S_out[0])

    plot(
         angle_arcsec, P1,
         angle_arcsec, P2,
         angle_arcsec, P3,
         -data[:, 0], S_out_1,
         -data[:, 0], S_out_2,
         -data[:, 0], S_out_3,
         xtitle="angle [arcsec]",
         legend=["P1 exp", "P2 exp", "P3 from exp", "P1-calc", "P2-calc", "Pc calc",], grid=1, yrange=[-1, 1], xrange=[-300, 300],
         title="COMPARISON (Calculation Angle Reversed; chi=%d deg)" % (np.degrees(chi)), linestyle=[None, None, None, ':', ':', ':'], color=['k','b','r','k','b','r'], show=0)

    plot(
         angle_arcsec, P1,
         angle_arcsec, P2,
         angle_arcsec, P3,
         xtitle="angle [arcsec]",
         legend=["P1 exp", "P2 exp", "P3 calculated from P1,P2"], grid=1, yrange=[-1, 1], xrange=[-300, 300],
         title="Experimental", linestyle=[None, None, None], color=['k','b','r'], show=0)


    from check_paper import _pc_giles, _pc_np
    G = _pc_giles(r_s, r_p, chi)
    plot(-data[:, 0], G,
         # -data[:, 0], S_out_3,
         -data[:, 0], _pc_np(r_s, r_p, np.pi/4, [1,1,0,0]),
         legend=['Giles','Now'], grid=1, yrange=[-1, 1], xrange=[-300, 300],
         title="Giles", linestyle=[None,":"], show=0)



    plot_show()
