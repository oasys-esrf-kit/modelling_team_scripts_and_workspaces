
def run_beamline(sigdiz=0.000007):
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    #
    #
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
    light_source = SourceGeometrical(name='Geometrical Source Gaussian', nrays=50000, seed=5676561)
    light_source.set_spatial_type_gaussian(sigma_h=0.000000 ,sigma_v=0.000000)
    light_source.set_depth_distribution_off()
    light_source.set_angular_distribution_gaussian(sigdix=0.000000 ,sigdiz=sigdiz)
    light_source.set_energy_distribution_uniform(value_min=16990.000000, value_max=17010.000000, unit='eV')
    light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
    beam = light_source.get_beam()

    beamline.set_light_source(light_source)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Generic Crystal',
                                     boundary_shape=boundary_shape, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=17000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0, # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0.01, angle_radial=1.45421465, angle_azimuthal=0, angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_x=1.74533e-05, rotation_y=0, rotation_z=0)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element ,coordinates=coordinates, movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Generic Crystal',
                                     boundary_shape=boundary_shape, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=17000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0, # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0.01, angle_radial=1.45421465, angle_azimuthal=3.141592654, angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_x=-1.74533e-05, rotation_y=0, rotation_z=0)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element ,coordinates=coordinates, movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Generic Crystal',
                                     boundary_shape=boundary_shape, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=17000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0, # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0.01, angle_radial=1.45421465, angle_azimuthal=0, angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_x=0, rotation_y=0, rotation_z=0)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element ,coordinates=coordinates, movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX
    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
    optical_element = S4PlaneCrystal(name='Generic Crystal',
                                     boundary_shape=boundary_shape, material='Si',
                                     miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                     f_bragg_a=False, asymmetry_angle=0.0,
                                     is_thick=1, thickness=0.001,
                                     f_central=1, f_phot_cent=0, phot_cent=17000.0,
                                     file_refl='bragg.dat',
                                     f_ext=0,
                                     material_constants_library_flag=0, # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                     method_efields_management=0, # 0=new in S4; 1=like in S3
                                     )
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0.01, angle_radial=1.45421465, angle_azimuthal=3.141592654, angle_radial_out=1.45421465)
    from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
    movements = S4BeamlineElementMovements(f_move=1, offset_x=0, offset_y=0, offset_z=0, rotation_x=0, rotation_y=0, rotation_z=0)
    from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
    beamline_element = S4PlaneCrystalElement(optical_element=optical_element ,coordinates=coordinates, movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # optical element number XX

    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')

    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=4.71238898, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)
    return beam



#
# main
#
import numpy
from srxraylib.plot.gol import plot, plot_image, plot_show




y_max = 165.697
nolost = 1
npoints = 13

divergences_sigma = numpy.linspace(0, 10e-6, npoints)
fwhm = numpy.zeros(npoints)
center = numpy.zeros(npoints)
col = 1
ref = 23
col_title = "X (col 1)"

CEN_x = numpy.zeros_like(divergences_sigma)
CEN_z = numpy.zeros_like(divergences_sigma)
CEN_xp = numpy.zeros_like(divergences_sigma)
CEN_zp = numpy.zeros_like(divergences_sigma)
SD_x = numpy.zeros_like(divergences_sigma)
SD_z = numpy.zeros_like(divergences_sigma)
SD_xp = numpy.zeros_like(divergences_sigma)
SD_zp = numpy.zeros_like(divergences_sigma)


for i in range(npoints):
    beam1 = run_beamline(sigdiz=divergences_sigma[i])
    beam1.retrace(y_max, resetY=True)

    x = beam1.get_column(1, nolost=1)
    z = beam1.get_column(3, nolost=1)
    xp = beam1.get_column(4, nolost=1)
    zp = beam1.get_column(6, nolost=1)
    w = beam1.get_column(23, nolost=1)
    titles = ['x', 'z', 'xp', 'zp']

    for j, array in enumerate([x, z, xp, zp]):
        average = numpy.average(array, weights=w)
        variance = numpy.average((array - average) ** 2, weights=w)
        print(divergences_sigma[i], titles[j], average, "+/-", numpy.sqrt(variance))
        if j == 0:
            CEN_x[i] = average
            SD_x[i] = numpy.sqrt(variance)
        elif j == 1:
            CEN_z[i] = average
            SD_z[i] = numpy.sqrt(variance)
        elif j == 2:
            CEN_xp[i] = average
            SD_xp[i] = numpy.sqrt(variance)
        elif j == 3:
            CEN_zp[i] = average
            SD_zp[i] = numpy.sqrt(variance)


#
# plots
#

print("Y: ", y_max)
print("divergences (sigma) ", divergences_sigma)
print("center: ", CEN_x)

# print('Result arrays X,Y (shapes): ', out_x.shape, tkt_x['bin_center'].shape, positions.shape)
# x = tkt_x['bin_center']
# y = positions
#
# # 2D
# plot_image(out_x.T, y, 1e6 * x,
#            title="", ytitle="%s [um] (%d pixels)" % (col_title, x.size),
#            xtitle="Y [m] (%d pixels)" % (y.size), aspect="auto" )
# # FWHM
# fwhm[fwhm == 0] = "nan"
# plot(y, 1e6 * fwhm, title="FWHM",
#      xtitle="y [m]", ytitle="FHWH [um]", marker=".")
# # I0
# nx, ny = out_x.shape
# I0 = out_x.T[:, nx // 2]
# plot(y, I0, title="I at central profile", xtitle="y [m]", ytitle="I0", marker=".")
# center
plot(2.355e6 * divergences_sigma, 1e6 * CEN_x, title="CENTER",
     xtitle="divergence FWHM [urad]", ytitle="CENTROID [um]", marker="None") # , yrange=[1e6 * x_min, 1e6 * x_max])



