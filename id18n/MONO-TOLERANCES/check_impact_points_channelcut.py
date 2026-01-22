# import shadow4
# print(dir(shadow4))
# from shadow4.beam.s4_beam import S4Beam
# from shadow4.beamline.s4_beamline import S4Beamline
# print(dir(shadow4))
import numpy
from shadow4.beamline.optical_elements.compound.s4_compound import S4Compound, S4CompoundElement


def get_optical_element_instance_channel_cut(crystal_separation=0):

    from shadow4.beamline.optical_elements.mirrors.s4_conic_mirror import S4ConicMirror, S4ConicMirrorElement

    try:
        name = self.getNode().title
    except:
        name = "Channel Cut Crystal Monochromator"

    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_conic_crystal import S4ConicCrystal

    ccc1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -0.5 * crystal_separation]
    ccc2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  1.0, -0.5 * crystal_separation]


    optical_element1 = S4ConicMirror(name='Generic Mirror', boundary_shape=boundary_shape,
                                     conic_coefficients=ccc1,
                                     f_reflec=0, f_refl=6, file_refl='<none>',
                                     refraction_index=0.99999 + 0.001j,
                                     coating_material='Ni', coating_density=8.902, coating_roughness=0)
    optical_element2 = S4ConicMirror(name='Generic Mirror', boundary_shape=boundary_shape,
                                     conic_coefficients=ccc2,
                                     f_reflec=0, f_refl=6, file_refl='<none>',
                                     refraction_index=0.99999 + 0.001j,
                                     coating_material='Ni', coating_density=8.902, coating_roughness=0)
    #
    #
    #
    return S4Compound(name=name, oe_list=[optical_element1, optical_element2])

def run_beamline(
    crystal_separation=0,
    q2=0.0,
    theta_bragg=0.0,
    ):

    angle_radial = numpy.pi / 2 - theta_bragg

    #
    #
    #
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()


    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

    light_source = SourceGeometrical(name='Geometrical Source', nrays=500, seed=5676561)
    light_source.set_spatial_type_point()
    light_source.set_depth_distribution_off()
    light_source.set_angular_distribution_flat(hdiv1=0.000000, hdiv2=0.000000, vdiv1=0.000000, vdiv2=0.000000)
    light_source.set_energy_distribution_uniform(value_min=16990.0, value_max=17010.0, unit='eV')
    light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)

    beam = light_source.get_beam()
    beam0 = beam.duplicate()

    beamline.set_light_source(light_source)

    #
    # compound
    #
    optical_element = get_optical_element_instance_channel_cut(crystal_separation=crystal_separation)

    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=0, q=0,
                                     angle_radial=angle_radial,
                                     angle_azimuthal=q2,
                                     angle_radial_out=numpy.pi / 2 + theta_bragg) # numpy.radians(90 + theta_bragg_deg))


    beamline_element = S4CompoundElement(
        optical_element=optical_element,
        coordinates=coordinates,
        movements=None,
        input_beam=beam)

    beam, footprints = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    return beam0, beam, footprints

def plot_2d(beam, footprints, irange=[-0.0001, 0.0001], frange=[-0.009,0.009]):
    from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show

    # image
    ticket = beam.histo2(1, 3, nbins_h=100, nbins_v=100,
                         xrange=irange, yrange=irange, nolost=1, ref=23)

    title = "BEAM I: %.1f " % ticket['intensity']
    if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
    if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

    plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                               title=title, xtitle="column 1", ytitle="column 3",
                               cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=0)
    # footprint 1
    ticket = footprints[0].histo2(2, 1, nbins_h=100, nbins_v=100,
                                  xrange=frange, yrange=frange, nolost=1, ref=23)

    title = "FOOTPRINT FIRST CRYSTAL: %.1f " % ticket['intensity']
    if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
    if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

    plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                               title=title, xtitle="column 2", ytitle="column 1",
                               cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=0)

    # footprint 2
    ticket = footprints[1].histo2(2, 1, nbins_h=100, nbins_v=100,
                                  xrange=frange, yrange=frange, nolost=1, ref=23)

    title = "FOOTPRINT SECOND CRYSTAL: %.1f " % ticket['intensity']
    if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
    if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

    plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                               title=title, xtitle="column 2", ytitle="column 1",
                               cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=1)

def centroids(beam1, title='', factor=1.0):
    x = beam1.get_column(1, nolost=1)
    y = beam1.get_column(2, nolost=1)
    z = beam1.get_column(3, nolost=1)
    xp = beam1.get_column(4, nolost=1)
    yp = beam1.get_column(5, nolost=1)
    zp = beam1.get_column(6, nolost=1)
    w = beam1.get_column(23, nolost=1)

    arrays = [x, y, z, xp, yp, zp]
    titles = ['x', 'y', 'z', 'xp', 'yp', 'zp']

    # arrays = [x, y, z]
    # titles = ['x', 'y', 'z']

    print("\n-------------------------", title, "factor= %f" % factor)
    centroids = []
    sd = []
    for j, array in enumerate(arrays):
        average = numpy.average(array, weights=w)
        variance = numpy.average((array - average) ** 2, weights=w)
        stdev = numpy.sqrt(variance)
        print(titles[j], factor * average, "+/-", factor * stdev)
        centroids.append(average)
        sd.append(stdev)
    print("intensity: ", beam1.intensity(nolost=1))
    return centroids, sd

if __name__ in ["__main__"]:
    from srxraylib.plot.gol import plot, plot_show

    #
    # main
    #
    # theta_bragg_deg = 6.679637445440589
    crystal_separation = 0.005
    q2 = 1

    # theta_bragg = numpy.radians(theta_bragg_deg)
    # q1 = crystal_separation / numpy.sin(theta_bragg)

    # print(">>> crystal_separation, crystal path: ", crystal_separation, q1)

    # print(">>> theory beam separation 2: ", crystal_separation / numpy.sin(theta_bragg) * numpy.sin(2 * theta_bragg),
    #       2 * crystal_separation * numpy.cos(theta_bragg))



    ENERGY_keV = numpy.linspace(7.0, 40, 25)
    # OFFSET = numpy.linspace(numpy.radians(1), numpy.radians(20), 25)
    import scipy.constants as codata
    wavelength = codata.h * codata.c / codata.e / (ENERGY_keV * 1000)
    d_spacing =  3.135416e-10
    OFFSET = numpy.asin( wavelength / 2 / d_spacing)
    print(wavelength, numpy.degrees(OFFSET))

    CEN1_y = numpy.zeros_like(OFFSET)
    CEN2_y = numpy.zeros_like(OFFSET)
    CENSD1_y = numpy.zeros_like(OFFSET)
    CENSD2_y = numpy.zeros_like(OFFSET)
    IMG_z = numpy.zeros_like(OFFSET)
    IMGSD_z = numpy.zeros_like(OFFSET)

    for i, rotation in enumerate(OFFSET):
        print ("iteration %d of %d" % (i, OFFSET.size))
        beam0, beam1, footprints = run_beamline(crystal_separation=crystal_separation,
                                                theta_bragg=rotation,
                                                q2=q2)
        print("################################################################# Bragg angle [deg]: ", numpy.degrees(rotation))
        centroid, sd = centroids(beam0, title='SOURCE', factor=1e6)

        centroid, sd = centroids(footprints[0], title='CRYSTAL 1', factor=1e3)
        CEN1_y[i] = centroid[1]
        CENSD1_y[i] = sd[1]

        centroid, sd = centroids(footprints[1], title='CRYSTAL 2', factor=1e3)
        CEN2_y[i] = centroid[1]
        CENSD2_y[i] = sd[1]

        centroid, sd = centroids(beam1, title='IMAGE', factor=1e6)
        IMG_z[i] = centroid[2]
        IMGSD_z[i] = sd[2]

    # print(CEN1_y,  (CEN1_y - CEN1_y[0]))
    plot(numpy.degrees(OFFSET), 1e3 * (CEN1_y ),
         numpy.degrees(OFFSET), 1e3 * (CEN2_y ),
         numpy.degrees(OFFSET), 1e3 * (IMG_z  ),
         title="d=%.3f m, q2=%.3f m" % (crystal_separation, q2),
         xtitle="Bragg angle [deg]", ytitle="centroid [mm]", legend=["y crystal 1", "y crystal 2", "z image"],
         linestyle=[None, None, None, '--', ':',':'],
         figsize=(10,6), grid=1, show=0)

    plot(ENERGY_keV, 1e3 * (CEN1_y ),
         ENERGY_keV, 1e3 * (CEN2_y ),
         ENERGY_keV, 1e3 * (IMG_z  ),
         title="d=%.3f m, q2=%.3f m" % (crystal_separation, q2),
         xtitle="Photon energy [keV]", ytitle="centroid [mm]", legend=["y crystal 1", "y crystal 2", "z image"],
         linestyle=[None, None, None, '--', ':',':'],
         figsize=(10,6), grid=1, show=1)

