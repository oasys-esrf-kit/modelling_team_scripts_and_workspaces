# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 20:52:08 2025

@author: BRUMUND
"""
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

def run_beamline(undulator, ml1_deform = 'none', ml2_deform = 'none'):
    print(f'running beamline with: ml1:{ml1_deform}, ml2:{ml2_deform}, id={undulator["id"]}, nrays={undulator["nrays"]}')
    
    from shadow4.beamline.s4_beamline import S4Beamline
    
    beamline = S4Beamline()
    
    
    beamline.set_light_source(undulator['light_source'])
    beam = undulator['beam']
    
    # optical element number 01 - ML1
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.000496, x_right=0.000496, y_bottom=-0.0496, y_top=0.0496)
       
    from shadow4.beamline.optical_elements.multilayers.s4_plane_multilayer import S4PlaneMultilayer
    optical_element = S4PlaneMultilayer(name='Plane Multilayer',boundary_shape=boundary_shape,
        f_refl=2,file_refl='C:/Users/brumund/Work Folders/Documents/Oasys/working_directory/Pd1.36_B4C0.8_C0.8-11keV-th1.12deg-shifted_60eV.txt', structure='[C,Pt]x30+Si', period=50.000000, Gamma=0.400000)
    
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=31.5, q=0.25, angle_radial=1.551248639, angle_azimuthal=1.570796327, angle_radial_out=1.551248639)
    
    if ml1_deform != 'none':  
        print('considering ml1 deformation: {ml1_deform}')
        
        ideal_multilayer = optical_element
        boundary_shape = None
                
        from shadow4.beamline.optical_elements.multilayers.s4_numerical_mesh_multilayer import S4NumericalMeshMultilayer
        optical_element = S4NumericalMeshMultilayer(name='Numerical Mesh Multilayer',boundary_shape=boundary_shape,
            xx=None,yy=None,zz=None,surface_data_file=f'C:/Users/brumund/Work Folders/Documents/Oasys/working_directory/{ml1_deform}',
            f_refl=0,file_refl='', structure='[B/W]x50+Si', period=25.000000, Gamma=0.500000)
        
        numerical_mesh_multilayer = optical_element
        
        from shadow4.beamline.optical_elements.multilayers.s4_additional_numerical_mesh_multilayer import S4AdditionalNumericalMeshMultilayer
        optical_element = S4AdditionalNumericalMeshMultilayer(name='ideal + error Multilayer', ideal_multilayer=ideal_multilayer, numerical_mesh_multilayer=numerical_mesh_multilayer)
            
        
        movements = None
        from shadow4.beamline.optical_elements.multilayers.s4_additional_numerical_mesh_multilayer import S4AdditionalNumericalMeshMultilayerElement
        beamline_element = S4AdditionalNumericalMeshMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)
    else:
        print('considering no ml1 deformation')
        movements = None
        from shadow4.beamline.optical_elements.multilayers.s4_plane_multilayer import S4PlaneMultilayerElement
        beamline_element = S4PlaneMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)
    beam, footprint = beamline_element.trace_beam()
    
    beamline.append_beamline_element(beamline_element)
    
    # optical element number 02 - ML2
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.000498018, x_right=0.000496046, y_bottom=-0.0282957, y_top=0.0495699)
       
    from shadow4.beamline.optical_elements.multilayers.s4_plane_multilayer import S4PlaneMultilayer
    optical_element = S4PlaneMultilayer(name='Plane Multilayer',boundary_shape=boundary_shape,
        f_refl=2,file_refl='C:/Users/brumund/Work Folders/Documents/Oasys/working_directory/Pd1.36_B4C0.8_C0.8-11keV-th1.12deg-shifted_60eV.txt', structure='[C,Pt]x30+Si', period=50.000000, Gamma=0.400000)
    
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.25, q=1, angle_radial=1.551248639, angle_azimuthal=3.141592654, angle_radial_out=1.551248639)
    
    if ml2_deform != 'none':  
        print('considering ml2 deformation: {ml2_deform}')
        ideal_multilayer = optical_element
        boundary_shape = None
                
        from shadow4.beamline.optical_elements.multilayers.s4_numerical_mesh_multilayer import S4NumericalMeshMultilayer
        optical_element = S4NumericalMeshMultilayer(name='Numerical Mesh Multilayer',boundary_shape=boundary_shape,
            xx=None,yy=None,zz=None,surface_data_file=f'C:/Users/brumund/Work Folders/Documents/Oasys/working_directory/{ml2_deform}',
            f_refl=0,file_refl='', structure='[B/W]x50+Si', period=25.000000, Gamma=0.500000)
        
        numerical_mesh_multilayer = optical_element
        from syned.beamline.shape import Rectangle
        boundary_shape = Rectangle(x_left=-0.000498018, x_right=0.000496046, y_bottom=-0.0282957, y_top=0.0495699)
        
        from shadow4.beamline.optical_elements.multilayers.s4_additional_numerical_mesh_multilayer import S4AdditionalNumericalMeshMultilayer
        optical_element = S4AdditionalNumericalMeshMultilayer(name='ideal + error Multilayer', ideal_multilayer=ideal_multilayer, numerical_mesh_multilayer=numerical_mesh_multilayer)
        
        
        movements = None
        from shadow4.beamline.optical_elements.multilayers.s4_additional_numerical_mesh_multilayer import S4AdditionalNumericalMeshMultilayerElement
        beamline_element = S4AdditionalNumericalMeshMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)
    else:
        print('considering no ml2 deformation')
        movements = None
        from shadow4.beamline.optical_elements.multilayers.s4_plane_multilayer import S4PlaneMultilayerElement
        beamline_element = S4PlaneMultilayerElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)
    beam, footprint = beamline_element.trace_beam()
    
    beamline.append_beamline_element(beamline_element)
																												
    # optical element number 03 - Empty element, rotate back
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')
    
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.570796327, angle_azimuthal=1.570796327, angle_radial_out=1.570796327)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)
    
    beam, footprint = beamline_element.trace_beam()
    
    beamline.append_beamline_element(beamline_element)
    
    # optical element number 04 - KB Slit
    beam_slit = beam
    beam_slit.retrace(72.630000)
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.00102, x_right=0.00102, y_bottom=-0.00045, y_top=0.00045)
    
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
    optical_element = S4Screen(name='Slit (Size KBs)', boundary_shape=boundary_shape,
        i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
        i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)
    
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=72.63, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
    beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)
    
    beam, footprint = beamline_element.trace_beam()
    
    beamline.append_beamline_element(beamline_element)
    
    # optical element number 05 - KB 1
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.0025, x_right=0.0025, y_bottom=-0.15, y_top=0.15)
            
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
        surface_calculation=0,
        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
        p_focus=105.630000, q_focus=1.370000, grazing_angle=0.003000,
        is_cylinder=1, cylinder_direction=0, convexity=1,
        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999+0.001j,
        coating_material='Si', coating_density=2.33, coating_roughness=0)
    
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=1, q=0.26, angle_radial=1.567796327, angle_azimuthal=0, angle_radial_out=1.567796327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)
    
    beam, footprint = beamline_element.trace_beam()
    
    beamline.append_beamline_element(beamline_element)
    
    # optical element number 06 - KB 2
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(x_left=-0.0025, x_right=0.0025, y_bottom=-0.34, y_top=0.34)
            
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
    optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
        surface_calculation=0,
        min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
        p_focus=106.150000, q_focus=0.850000, grazing_angle=0.003000,
        is_cylinder=1, cylinder_direction=0, convexity=1,
        f_reflec=0, f_refl=0, file_refl='<none>', refraction_index=0.99999+0.001j,
        coating_material='Si', coating_density=2.33, coating_roughness=0)
    
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0.26, q=0.85, angle_radial=1.567796327, angle_azimuthal=1.570796327, angle_radial_out=1.567796327)
    movements = None
    from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
    beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)
    
    beam, footprint = beamline_element.trace_beam()
    
    beamline.append_beamline_element(beamline_element)
    
    # optical element number 07 - Empty element, rotate back
    
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
    optical_element = S4Empty(name='Empty Element')
    
    from syned.beamline.element_coordinates import ElementCoordinates
    coordinates = ElementCoordinates(p=0, q=0, angle_radial=1.570796327, angle_azimuthal=1.570796327, angle_radial_out=1.570796327)
    from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
    beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)
    
    beam, footprint = beamline_element.trace_beam()
    
    beamline.append_beamline_element(beamline_element)
    return beam, beam_slit

def run_source(undulator = 'cpmu16', nrays = 1e6, ng_e = 31):
    nrays = int(nrays)
    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam
    electron_beam = S4ElectronBeam(energy_in_GeV=6,energy_spread=0.001,current=0.2)
    electron_beam.set_sigmas_all(sigma_x=3.01836e-05, sigma_y=3.63641e-06, sigma_xp=4.36821e-06, sigma_yp=1.37498e-06)
    
    if undulator == 'cpmu16':
        # magnetic structure
        from shadow4.sources.undulator.s4_undulator import S4Undulator
        source = S4Undulator(
            K_vertical        = 1.325129, # syned Undulator parameter
            period_length     = 0.0164, # syned Undulator parameter
            number_of_periods = 121.951, # syned Undulator parameter
            emin              = 10850.0, # Photon energy scan from energy (in eV)
            emax              = 11200.0, # Photon energy scan to energy (in eV)
            ng_e              = ng_e, # Photon energy scan number of points
            # ng_e              = 3, # Photon energy scan number of points
            maxangle          = 3e-05, # Maximum radiation semiaperture in RADIANS
            ng_t              = 100, # Number of points in angle theta
            ng_p              = 21, # Number of points in angle phi
            ng_j              = 20, # Number of points in electron trajectory (per period) for internal calculation only
            code_undul_phot   = 'internal', # internal, pysru, srw
            flag_emittance    = 1, # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_size         = 0, # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
            distance          = 100.0, # distance to far field plane
            srw_range         = 0.05, # for SRW backpropagation, the range factor
            srw_resolution    = 50, # for SRW backpropagation, the resolution factor
            srw_semianalytical= 0, # for SRW backpropagation, use semianalytical treatement of phase
            magnification     = 0.05, # for internal/wofry backpropagation, the magnification factor
            flag_backprop_recalculate_source      = 0, # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
            flag_backprop_weight = 0, # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
            weight_ratio         = 0.5, # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
            flag_energy_spread   = 0, # for monochromatod sources, apply (1) or not (0) electron energy spread correction
            )
    elif undulator == 'ivu21':
    # magnetic structure
        from shadow4.sources.undulator.s4_undulator import S4Undulator
        source = S4Undulator(
            K_vertical        = 0.966, # syned Undulator parameter
            period_length     = 0.021, # syned Undulator parameter
            number_of_periods = 95.238, # syned Undulator parameter
            emin              = 10850.0, # Photon energy scan from energy (in eV)
            emax              = 11200.0, # Photon energy scan to energy (in eV)
            ng_e              = ng_e, # Photon energy scan number of points
            # ng_e              = 3, # Photon energy scan number of points
            maxangle          = 3e-05, # Maximum radiation semiaperture in RADIANS
            ng_t              = 100, # Number of points in angle theta
            ng_p              = 21, # Number of points in angle phi
            ng_j              = 20, # Number of points in electron trajectory (per period) for internal calculation only
            code_undul_phot   = 'internal', # internal, pysru, srw
            flag_emittance    = 1, # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_size         = 0, # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
            distance          = 100.0, # distance to far field plane
            srw_range         = 0.05, # for SRW backpropagation, the range factor
            srw_resolution    = 50, # for SRW backpropagation, the resolution factor
            srw_semianalytical= 0, # for SRW backpropagation, use semianalytical treatement of phase
            magnification     = 0.05, # for internal/wofry backpropagation, the magnification factor
            flag_backprop_recalculate_source      = 0, # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
            flag_backprop_weight = 0, # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
            weight_ratio         = 0.5, # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
            flag_energy_spread   = 0, # for monochromatod sources, apply (1) or not (0) electron energy spread correction
            )
    
    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,nrays=nrays,seed=5676561)
    beam = light_source.get_beam()
    
    output = {'beam' : beam, 'light_source' : light_source, 'id' : undulator, 'nrays' : nrays}
    return output

def plot_histo2(beam, nbinsh = 100, nbinsv = 100):
    ticket = beam.histo2(1, 3, nbins_h=nbinsh, nbins_v=nbinsv, nolost=1, ref=22)
    
    title = "I: %.1f " % ticket['intensity']
    if ticket['fwhm_h'] is not None: title += "FWHM H: %f um " % (ticket['fwhm_h']*1e6)
    if ticket['fwhm_v'] is not None: title += "FWHM V: %f um" % (ticket['fwhm_v']*1e6)
    
    plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
        title=title, xtitle="horizontal (mm)", ytitle="vertical (mm)",
        cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=1)
    return 

def run_dataFrame(design_df, nrays = 1e6):
    import numpy as np
    
    cpmu16 = run_source('cpmu16', nrays)
    ivu21 = run_source('ivu21', nrays)
    design_df['int_slit'] = 0
    design_df['fwhmH_slit'] = 0
    design_df['fwhmV_slit'] = 0
    design_df['centrH_slit'] = 0
    design_df['centrV_slit'] = 0
    design_df['int'] = 0
    design_df['fwhmH'] = 0
    design_df['fwhmV'] = 0
    design_df['centrH'] = 0
    design_df['centrV'] = 0
    beams_slit = []
    beams = []
    for i, row in design_df.iterrows():
        if row['id'] == 'cpmu16': current_id = cpmu16
        if row['id'] == 'ivu21': current_id = ivu21
        beam, beam_slit = run_beamline(current_id,
                                       ml1_deform=row['ml1_deform'], 
                                       ml2_deform=row['ml2_deform'])
        ticket_slit = beam_slit.histo2(1, 3, nbins=500, nolost=1, ref=22)
        ticket = beam.histo2(1, 3, nbins=500, nolost=1, ref=22)
        design_df['int_slit'][i] = ticket_slit['intensity']
        design_df['fwhmH_slit'][i] = ticket_slit['fwhm_h']
        design_df['fwhmV_slit'][i] = ticket_slit['fwhm_v']
        design_df['centrH_slit'] = np.average(ticket_slit['fwhm_coordinates_h'])
        design_df['centrV_slit'] = np.average(ticket_slit['fwhm_coordinates_v'])
        design_df['int'][i] = ticket['intensity']
        design_df['fwhmH'][i] = ticket['fwhm_h']
        design_df['fwhmV'][i] = ticket['fwhm_v']
        design_df['centrH'] = np.average(ticket['fwhm_coordinates_h'])
        design_df['centrV'] = np.average(ticket['fwhm_coordinates_v'])
        
        beams.append(beam)
        beams_slit.append(beam_slit)
    return design_df, beams, beams_slit
#
# main 
#
from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show


df = pd.read_csv('./design_table.csv')
out, beams, beams_slit = run_dataFrame(df)
# und = run_source(nrays = 1e4, ng_e=3)
# beam, beam_slit = run_beamline(und)