# -*- coding: utf-8 -*-
"""
Created on Thu May 15 16:35:21 2025

@author: BRUMUND
"""


#
# import undulator radiation from file
#
import contextlib
import h5py
import numpy as np
import scipy.constants as con
from xoppylib.power.power3d import calculate_component_absorbance_and_transmittance
from xoppylib.power.power3d import apply_transmittance_to_incident_beam



def read_file(file = 'merged-cpmu18.h5'):
    print('Reading file:', file)
    path = 'C:/Users/brumund/Work Folders/Documents/Python/OASYS-und-rad-multiproc-id11/'
    fname = path + file
    code = 'SRW'
    hf = h5py.File('C:/Users/brumund/Work Folders/Documents/Python/OASYS-und-rad-multiproc-id11/merged-cpmu18.h5','r')
    flux3D = hf["/XOPPY_RADIATION/Radiation/stack_data"][:]
    energy = hf["/XOPPY_RADIATION/Radiation/axis0"][:]
    horizontal = hf["/XOPPY_RADIATION/Radiation/axis1"][:]
    vertical = hf["/XOPPY_RADIATION/Radiation/axis2"][:]
    hf.close()
    
    return energy, horizontal, vertical, flux3D

def run_bl(energy, horizontal, vertical, flux3D, mat = 'Be', n_lenses = 1):
    # TODO: Add Al lenses (if Al, use 60 Be upstream)
    
    # compute local transmittance and absorbance
    e0, h0, v0, f0  = energy, horizontal, vertical, flux3D
    transmittance, absorbance, E, H, V, txt = calculate_component_absorbance_and_transmittance(
                    e0, # energy in eV
                    h0, # h in mm
                    v0, # v in mm
                    substance='C',
                    thick=0.3,
                    angle=3.0,
                    defection=1,
                    dens='3.51',
                    roughness=0.0,
                    flags=0, # 0=Filter 1=Mirror 2=Aperture 3=magnifier, 4=Screen rotation  5=Thin object  6=Multilayer 7=External file
                    gapshape=1,     # 0=rectangle, 1=circle  # used if flag in (0, 1, 2, 6, 7)
                    hgap=0.8,
                    vgap=0.8,
                    hgapcenter=0.0,
                    vgapcenter=0.0,
                    hmag=1.0,
                    vmag=1.0,
                    hrot=0.0,
                    vrot=0.0,
                    thin_object_file='',
                    thin_object_thickness_outside_file_area=0.0,
                    thin_object_back_profile_flag=0,
                    thin_object_back_profile_file='',
                    multilayer_file='',
                    external_reflectivity_file='',
                    )
    
    # apply transmittance to incident beam 
    f_transmitted, e, h, v = apply_transmittance_to_incident_beam(transmittance, f0, e0, h0, v0,
                    flags = 0,
                    hgap = 0.8,
                    vgap = 0.8,
                    hgapcenter = 0.0,
                    vgapcenter = 0.0,
                    hmag = 1.0,
                    vmag = 1.0,
                    interpolation_flag     = 0,
                    interpolation_factor_h = 1.0,
                    interpolation_factor_v = 1.0,
                    slit_crop = 0,
                    )
    
    f_absorbed = f0 * absorbance / (H[0] / h0[0]) / (V[0] / v0[0])
    
    # data to pass
    energy, horizontal, vertical, flux3D = e, h, v, f_transmitted.copy()
    
    P_lenses = []
    
    # if Al run 60 Be lenses
    if mat == 'Al':
        # compute local transmittance and absorbance
        e0, h0, v0, f0  = energy, horizontal, vertical, flux3D
        transmittance, absorbance, E, H, V, txt = calculate_component_absorbance_and_transmittance(
                        e0, # energy in eV
                        h0, # h in mm
                        v0, # v in mm
                        substance='Be',
                        thick=0.6,
                        angle=3.0,
                        defection=1,
                        dens='?',
                        roughness=0.0,
                        flags=5, # 0=Filter 1=Mirror 2=Aperture 3=magnifier, 4=Screen rotation  5=Thin object  6=Multilayer 7=External file
                        hgap=0.8,
                        vgap=0.8,
                        hgapcenter=0.0,
                        vgapcenter=0.0,
                        hmag=1.0,
                        vmag=1.0,
                        hrot=0.0,
                        vrot=0.0,
                        thin_object_file=r'c:\Users\brumund\Work Folders\Documents\Oasys\working_directory\lens-60-Be.h5',
                        thin_object_thickness_outside_file_area=0.0,
                        thin_object_back_profile_flag=0,
                        thin_object_back_profile_file='<none>',
                        multilayer_file='',
                        external_reflectivity_file='',
                        )
        
        # apply transmittance to incident beam 
        f_transmitted, e, h, v = apply_transmittance_to_incident_beam(transmittance, f0, e0, h0, v0,
                        flags = 5,
                        hgap = 0.8,
                        vgap = 0.8,
                        hgapcenter = 0.0,
                        vgapcenter = 0.0,
                        hmag = 1.0,
                        vmag = 1.0,
                        interpolation_flag     = 0,
                        interpolation_factor_h = 1.0,
                        interpolation_factor_v = 1.0,
                        slit_crop = 0,
                        )
        
        f_absorbed = f0 * absorbance / (H[0] / h0[0]) / (V[0] / v0[0])
        
        # data to pass
        energy, horizontal, vertical, flux3D = e, h, v, f_transmitted
    
    for i in range(n_lenses):
        print('Run material ', mat,' lens i=', i)
        with contextlib.redirect_stdout(None):
            # compute local transmittance and absorbance
            e0, h0, v0, f0  = energy, horizontal, vertical, flux3D
            transmittance, absorbance, E, H, V, txt = calculate_component_absorbance_and_transmittance(
                            e0, # energy in eV
                            h0, # h in mm
                            v0, # v in mm
                            substance=mat,
                            thick=0.6,
                            angle=3.0,
                            defection=1,
                            dens='?',
                            roughness=0.0,
                            flags=5, # 0=Filter 1=Mirror 2=Aperture 3=magnifier, 4=Screen rotation  5=Thin object  6=Multilayer 7=External file
                            hgap=0.8,
                            vgap=0.8,
                            hgapcenter=0.0,
                            vgapcenter=0.0,
                            hmag=1.0,
                            vmag=1.0,
                            hrot=0.0,
                            vrot=0.0,
                            thin_object_file=r'c:\Users\brumund\Work Folders\Documents\Oasys\working_directory\lens-1-Be.h5',
                            thin_object_thickness_outside_file_area=0.0,
                            thin_object_back_profile_flag=0,
                            thin_object_back_profile_file='<none>',
                            multilayer_file='',
                            external_reflectivity_file='',
                            )
            
            # apply transmittance to incident beam 
            f_transmitted, e, h, v = apply_transmittance_to_incident_beam(transmittance, f0, e0, h0, v0,
                            flags = 5,
                            hgap = 0.8,
                            vgap = 0.8,
                            hgapcenter = 0.0,
                            vgapcenter = 0.0,
                            hmag = 1.0,
                            vmag = 1.0,
                            interpolation_flag     = 0,
                            interpolation_factor_h = 1.0,
                            interpolation_factor_v = 1.0,
                            slit_crop = 0,
                            )
        
        f_absorbed = f0 * absorbance / (H[0] / h0[0]) / (V[0] / v0[0])
        
        # data to pass
        energy, horizontal, vertical, flux3D = e, h, v, f_transmitted
    
        P_abs = np.trapz(np.trapz(np.trapz(f_absorbed,v),h),e)*con.e*1000
        P_lenses.append(P_abs)
        
    return np.array(P_lenses)

if __name__ == '__main__':
    from XpyBru.PlotTools import pol
    
    energy, horizontal, vertical, flux3D = read_file('merged-cpmu18.h5')
    P_Be_cpmu18 = run_bl(energy, horizontal, vertical, flux3D, n_lenses = 16)
    P_Al_cpmu18 = run_bl(energy, horizontal, vertical, flux3D, n_lenses = 16, mat='Al')
    
    
    energy, horizontal, vertical, flux3D = read_file('merged-cpmu20.h5')
    P_Be_cpmu20 = run_bl(energy, horizontal, vertical, flux3D, n_lenses = 16)
    P_Al_cpmu20 = run_bl(energy, horizontal, vertical, flux3D, n_lenses = 16, mat='Al')
    
    P_Be = P_Be_cpmu18 + P_Be_cpmu20
    P_Al = P_Al_cpmu18 + P_Al_cpmu20
    
    i = list(range(16))
    fig, ax = pol([i,i],[P_Be,P_Al], x_label='Lens i (-)', y_label='Abs. power (W)', 
        legend=[f'Beryllium (P_tot={sum(P_Be):.1f} W)', f'Aluminum (P_tot={sum(P_Al):.1f} W)'])
    