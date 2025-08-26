
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import sys
sys.path.insert(0,'n:/GitLab/phils-Xray-tools/XpyBru')

def oasys_wiggler_radiation(H_mesh=501, V_mesh=501, e_points=401, file='SW_2PA_neg.txt'):
    #
    # script to make the calculations (created by XOPPY:wiggler_radiation)
    #
    
    from xoppylib.sources.xoppy_bm_wiggler import xoppy_calc_wiggler_radiation
    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"]          = 6.0
    h5_parameters["ELECTRONCURRENT"]         = 0.2
    h5_parameters["PERIODID"]                = 0.12
    h5_parameters["NPERIODS"]                = 37.0
    h5_parameters["KV"]                      = 22.416
    h5_parameters["FIELD"]                   = 1   # 0= sinusoidal, 1=from file
    h5_parameters["FILE"]                    = 'C:/Users/brumund/Work Folders/Documents/Oasys/working_directory/' + file
    h5_parameters["POLARIZATION"]            = 0  # 0=total, 1=s, 2=p
    h5_parameters["DISTANCE"]                = 22.8
    h5_parameters["HSLITPOINTS"]             = H_mesh
    h5_parameters["VSLITPOINTS"]             = V_mesh
    h5_parameters["PHOTONENERGYMIN"]         = 100.0
    h5_parameters["PHOTONENERGYMAX"]         = 100100.0
    h5_parameters["PHOTONENERGYPOINTS"]      = e_points
    h5_parameters["SHIFT_X_FLAG"]            = 0
    h5_parameters["SHIFT_X_VALUE"]           = 0.0
    h5_parameters["SHIFT_BETAX_FLAG"]        = 0
    h5_parameters["SHIFT_BETAX_VALUE"]       = 0
    h5_parameters["CONVOLUTION"]             = 1
    h5_parameters["PASSEPARTOUT"]            = 3.0
    
    e, h, v, p, traj = xoppy_calc_wiggler_radiation(
            ELECTRONENERGY           = h5_parameters["ELECTRONENERGY"]         ,
            ELECTRONCURRENT          = h5_parameters["ELECTRONCURRENT"]        ,
            PERIODID                 = h5_parameters["PERIODID"]               ,
            NPERIODS                 = h5_parameters["NPERIODS"]               ,
            KV                       = h5_parameters["KV"]                     ,
            FIELD                    = h5_parameters["FIELD"]                  ,
            FILE                     = h5_parameters["FILE"]                   ,
            POLARIZATION             = h5_parameters["POLARIZATION"]           ,
            DISTANCE                 = h5_parameters["DISTANCE"]               ,
            HSLITPOINTS              = h5_parameters["HSLITPOINTS"]            ,
            VSLITPOINTS              = h5_parameters["VSLITPOINTS"]            ,
            PHOTONENERGYMIN          = h5_parameters["PHOTONENERGYMIN"]        ,
            PHOTONENERGYMAX          = h5_parameters["PHOTONENERGYMAX"]        ,
            PHOTONENERGYPOINTS       = h5_parameters["PHOTONENERGYPOINTS"]     ,
            SHIFT_X_FLAG             = h5_parameters["SHIFT_X_FLAG"]           ,
            SHIFT_X_VALUE            = h5_parameters["SHIFT_X_VALUE"]          ,
            SHIFT_BETAX_FLAG         = h5_parameters["SHIFT_BETAX_FLAG"]       ,
            SHIFT_BETAX_VALUE        = h5_parameters["SHIFT_BETAX_VALUE"]      ,
            CONVOLUTION              = h5_parameters["CONVOLUTION"]            ,
            PASSEPARTOUT             = h5_parameters["PASSEPARTOUT"]            ,
            h5_file                  = "wiggler_radiation.h5"                  ,
            h5_entry_name            = "XOPPY_RADIATION"                       ,
            h5_initialize            = False                                    ,
            h5_parameters            = h5_parameters                           ,
            )
    
    # h5_parameters = dict()
    # h5_parameters["ELECTRONENERGY"]          = 6.0
    # h5_parameters["ELECTRONCURRENT"]         = 0.2
    # h5_parameters["PERIODID"]                = 0.12
    # h5_parameters["NPERIODS"]                = 37.0
    # h5_parameters["KV"]                      = 22.416
    # h5_parameters["FIELD"]                   = 1   # 0= sinusoidal, 1=from file
    # h5_parameters["FILE"]                    = file
    # h5_parameters["POLARIZATION"]            = 0  # 0=total, 1=s, 2=p
    # h5_parameters["DISTANCE"]                = 23.0
    # h5_parameters["HSLITPOINTS"]             = H_mesh
    # h5_parameters["VSLITPOINTS"]             = V_mesh
    # h5_parameters["PHOTONENERGYMIN"]         = 50.0
    # h5_parameters["PHOTONENERGYMAX"]         = 100000.0
    # h5_parameters["PHOTONENERGYPOINTS"]      = e_points
    # h5_parameters["SHIFT_X_FLAG"]            = 0 # 4 for value at zero
    # h5_parameters["SHIFT_X_VALUE"]           = 0.0
    # h5_parameters["SHIFT_BETAX_FLAG"]        = 0 # 4 for value at zero
    # h5_parameters["SHIFT_BETAX_VALUE"]       = 0
    # h5_parameters["CONVOLUTION"]             = 1
    # h5_parameters["PASSEPARTOUT"]            = 3.0
    
    # e, h, v, p, traj = xoppy_calc_wiggler_radiation(
    #         ELECTRONENERGY           = h5_parameters["ELECTRONENERGY"]         ,
    #         ELECTRONCURRENT          = h5_parameters["ELECTRONCURRENT"]        ,
    #         PERIODID                 = h5_parameters["PERIODID"]               ,
    #         NPERIODS                 = h5_parameters["NPERIODS"]               ,
    #         KV                       = h5_parameters["KV"]                     ,
    #         FIELD                    = h5_parameters["FIELD"]                  ,
    #         FILE                     = h5_parameters["FILE"]                   ,
    #         POLARIZATION             = h5_parameters["POLARIZATION"]           ,
    #         DISTANCE                 = h5_parameters["DISTANCE"]               ,
    #         HSLITPOINTS              = h5_parameters["HSLITPOINTS"]            ,
    #         VSLITPOINTS              = h5_parameters["VSLITPOINTS"]            ,
    #         PHOTONENERGYMIN          = h5_parameters["PHOTONENERGYMIN"]        ,
    #         PHOTONENERGYMAX          = h5_parameters["PHOTONENERGYMAX"]        ,
    #         PHOTONENERGYPOINTS       = h5_parameters["PHOTONENERGYPOINTS"]     ,
    #         SHIFT_X_FLAG             = h5_parameters["SHIFT_X_FLAG"]           ,
    #         SHIFT_X_VALUE            = h5_parameters["SHIFT_X_VALUE"]          ,
    #         SHIFT_BETAX_FLAG         = h5_parameters["SHIFT_BETAX_FLAG"]       ,
    #         SHIFT_BETAX_VALUE        = h5_parameters["SHIFT_BETAX_VALUE"]      ,
    #         CONVOLUTION              = h5_parameters["CONVOLUTION"]            ,
    #         PASSEPARTOUT             = h5_parameters["PASSEPARTOUT"]            ,
    #         h5_file                  = ""                  ,
    #         h5_entry_name            = "XOPPY_RADIATION"                       ,
    #         h5_initialize            = False                                    ,
    #         h5_parameters            = h5_parameters                           ,
    #         )
    
    return p, h, v, e

def compress_3d_results_to_numpy(p, h, v, e):
    np.savez_compressed('./2PWA_d23m_3dradiation.npz', p=p, h=h, v=v, e=e)
    return


def read_data(file = './2PWA_d23m_3dradiation.npz', flip=True, h_corr=265):    
    print(f'Reading data: {file}')
    with np.load(file) as data:
        p, h, v, e = data['p'], data['h'], data['v'], data['e']
    print('file reading finished!')
    
    # Flipping p data along second axis (1: horizontal)
    if flip:
        p = np.flip(p,1)
        
    h = h+h_corr
    
    print('Calculating spectral power density and total power density')
    
    sp_pow  = p*1000*const.e            # Spectral power density (W/ev/mm^2)
    pd      = np.trapz(sp_pow, e, axis=0)    # Power density          (W/mm^2)
    
    print(f'Max. flux density: {p.max():.2e}ph/s/mm^2/0.1%bw, max. power density: {pd.max():.2f}W/mm^2')
    return p, h, v, e, sp_pow, pd

def get_filter_info(e, mat = 'Be', dens = '?', t = 0.5):
    from xoppylib.power.xoppy_calc_power import xoppy_calc_power
    import xraylib

    out = xoppy_calc_power(e, 
                np.ones(e.shape), 
                substance = [mat], 
                flags = [0], 
                dens = [dens], 
                nelements = 1,
                thick = [t],
                material_constants_library = xraylib)
    
    filter_abs = out['data'][1]-out['data'][-1]
    filter_trans = out['data'][-1]

    return filter_abs, filter_trans

def apply_absorption_filter(e, sp_pow, filter_abs):
    print('applying filter to each energy entry!')
    abs_3d = (sp_pow * filter_abs[:,None,None])
    
    print('done! Integrating over energy to obtain absorbed power density!')
    abs_2d = np.trapz(abs_3d, e, axis = 0)
    
    return abs_2d

def set_power_zero_inside(h, v, pd, hrange=(-3, 73), vrange=(-3, 3)):
    
    if (type(hrange) == tuple) and (len(hrange) == 2):
        hi1 = np.argmin(np.abs(h-hrange[0]))
        hi2 = np.argmin(np.abs(h-hrange[1]))+1
    else:
        print('Info none or no correct tuple provided for hrange!')
        hi1 = 0
        hi2 = -1

    if (type(vrange) == tuple) and (len(vrange) == 2):
        vi1 = np.argmin(np.abs(v-vrange[0]))
        vi2 = np.argmin(np.abs(v-vrange[1]))+1
    else:
        print('Info none or no correct tuple provided for vrange!')
        vi1 = 0
        vi2 = -1
     
    print(f'hi1: {hi1}, hi2: {hi2}, vi1: {vi1}, vi2: {vi2}')    
    pd[hi1:hi2,vi1:vi2] = np.zeros(pd[hi1:hi2,vi1:vi2].shape)
    
    return pd

def set_power_zero_outside(h, v, pd, hrange=(-3, 73), vrange=(-3, 3)):
    
    pd_neg = set_power_zero_inside(h, v, pd.copy(), hrange=hrange, vrange=vrange)
    pd -=  pd_neg 
    
    return pd

def get_total_power(h, v, pd):
    
    P = np.trapz(np.trapz(pd, v),h)
    
    return P


def create_inset_plot(ax, h, v, pd, pos=[0.5, -1, 0.5, 0.5], hrange=(-3, 73), vrange=(-3, 3)):
    
    ax_in = ax.inset_axes(pos)
    
    ax_in.imshow(pd.T, extent=[h[0],h[-1],v[0],v[-1]])

    ax_in.set_xlim(hrange[0], hrange[1])
    ax_in.set_ylim(vrange[0], vrange[1])
    ax_in.set_aspect('auto')
    ax.indicate_inset_zoom(ax_in)

    return ax_in


if __name__ == '__main__':
  
    # example plot
    from srxraylib.plot.gol import plot_image
    # import matplotlib.pyplot as plt
    
    savefig = True
    
    # p, h, v, e = oasys_wiggler_radiation()
    # compress_3d_results_to_numpy(p, h, v, e)
    # Read saved Wiggler data, flip and correct x data
    p, h, v, e, sp_pow, pd = read_data(file = './2PWA_d23m_3dradiation.npz', flip=False, h_corr=-105)
    
    # =============================================================================
    #   Plot wiggler power density in 23m  
    # =============================================================================
    P_tot = get_total_power(h, v, pd)
    fig, ax = plot_image(pd ,h ,v , 
                         title=f"Wiggler radiation, power density, Total P={P_tot:.2f}W", 
                         xtitle="H [mm]",
                         ytitle="V [mm]", aspect=1.5)
    
    create_inset_plot(ax, h, v, pd, pos=[0.35, -0.75, 0.5, 0.4])
    
    if savefig:
        fig.savefig('fig1-wiggler-pd-d23m.png', dpi = 300)
        # fig.savefig('fig1-wiggler-pd-d23m.svg')
    # =============================================================================
    #   Plot absorbed power density by support
    # =============================================================================
    pd_cu = set_power_zero_inside(h, v, pd.copy())
    
    P_tot = get_total_power(h, v, pd_cu)
    fig, ax = plot_image(pd_cu ,h ,v , 
                         title=f"Power density outside of Be window, Total P={P_tot:.2f}W", 
                         xtitle="H [mm]",
                         ytitle="V [mm]", aspect='auto')
    
    if savefig:
        fig.savefig('fig2-pd-onCu-d23m.png', dpi = 300)
        # fig.savefig('fig2-pd-onCu-d23m.svg')
    
    # =============================================================================
    #     Plot Absorbed power density by window
    # =============================================================================
    # # Beryllium
    # filter_abs, filter_trans = get_filter_info(e, t=0.5)
    # abs_2d_500um = apply_absorption_filter(e, sp_pow, filter_abs)
    # # abs_2d_500um = set_power_zero_outside(h, v, abs_2d_500um.copy())
    
    # P_tot = get_total_power(h, v, abs_2d_500um)
    # fig, ax = plot_image(abs_2d_500um,h,v,
    #                      title=f"Absorbed in t=0.50mm Be, Total P={P_tot:.2f}W",
    #                      xtitle="H [mm]",
    #                      ytitle="V [mm]", xrange=(-3, 73), yrange=(-3,3), aspect=4)
    # if savefig:
    #     fig.savefig('fig3-pd-absBe-500um-d23m.png', dpi = 300)
    #     # fig.savefig('fig3-pd-absBe-500um-d23m.svg')
        
    
    # filter_abs, filter_trans = get_filter_info(e, t=0.25)
    # abs_2d_250um = apply_absorption_filter(e, sp_pow, filter_abs)
    # # abs_2d_250um = set_power_zero_outside(h, v, abs_2d_250um.copy())
    
    # P_tot = get_total_power(h, v, abs_2d_250um)
    # fig, ax = plot_image(abs_2d_250um,h,v,
    #                      title=f"Absorbed in t=0.25mm Be, Total P={P_tot:.2f}W",
    #                      xtitle="H [mm]",
    #                      ytitle="V [mm]", xrange=(-3, 73), yrange=(-3,3), aspect=4)
    
    # if savefig:
    #     fig.savefig('fig4-pd-absBe-250um-d23m.png', dpi = 300)
    #     # fig.savefig('fig4-pd-absBe-250um-d23m.svg')
        
    
    # Boron Carbide  
    filter_abs, filter_trans = get_filter_info(e, mat = 'B4C', dens = 2.52, t=0.5)
    abs_2d_500um = apply_absorption_filter(e, sp_pow, filter_abs)
    abs_2d_500um = set_power_zero_outside(h, v, abs_2d_500um.copy())
    
    P_tot = get_total_power(h, v, abs_2d_500um)
    fig, ax = plot_image(abs_2d_500um,h,v,
                         title=f"Absorbed in t=0.50mm B4C, Total P={P_tot:.2f}W",
                         xtitle="H [mm]",
                         ytitle="V [mm]", xrange=(-3, 73), yrange=(-3,3), aspect=4)
    if savefig:
        fig.savefig('fig3-pd-absB4C-500um-d23m.png', dpi = 300)
        # fig.savefig('fig3-pd-absBe-500um-d23m.svg')
   
    filter_abs, filter_trans = get_filter_info(e, mat = 'B4C', dens = 2.52, t=0.25)
    abs_2d_250um = apply_absorption_filter(e, sp_pow, filter_abs)
    abs_2d_250um = set_power_zero_outside(h, v, abs_2d_250um.copy())
    
    P_tot = get_total_power(h, v, abs_2d_250um)
    fig, ax = plot_image(abs_2d_250um,h,v,
                         title=f"Absorbed in t=0.25mm B4C, Total P={P_tot:.2f}W",
                         xtitle="H [mm]",
                         ytitle="V [mm]", xrange=(-3, 73), yrange=(-3,3), aspect=4)
    
    if savefig:
        fig.savefig('fig4-pd-absB4C-250um-d23m.png', dpi = 300)
        # fig.savefig('fig4-pd-absBe-250um-d23m.svg')
    
    