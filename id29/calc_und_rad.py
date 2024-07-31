import numpy as np
import scipy.constants as codata
from syned.util.json_tools import load_from_json_file
from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_radiation

def run_calcualtions(syned_json_file, k_value, e_min=1000, e_max=50000, e_step=5,
                    slit_distance=30, h_aperture=0.001, v_aperture=0.001,
                    harmonic=1, use_emittances=True, save_file=True, prefix='Test'):
    
    #load syned JSON
    
    syned_obj = load_from_json_file(syned_json_file)
    
    e_beam = syned_obj.get_electron_beam() 
    
    (sigma_x, sig_div_x, sigma_y, sig_div_y) = e_beam.get_sigmas_all()
    
    und = syned_obj.get_magnetic_structure()

    h5_name = f"{prefix}_und_rad_{syned_json_file.replace('.json','')}.h5"
    
    #k_value = und.get_K_from_photon_energy(e_photon_res, e_beam.gamma(),
    #                                            harmonic=harmonic)	    
    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"]          = e_beam.energy()
    h5_parameters["ELECTRONENERGYSPREAD"]    = e_beam._energy_spread
    h5_parameters["ELECTRONCURRENT"]         = e_beam.current()
    h5_parameters["ELECTRONBEAMSIZEH"]       = sigma_x
    h5_parameters["ELECTRONBEAMSIZEV"]       = sigma_y
    h5_parameters["ELECTRONBEAMDIVERGENCEH"] = sig_div_x
    h5_parameters["ELECTRONBEAMDIVERGENCEV"] = sig_div_y
    h5_parameters["PERIODID"]                = und.period_length()
    h5_parameters["NPERIODS"]                = und.number_of_periods()
    h5_parameters["KV"]                      = k_value
    h5_parameters["KH"]                      = 0.0
    h5_parameters["KPHASE"]                  = 0.0
    h5_parameters["DISTANCE"]                = slit_distance
    h5_parameters["SETRESONANCE"]            = 0
    h5_parameters["HARMONICNUMBER"]          = harmonic
    h5_parameters["GAPH"]                    = h_aperture
    h5_parameters["GAPV"]                    = v_aperture
    h5_parameters["HSLITPOINTS"]             = 151
    h5_parameters["VSLITPOINTS"]             = 151
    h5_parameters["METHOD"]                  = 2 #SRW
    h5_parameters["PHOTONENERGYMIN"]         = e_min
    h5_parameters["PHOTONENERGYMAX"]         = e_max
    h5_parameters["PHOTONENERGYPOINTS"]      = int((e_max - e_min)/e_step)
    h5_parameters["USEEMITTANCES"]           = use_emittances

    energy, horizontal, vertical, flux3D, code = xoppy_calc_undulator_radiation(
        ELECTRONENERGY           = h5_parameters["ELECTRONENERGY"]         ,
        ELECTRONENERGYSPREAD     = h5_parameters["ELECTRONENERGYSPREAD"]   ,
        ELECTRONCURRENT          = h5_parameters["ELECTRONCURRENT"]        ,
        ELECTRONBEAMSIZEH        = h5_parameters["ELECTRONBEAMSIZEH"]      ,
        ELECTRONBEAMSIZEV        = h5_parameters["ELECTRONBEAMSIZEV"]      ,
        ELECTRONBEAMDIVERGENCEH  = h5_parameters["ELECTRONBEAMDIVERGENCEH"],
        ELECTRONBEAMDIVERGENCEV  = h5_parameters["ELECTRONBEAMDIVERGENCEV"],
        PERIODID                 = h5_parameters["PERIODID"]               ,
        NPERIODS                 = h5_parameters["NPERIODS"]               ,
        KV                       = h5_parameters["KV"]                     ,
        KH                       = h5_parameters["KH"]                     ,
        KPHASE                   = h5_parameters["KPHASE"]                 ,
        DISTANCE                 = h5_parameters["DISTANCE"]               ,
        SETRESONANCE             = h5_parameters["SETRESONANCE"]           ,
        HARMONICNUMBER           = h5_parameters["HARMONICNUMBER"]         ,
        GAPH                     = h5_parameters["GAPH"]                   ,
        GAPV                     = h5_parameters["GAPV"]                   ,
        HSLITPOINTS              = h5_parameters["HSLITPOINTS"]            ,
        VSLITPOINTS              = h5_parameters["VSLITPOINTS"]            ,
        METHOD                   = h5_parameters["METHOD"]                 ,
        PHOTONENERGYMIN          = h5_parameters["PHOTONENERGYMIN"]        ,
        PHOTONENERGYMAX          = h5_parameters["PHOTONENERGYMAX"]        ,
        PHOTONENERGYPOINTS       = h5_parameters["PHOTONENERGYPOINTS"]     ,
        USEEMITTANCES            = h5_parameters["USEEMITTANCES"]          ,
        h5_file                  = h5_name,
        h5_entry_name            = "XOPPY_RADIATION",
        h5_initialize            = True,
        h5_parameters            = h5_parameters, 
        )    
    

if __name__=="__main__":
    #pass
    run_calcualtions('ESRF_ID29_EBS_CPMU16_4.json', 1.255, e_min=9000, e_max=100000, e_step=5,
                    slit_distance=30, h_aperture=0.002, v_aperture=0.001,
                    harmonic=1, use_emittances=True, save_file=True, prefix='11.56keV_test')