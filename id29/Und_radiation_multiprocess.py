# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:28:10 2024

@author: BRUMUND
"""

from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_radiation
import os
import numpy as np
import multiprocessing

def calc_und_rad(K = 2.277, n = 10, E = 10):
    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"]          = 6.0
    h5_parameters["ELECTRONENERGYSPREAD"]    = 0.001
    h5_parameters["ELECTRONCURRENT"]         = 0.2
    h5_parameters["ELECTRONBEAMSIZEH"]       = 3.01836e-05
    h5_parameters["ELECTRONBEAMSIZEV"]       = 3.63641e-06
    h5_parameters["ELECTRONBEAMDIVERGENCEH"] = 4.36821e-06
    h5_parameters["ELECTRONBEAMDIVERGENCEV"] = 1.37498e-06
    h5_parameters["PERIODID"]                = 0.0164
    h5_parameters["NPERIODS"]                = 121.951
    h5_parameters["KV"]                      = K
    h5_parameters["KH"]                      = 0.0
    h5_parameters["KPHASE"]                  = 0.0
    h5_parameters["DISTANCE"]                = 31.5
    h5_parameters["SETRESONANCE"]            = 0
    h5_parameters["HARMONICNUMBER"]          = 1
    h5_parameters["GAPH"]                    = 0.0021
    h5_parameters["GAPV"]                    = 0.0012
    h5_parameters["HSLITPOINTS"]             = 101
    h5_parameters["VSLITPOINTS"]             = 51
    h5_parameters["METHOD"]                  = 2
    h5_parameters["PHOTONENERGYMIN"]         = 8000.0
    h5_parameters["PHOTONENERGYMAX"]         = 105000.0
    h5_parameters["PHOTONENERGYPOINTS"]      = n
    h5_parameters["USEEMITTANCES"]           = 1
            
    path = f'./und_rad-K{K:5.3f}-E{E:4.1f}keV'
    try:
        os.mkdir(path)
        os.chdir(path)
    except OSError as error:
        print(error, 'changing directory')
        os.chdir(path)
    
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
            h5_file                  = "undulator_radiation.h5",
            h5_entry_name            = "XOPPY_RADIATION",
            h5_initialize            = True,
            h5_parameters            = h5_parameters, 
            )
    os.rename('undulator_radiation.h5', f'und_rad-K{K:5.3f}-n{n:.0f}-H101-V51.h5')
    os.remove('undulator_radiation.spec')
    pass

if __name__ == "__main__":
    # printing main program process id
    
    mp_context = multiprocessing.get_context('spawn')
    
    print("ID of main process: {}".format(os.getpid()))
    # Ks=np.linspace(0.1,1,10)
    Ks=[1.46,1.325,1.2,0.97,0.76,0.24]
    Es = [10,11,12,14,16,20]
    # 
    # Ks=[1.46,1.325]
    # Es = [10,11]
    ps = []
    
    for i, (Ki,Ei) in enumerate(zip(Ks,Es)):
        ps.append(mp_context.Process(target=calc_und_rad, args=(Ki,4851,Ei)))
        ps[i].start()
        
