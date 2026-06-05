#
# script to make the calculations (created by XOPPY:undulator_spectrum)
#
from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_spectrum
from xoppylib.xoppy_run_binaries import xoppy_calc_ws
import numpy, h5py

def calc_power( h5name, name, KV, PERIODID, NPERIODS ):
    args = {
            'ELECTRONENERGY':6.0,
            'ELECTRONENERGYSPREAD':0.001,
            'ELECTRONCURRENT':0.2,
            # 'ELECTRONBEAMSIZEH':3.01836e-05,
            # 'ELECTRONBEAMSIZEV':3.63641e-06,
            # 'ELECTRONBEAMDIVERGENCEH':4.36821e-06,
            # 'ELECTRONBEAMDIVERGENCEV':1.37498e-06,
            "ELECTRONBEAMSIZEH" : 3.34281e-05,
            "ELECTRONBEAMSIZEV" : 7.28139e-06,
            "ELECTRONBEAMDIVERGENCEH" : 4.51097e-06,
            "ELECTRONBEAMDIVERGENCEV" : 1.94034e-06,
            'PERIODID':PERIODID,
            'NPERIODS':NPERIODS,
            'KV':KV,
            'KH':0.0,
            'KPHASE':0.0,
            'DISTANCE':23.0,
            'GAPH':0.0018,
            'GAPV':0.001,
            'GAPH_CENTER':0.0,
            'GAPV_CENTER':0.0,
            'PHOTONENERGYMIN':1000.0,
            'PHOTONENERGYMAX':200000.0,
            'PHOTONENERGYPOINTS':2000,
            'METHOD':2,
            'USEEMITTANCES':1
    }

    energy, flux, spectral_power, cumulated_power = xoppy_calc_undulator_spectrum( **args )
    with h5py.File( h5name, 'a') as hout:
            g = hout.require_group(name)
            g['energy'] = energy
            g['flux'] = flux
            g['spectral_power'] = spectral_power
            g['cumulated_power'] = cumulated_power
            g.attrs.update( args )

    args = {
        'ENERGY' : 6.0,
        'CUR'    : 200.0,
        'PERIOD' : PERIODID*100, # cm
        'N'      : NPERIODS,
        'KX'     : 0.0,
        'KY'     : KV,
        'EMIN'   : 1000.0,
        'EMAX'   : 200000.0,
        'NEE'    : 500,
        'D'      : 23.0,
        'XPC'    : 0.0,
        'YPC'    : 0.0,
        'XPS'    : 1.8,
        'YPS'    : 1.0,
        'NXP'    : 10,
        'NYP'    : 10,
    }

    out_file =  xoppy_calc_ws( ** args )

    # data to pass to power
    data = numpy.loadtxt(out_file)
    energy = data[:,0]             # eV
    flux = data[:,1]               # "Flux [photons/s/0.1%bw]"
    spectral_power = data[:,2]     # "Power [W/eV]"
    cumulated_power = data[:,3]

    with h5py.File( h5name, 'a') as hout:
        hout.require_group(name+"_WS")
        hout[name+"_WS"]['energy'] = energy
        hout[name+"_WS"]['flux'] = flux
        hout[name+"_WS"]['spectral_power'] = spectral_power
        hout[name+"_WS"]['cumulated_power'] = cumulated_power
        hout[name+"_WS"].attrs.update( args )




if __name__ == "__main__":

    import sys
    sys.path.insert(0,'.')
    # import calcspectra, numpy

    # srio: cpmu20 of id16 has K=2.334 at 6mm gap
    IDS = {'cpmu18': { 'KV' : 1.6563, 'gap_min' : 6.0,  'nperiods' : 111, 'period' : 0.018  },
           'u22'   : { 'KV' : 1.5426, 'gap_min' : 6.8,  'nperiods' : 91, 'period' : 0.022 },
           'cpmu20': {'KV': 2.334, 'gap_min': 6.0, 'nperiods': 98, 'period': 0.0205}}
           # 'cpmu20': { 'KV' : 2.2338, 'gap_min' : 6.0,  'nperiods' : 98, 'period' : 0.0205 } }

    h5name = 'WattDog.h5'

    for device in IDS:
        KV = IDS[device]['KV']
        NPERIODS = IDS[device]['nperiods']
        PERIODID = IDS[device]['period']
        groupname = f'{device}_{KV}'
        calc_power( h5name, groupname, KV, PERIODID, NPERIODS )
        for KV in numpy.arange( 0.5, int(KV*10)/10, 0.1):
            groupname = f'{device}_{KV:.1f}'
            calc_power( h5name, groupname, KV, PERIODID, NPERIODS )


    #
    # attenuators
    #

    import collections
    import pprint

    import numpy
    from xoppylib.power.xoppy_calc_power import xoppy_calc_power
    import xraylib
    from dabax.dabax_xraylib import DabaxXraylib
    import h5py


    def flatten(xss):
        return [x for xs in xss for x in xs]


    def sumdisks(l):
        totals = collections.defaultdict(int)
        for element, thickness in l:
            totals[element] += thickness
        return totals


    cvd = ('C', 0.3)
    D = cvd,
    Cu = ('Cu', 0.0025), cvd, ('Cu', 0.0025)
    d = ('Au', 0.0038), cvd, ('Au', 0.0038)
    c = ('Au', 0.0013), cvd, ('Au', 0.0038)
    Al1 = ('Al', 1.0),

    attenuators = {'2026': flatten((D, Cu, Cu, d, c, c, Al1, Al1))}

    ## Post Refurbishment attenuators
    cvd5 = ('C', 0.5)

    ax1 = cvd5, ('Cu', 0.0015), cvd, ('Cu', 0.0015)
    ax2 = ('Cu', 0.0025), cvd, ('Cu', 0.0025), ('Cu', 0.005), cvd, ('Cu', 0.005)
    ax3 = ('Cu', 0.0075), cvd, ('Cu', 0.0075), ('Cu', 0.010), cvd, ('Cu', 0.010)
    ax4 = ('Cu', 0.0150), cvd, ('Cu', 0.0150), ('Ag', 0.005), cvd, ('Ag', 0.005)
    ax5 = ('Cu', 0.0075), cvd, ('Cu', 0.0075), ('Ag', 0.010), cvd, ('Ag', 0.010)
    ax6 = ('Cu', 0.0150), cvd, ('Cu', 0.0150), ('Ag', 0.025), cvd, ('Ag', 0.025)

    axes = [ax1, ax2, ax3, ax4, ax5, ax6]

    density = {'C': 3.52,
               'Cu': 8.96,
               'Au': 19.32,
               'Ag': 10.49,
               'Al': 2.70,
               'Si': 2.33}

    for i in range(6):
        disks = [axes[j] for j in range(i)]
        attenuators[f'2028.{i}'] = flatten(disks)

    with h5py.File('WattDog.h5', 'r+') as hin:
        for name in list(hin):
            grp = hin[name]
            energy = grp['energy'][:]
            spectral_power = grp['spectral_power'][:]

            for k in attenuators:
                # print(k, attenuators[k])
                # print(k, sumdisks( attenuators[k] ) )
                att = attenuators[k]
                sum_att = sumdisks(att)
                print(sum_att)
                substance = [element for element in sum_att]
                thickness = [sum_att[element] for element in substance]
                dens = [density[element] for element in substance]
                n = len(substance) + 2  # two extra silicon behind for mono
                out_dictionary = xoppy_calc_power(
                    energy.copy(),
                    spectral_power.copy(),
                    substance=substance + ['Si', 'Si'],
                    thick=thickness + [2.5, 90],  # LL and Bragg mono's
                    angle=[3, ] * n,  # in mrad (for mirrors)
                    dens=dens + [density['Si'], density['Si']],
                    roughness=[0, ] * n,  # in A (for mirrors)
                    flags=[0, ] * n,  # 0=Filter, 1=Mirror
                    nelements=n,
                    FILE_DUMP=0,
                    material_constants_library=xraylib,
                )
                #            print('Here! input power', spectral_power.sum() * (energy[1]-energy[0]) )
                #            print(out_dictionary['info'])
                #            1/0
                gout = grp.require_group(k)
                for k in out_dictionary:
                    if k in gout:
                        gout[k][()] = out_dictionary[k]
                    else:
                        gout[k] = out_dictionary[k]


