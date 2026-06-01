import h5py

with h5py.File('WattDog.h5','r') as hin:
    data = {}
    labels = {}
    info = {}
    energy = {}
    spectral_power = {}
#    hin.visititems(print)
    for src in list(hin):
        energy[src] = hin[src]['energy'][()]
        spectral_power[src] = hin[src]['spectral_power'][()]
        data[src] = {}
        labels[src] = {}
        info[src] = {}
        for name in list(hin[src]):
            grp = hin[src][name]
            if isinstance( grp, h5py.Group ):
                data[src][name] = grp['data'][:]
                labels[src][name] = grp['labels'][()]
                info[src][name] = grp['info'][()]


# re-arrange into dicts of dicts
results = {}
for name in energy:
    items = name.split("_")
    method = 'SRW'
    if len(items)==2:
        device, K = items
    else:
        device, K, method = items
    for att in list( data[name] ):
        if device not in results:
            results[device] = {}
        if method not in results[device]:
            results[device][method] = {}
        if K not in results[device][method]:
            results[device][method][K] = {}
        results[device][method][K][att] = data[name][att]

import pylab as plt

f, a = plt.subplots(3, 7, constrained_layout=True, figsize=(15, 8))

for j, device in enumerate(('cpmu18', 'cpmu20', 'u22')):
    Kvs = list(results[device]['SRW'])
    attenuators = list(results[device]['SRW'][Kvs[0]])
    attenuators.append(attenuators.pop(0))  # 2026 last
    for i, attenuator in enumerate(attenuators):
        ll = []
        bb = []
        for K in Kvs:
            data = results[device]['SRW'][K][attenuator]
            e = data[0]
            de = e[1] - e[0]
            ll.append(data[-8].sum() * de)
            bb.append((data[-8].sum() + data[-2].sum()) * de)
        llw = []
        bbw = []
        for K in Kvs:
            data = results[device]['WS'][K][attenuator]
            e = data[0]
            de = e[1] - e[0]
            llw.append(data[-8].sum() * de)
            bbw.append((data[-8].sum() + data[-2].sum()) * de)
        a[j, i].plot([float(k) for k in Kvs], ll, "o-", label='ll srw')
        a[j, i].plot([float(k) for k in Kvs], llw, "o-", label='ll ws')
        a[j, i].plot([float(k) for k in Kvs], bb, "o-", label='bb srw')
        a[j, i].plot([float(k) for k in Kvs], bbw, "o-", label='bb ws')
        a[j, i].set(title=f'{device} {attenuator}', ylim=(0, None), xlim=(0, None))
        a[j, i].legend(loc='upper left')

for i in range(a.shape[0]):
    a[i, 0].set(ylabel='Power/W')
for i in range(a.shape[1]):
    a[-1, i].set(xlabel='K')

plt.show()



# https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z14.html
# Energy       μ/ρ        μen/ρ
#      (MeV)      (cm2/g)     (cm2/g)
#____________________________________
#
import numpy as np
nist = np.transpose([[float(x) for x in line.split()] for line in """1.00000E-03  1.570E+03  1.567E+03 
   1.50000E-03  5.355E+02  5.331E+02 
   1.83890E-03  3.092E+02  3.070E+02 
   1.83890E-03  3.192E+03  3.059E+03 
   2.00000E-03  2.777E+03  2.669E+03 
   3.00000E-03  9.784E+02  9.516E+02 
   4.00000E-03  4.529E+02  4.427E+02 
   5.00000E-03  2.450E+02  2.400E+02 
   6.00000E-03  1.470E+02  1.439E+02 
   8.00000E-03  6.468E+01  6.313E+01 
   1.00000E-02  3.389E+01  3.289E+01 
   1.50000E-02  1.034E+01  9.794E+00 
   2.00000E-02  4.464E+00  4.076E+00 
   3.00000E-02  1.436E+00  1.164E+00 
   4.00000E-02  7.012E-01  4.782E-01 
   5.00000E-02  4.385E-01  2.430E-01 
   6.00000E-02  3.207E-01  1.434E-01 
   8.00000E-02  2.228E-01  6.896E-02 
   1.00000E-01  1.835E-01  4.513E-02 
   1.50000E-01  1.448E-01  3.086E-02 
   2.00000E-01  1.275E-01  2.905E-02 
   3.00000E-01  1.082E-01  2.932E-02 
   4.00000E-01  9.614E-02  2.968E-02 
   5.00000E-01  8.748E-02  2.971E-02 
   6.00000E-01  8.077E-02  2.951E-02 
   8.00000E-01  7.082E-02  2.875E-02 
   1.00000E+00  6.361E-02  2.778E-02 
   1.25000E+00  5.688E-02  2.652E-02 
   1.50000E+00  5.183E-02  2.535E-02 
   2.00000E+00  4.480E-02  2.345E-02 
   3.00000E+00  3.678E-02  2.101E-02 
   4.00000E+00  3.240E-02  1.963E-02 
   5.00000E+00  2.967E-02  1.878E-02 
   6.00000E+00  2.788E-02  1.827E-02 
   8.00000E+00  2.574E-02  1.773E-02 
   1.00000E+01  2.462E-02  1.753E-02 
   1.50000E+01  2.352E-02  1.746E-02 
   2.00000E+01  2.338E-02  1.757E-02 """.split('\n')])

f,a = plt.subplots(1,2)
a[0].plot( nist[0], nist[1], '-', label='mu')
a[0].plot( nist[0], nist[2], '-', label='mu_en')
a[0].legend()
a[0].set( xscale='log', yscale='log')
a[1].plot( nist[0]*1e3, nist[2]/nist[1], '-', label='mu_en/mu')
a[1].set( xscale='log', title='Fraction')
a[1].legend()

plt.show()


import xraylib

Z_Si = xraylib.SymbolToAtomicNumber( 'Si' )

def nistcalc( data_in, plot=True ):
    data = data_in[:-6] # snip off bragg
    dE = (data[0][1:]-data[0][:-1]).mean() # constant step in energy
    E_keV = data[0]/1e3
    E_MeV = E_keV /1e3
    mu_nist = np.interp( E_MeV, nist[0], nist[1] )   # might not be optimal interpolation scheme (log log?)
    muen_nist = np.interp( E_MeV, nist[0], nist[2] )

    if plot:
        f, a = plt.subplots(1, 2,  figsize=(16,5), constrained_layout=True)
        phot = [ xraylib.CS_Photo(Z_Si, E) for E in E_keV ]
        comp = [ xraylib.CS_Compt(Z_Si, E) for E in E_keV ]
        rayl = [ xraylib.CS_Rayl(Z_Si, E) for E in E_keV ]
        total = [ xraylib.CS_Total(Z_Si, E) for E in E_keV ]
        a[0].plot( E_keV, data[-7], 'r-', label='Incident on crystal')
        a[0].plot( E_keV, data[-1], 'g-', label='Transmitted')
        a[0].plot( E_keV, data[-2], 'k-', label='Power lost in crystal (CS Total)')
        a[0].legend()
        a[0].set( ylabel='Power W/eV', xlabel='E/keV', title='Cumulative power through 2.5mm Si')

        a[1].plot( E_keV, np.cumsum(data[-2])*dE, 'k-', label='CS Total')
        a[1].stackplot( E_keV, ( np.cumsum(data[-2]*phot/total)*dE, np.cumsum(data[-2]*comp/total)*dE, np.cumsum(data[-2]*rayl/total)*dE),
               labels = ('CS Photo', 'CS Compton', 'CS Rayleigh'),
               colors= ('gray', 'skyblue', 'blue'))
        a[1].set( ylabel='Cumulative Power in 2.5mm Si / W',
          xlabel='E/keV',
          title='Xraylib: losses in crystal, 2.5mm Si')

        a[1].plot( E_keV, np.cumsum(data[-2]*muen_nist/mu_nist)*dE, "m-", label='NIST mass energy abs. (mu_en/mu)')
        a[1].legend()

    e_deposited = np.sum(data[-2]*muen_nist/mu_nist)*dE
    return e_deposited

nistcalc( results['cpmu18']['SRW']['1.6563']['2028.5'] )

plt.show()

nistcalc( results['cpmu18']['WS']['1.6563']['2028.5'] )

plt.show()

f, a = plt.subplots( 3, 7, constrained_layout=True, figsize=(15,8))

for j, device in enumerate( ('cpmu18' , 'cpmu20', 'u22' )):
    Kvs = list(results[device]['SRW'])
    attenuators = list(results[device]['SRW'][Kvs[0]])
    attenuators.append( attenuators.pop(0) ) # 2026 last
    for i, attenuator in enumerate( attenuators ):
        ll = []
        lln = []
        for K in Kvs:
            data = results[device]['SRW'][K][attenuator]
            e = data[0]
            de = e[1]-e[0]
            ll.append( data[-8].sum() * de )
            lln.append( nistcalc( results[device]['SRW'][K][attenuator], plot=False))
        llw = []
        llwn = []
        for K in Kvs:
                data = results[device]['WS'][K][attenuator]
                e = data[0]
                de = e[1]-e[0]
                llw.append( data[-8].sum() * de )
                llwn.append( nistcalc( results[device]['WS'][K][attenuator], plot=False))
        a[j,i].plot( [ float(k) for k in Kvs ], ll, "o-", label='ll srw' )
        a[j,i].plot( [ float(k) for k in Kvs ], llw, "o-", label='ll ws' )
        a[j,i].plot( [ float(k) for k in Kvs ], lln, "o-", label='ll nist srw' )
        a[j,i].plot( [ float(k) for k in Kvs ], llwn, "o-", label='ll nist ws' )
        a[j,i].set( title=f'{device} {attenuator}', ylim=(0,None), xlim=(0,None))
        a[j,i].legend(loc='upper left')

for i in range(a.shape[0]):
    a[i,0].set(ylabel='Power/W')
for i in range(a.shape[1]):
    a[-1,i].set(xlabel='K')


plt.show()