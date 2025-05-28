from srxraylib.plot.gol import plot
from oasys.util.oasys_util import get_fwhm


data40 = in_object_1.get_contents('xoppy_data')
data80 = in_object_2.get_contents('xoppy_data')

fwhm40, _, _ = get_fwhm(data40[:,-1], data40[:,0])
fwhm80, _, _ = get_fwhm(data80[:,-1], data80[:,0])

plot(data40[:,0], data40[:,-1], data80[:,0], data80[:,-1], xtitle='$\Theta - \Theta_{Bragg} \, (\mu rad)$', ytitle='reflectivity (-)', legend = [f'40keV fwhm={fwhm40:.1f} urad', f'80keV fwhm={fwhm80:.1f} urad'])
