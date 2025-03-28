import numpy
from srxraylib.plot.gol import plot, plot_show
import scipy.constants as codata

def calc_T(theta1=0.2):
    from xoppylib.mlayer import MLayer

    out = MLayer.initialize_from_bilayer_stack(
        material_S="Si", density_S=2.33, roughness_S=0.0,
        material_E="C", density_E=2.2, roughness_E=4.0,
        material_O="Ti", density_O=4.55, roughness_O=4.0,
        bilayer_pairs=10,
        bilayer_thickness=173.0,
        bilayer_gamma=0.6879,
    )

    for key in out.pre_mlayer_dict.keys():
        print(key, out.pre_mlayer_dict[key])
    # reflectivity is for amplitude
    rs, rp, e, t = out.scan(h5file="multilayerTiC.h5",
                            energyN=25000, energy1=1000.0, energy2=300000.0,
                            thetaN=1, theta1=theta1, theta2=6.0)
    return e, rp[:, 0] ** 2

def calc_P(theta1=0.2):
    from xoppylib.mlayer import MLayer

    out = MLayer.initialize_from_bilayer_stack(
        material_S="Si", density_S=2.33, roughness_S=0.0,
        material_E="C", density_E=2.2, roughness_E=3.0,
        material_O="Pd", density_O=10.7, roughness_O=3.0,
        bilayer_pairs=15,
        bilayer_thickness=89.0,
        bilayer_gamma=0.5618,
    )

    for key in out.pre_mlayer_dict.keys():
        print(key, out.pre_mlayer_dict[key])
    # reflectivity is for amplitude
    rs, rp, e, t = out.scan(h5file="",
                            energyN=5000, energy1=1000.0, energy2=100000.0,
                            thetaN=1, theta1=theta1, theta2=6.0)

    return e, rp[:, 0] ** 2

def calc_W(theta1=0.2):
    from xoppylib.mlayer import MLayer

    out = MLayer.initialize_from_bilayer_stack(
        material_S="Si", density_S=2.33, roughness_S=0.0,
        material_E="C", density_E=2.9, roughness_E=3.0,
        material_O="W", density_O=16.6, roughness_O=3.0,
        bilayer_pairs=30,
        bilayer_thickness=48.0,
        bilayer_gamma=0.5833,
    )

    for key in out.pre_mlayer_dict.keys():
        print(key, out.pre_mlayer_dict[key])
    # reflectivity is for amplitude
    rs, rp, e, t = out.scan(h5file="",
                            energyN=1000, energy1=5000.0, energy2=100000.0,
                            thetaN=1, theta1=theta1, theta2=6.0)

    return e, rp[:, 0] ** 2

Timd = numpy.loadtxt('/users/srio/Ti_C.txt', skiprows=16)

Tp2 = numpy.loadtxt('../MORAWE/ID18_DMM_Flux_TiC_E_0.2deg.txt', skiprows=1)
Tp4 = numpy.loadtxt('../MORAWE/ID18_DMM_Flux_TiC_E_0.4deg.txt', skiprows=1)

Pp2 = numpy.loadtxt('../MORAWE/ID18_DMM_Flux_PdC_E_0.2deg.txt', skiprows=1)
Pp4 = numpy.loadtxt('../MORAWE/ID18_DMM_Flux_PdC_E_0.4deg.txt', skiprows=1)

Wp2 = numpy.loadtxt('../MORAWE/ID18_DMM_Flux_WB4C_E_0.2deg.txt', skiprows=1)
Wp4 = numpy.loadtxt('../MORAWE/ID18_DMM_Flux_WB4C_E_0.4deg.txt', skiprows=1)


xTp2, yTp2 = calc_T(theta1=0.2)
xTp4, yTp4 = calc_T(theta1=0.4)

xPp2, yPp2 = calc_P(theta1=0.2)
xPp4, yPp4 = calc_P(theta1=0.4)

xWp2, yWp2 = calc_W(theta1=0.2)
xWp4, yWp4 = calc_W(theta1=0.4)


# [Ti(5.4)/C(11.9)]10  E/O
print("===================== [Ti(5.4)/C(11.9)]10  O/E ====================")
thickness = 5.4 + 11.9
gamma = 11.9 / thickness # E/(E+O)
print("Thickness, gamma: ", thickness, gamma)
lam_p2 = 2 * thickness * 1e-9 * numpy.sin(numpy.radians(0.2))
print("lambda [A}", lam_p2 * 1e10)
print("0.2 deg energy [keV]", 1e-3 * codata.h * codata.c / codata.e / lam_p2)
lam_p4 = 2 * thickness * 1e-9 * numpy.sin(numpy.radians(0.4))
print("lambda [A}", lam_p4 * 1e10)
print("0.4 deg energy [keV]", 1e-3 * codata.h * codata.c / codata.e / lam_p4)

print("===================== [Pd(3.9)/C(5.0)]15  O/E  ====================")
thickness = 3.9 + 5.0
gamma = 5.0 / thickness # E/(E+O)
print("Thickness, gamma: ", thickness, gamma)
lam_p2 = 2 * thickness * 1e-9 * numpy.sin(numpy.radians(0.2))
print("lambda [A}", lam_p2 * 1e10)
print("0.2 deg energy [keV]", 1e-3 * codata.h * codata.c / codata.e / lam_p2)
lam_p4 = 2 * thickness * 1e-9 * numpy.sin(numpy.radians(0.4))
print("lambda [A}", lam_p4 * 1e10)
print("0.4 deg energy [keV]", 1e-3 * codata.h * codata.c / codata.e / lam_p4)

print("===================== [W(2.0)/B4C(2.8)]30  O/E  ====================")
thickness = 2.0 + 2.8
gamma = 2.8 / thickness # E/(E+O)
print("Thickness, gamma: ", thickness, gamma)
lam_p2 = 2 * thickness * 1e-9 * numpy.sin(numpy.radians(0.2))
print("lambda [A}", lam_p2 * 1e10)
print("0.2 deg energy [keV]", 1e-3 * codata.h * codata.c / codata.e / lam_p2)
lam_p4 = 2 * thickness * 1e-9 * numpy.sin(numpy.radians(0.4))
print("lambda [A}", lam_p4 * 1e10)
print("0.4 deg energy [keV]", 1e-3 * codata.h * codata.c / codata.e / lam_p4)


plot(Tp2[:,0], Tp2[:, 1],
     xTp2 * 1e-3, yTp2**2 ,
     Timd[:,0], Timd[:,2]**2 ,
     # Tp4[:,0], Tp4[:, 1],
     # xTp4 * 1e-3, yTp4,
     title="[Ti(5.4)/C(11.9)]10 ",
     xtitle="E [keV]", ytitle="Reflectivity", xrange=[1,100], ylog=0,
     legend=["Morawe 0.2 deg", "Oasys 0.2 deg", "IMD 2.0 deg", "0.4 deg",  "Oasys 0.4 deg"], show=0,
     )

plot(Pp2[:,0], Pp2[:, 1],
     xPp2 * 1e-3, yPp2**2,
     # Tp4[:,0], Tp4[:, 1],
     # xPp4 * 1e-3, yPp4,
     title="[Pd(3.9)/C(5.0)]15",
     xtitle="E [keV]", ytitle="Reflectivity", xrange=[1,100], ylog=0,
     legend=["0.2 deg", "Oasys 0.2 deg", "0.4 deg",  "Oasys 0.4 deg"], show=0,
     )

plot(Wp2[:,0], Wp2[:, 1],
     xWp2 * 1e-3, yWp2**2,
     # Tp4[:,0], Tp4[:, 1],
     # xWp4 * 1e-3, yWp4,
     title="[W(2.0)/B4C(2.8)]30",
     xtitle="E [keV]", ytitle="Reflectivity", xrange=[1,100], ylog=0,
     legend=["0.2 deg", "Oasys 0.2 deg", "0.4 deg",  "Oasys 0.4 deg"], show=0,
     )


plot_show()

