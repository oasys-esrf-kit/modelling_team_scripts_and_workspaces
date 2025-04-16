# calculation of the multilayer shape (see https://confluence.esrf.fr/display/id18nanotomo/Multilayer+Monochromator )

import numpy
from srxraylib.plot.gol import plot, plot_show

def h(x, Alpha=1.0, c1=0, c2=0):
    return Alpha * x * (numpy.log(x) - 1) + c1 * x + c2

if __name__ == "__main__":

    offset = 0.006
    offset_crystals = offset / 2
    SecS = 54.4
    p = 29.5
    q = SecS - p # 24.9
    M1_length = 0.300
    print('p, q, offset_crystals: ', p, q, offset_crystals)

    f = 1.0 / (1 / p + 1 / q)
    print("f: ", f )

    theta_min = 3.5e-3
    theta_max = 7e-3

    print("angle limits: ", theta_min, theta_max, "in degrees: ", numpy.degrees(theta_min), numpy.degrees(theta_max))

    x_min = offset_crystals / numpy.sin(theta_max) * numpy.cos(theta_max)
    x_max = offset_crystals / numpy.sin(theta_min) * numpy.cos(theta_min)
    x_0 = 0.5 * (x_min + x_max)
    M2_length = (x_max - x_min) + M1_length
    print("x limits: ", x_min, x_max, x_0, M2_length)

    theta_0 = numpy.arctan(offset_crystals / x_0)
    print("theta_0 = ", theta_0)

    npoints = 100
    Alpha = offset / 4 / f
    x = numpy.linspace(x_0 - 0.5 * M2_length, x_0 + 0.5 * M2_length, npoints)

    #
    # get c1, c2
    #
    print(h(x_min, Alpha=Alpha))
    print(h(x_max, Alpha=Alpha))
    slope = (h(x_max, Alpha=Alpha) - h(x_min, Alpha=Alpha)) / (x_max - x_min)
    c1 = -numpy.tan(slope)

    c2 = -numpy.min(h(x, Alpha=Alpha, c1=c1))
    print('Alpha: ', Alpha)
    print('c1, c2: ', c1, c2)


    #
    # shadow parameters:
    #
    slit = 27.3
    p1 = (p - slit)
    D = x_0 / numpy.cos(theta_0)
    q1 = 0.5 * D
    p2 = q1
    q2 = SecS - (slit + p1 + q1 + p2)
    print("==============SHADOW parameters: ==================")
    print("entrance slit at: ", slit)
    print("p1, q1: ", p1, q1)
    print("D (inter-mirror distance along the beam): ", D, x_0)
    print("p2, q2: ", p2, q2)
    print("theta_0 [mrad]", 1e3 * theta_0)
    print("Radius: ", 2 * f / numpy.sin(theta_0))
    print("Demagnification: ", q / p)
    print("M2 focal distances: ", slit + p1 + q1 + p2, q2)
    print("Magnification: ", q / p, q2 / (slit + p1 + q1 + p2))
    print("Size at SecS: ", 64.8 * q / p, 64.8 * q2 / (slit + p1 + q1 + p2))
    print("Limits in y: ", x.min() - x_0, x.max() - x_0)


    yy = (x - x_0)
    zz = h(x, Alpha=Alpha, c1=c1, c2=c2)

    #
    # write oasys surface file
    #

    npoints_x = 21
    xx = numpy.linspace(-0.005, 0.005, npoints_x)
    Z = numpy.outer(numpy.ones_like(xx), zz)

    from shadow4.optical_surfaces.s4_mesh import S4Mesh
    m = S4Mesh(surface = None,
                 mesh_x=xx,
                 mesh_y=yy,
                 mesh_z=Z )
    m.write_mesh_h5file(xx, yy, filename="M2_shape.h5")

    #
    # plot final
    #
    plot(yy, 1e6 * zz, grid=1, xtitle="x [m]", ytitle="Height [um]", show=0)
    plot(yy, 1e6 * numpy.gradient(zz, yy), grid=1, xtitle="x [m]", ytitle="Slope [urad]", show=0)
    Curv =  numpy.gradient(numpy.gradient(zz, yy), yy)
    plot(yy, Curv, grid=1, xtitle="x [m]", ytitle="Curvature [m^-1]", show=0)
    plot(yy, 1.0 / Curv, grid=1, xtitle="x [m]", ytitle="Radius [m]", show=0)

    # plot(x - x_0, h(x, Alpha=Alpha),
    #      x - x_0, h(x, Alpha=Alpha, c1=c1),
    #      x - x_0, h(x, Alpha=Alpha, c1=c1, c2=c2),
    #      grid=1, show=0)

    plot_show()



