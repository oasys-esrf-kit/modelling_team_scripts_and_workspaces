import numpy
from srxraylib.plot.gol import plot, plot_show


def pos(theta, theta0=0.0, d=1.0):
    return d / numpy.sin(theta) * numpy.cos(theta) - d / numpy.sin(theta0) * numpy.cos(theta0)

def angle_from_pos(pos, theta0=0.0, d=1.0):
    return numpy.arctan(1.0 / (pos / d + 1.0 / numpy.tan(theta0)))


def calc1():
    theta0 = numpy.radians(0.3)
    theta1 = numpy.radians(0.2)
    theta2 = numpy.radians(0.4)

    d = 0.004

    M2 = 31.5

    d_M1M2 = d / numpy.sin(theta0)
    M1 = 31.5 - d_M1M2

    pos0 = pos(theta0, theta0=theta0, d=d)  # d / numpy.sin(theta0) * numpy.cos(theta0) -  d_M1M2 * numpy.cos(theta0)
    pos1 = pos(theta1, theta0=theta0, d=d)  # d / numpy.sin(theta1) * numpy.cos(theta1) -  d_M1M2 * numpy.cos(theta0)
    pos2 = pos(theta2, theta0=theta0, d=d)  # d / numpy.sin(theta2) * numpy.cos(theta2) -  d_M1M2 * numpy.cos(theta0)

    print("d [m] = ", d)
    print("M1 [m] = ", M1)
    print("M2 [m] = ", M2)
    print("d_M1M2 [m] = ", d_M1M2)
    print("pos0 [m] = ", pos0)
    print("pos1 [m] = ", pos1)
    print("pos2 [m] = ", pos2)

    q0 = 54.0 - 31.5
    q1 = q0 - pos1 / numpy.cos(theta1)
    q2 = q0 - pos2 / numpy.cos(theta2)

    R0 = 2 / ((1/M2 + 1/q0) * numpy.sin(theta0))
    R1 = 2 / ((1/M2 + 1/q1) * numpy.sin(theta1))
    R2 = 2 / ((1/M2 + 1/q2) * numpy.sin(theta2))


    print("q0 [m] = ", q0)
    print("q1 [m] = ", q1)
    print("q2 [m] = ", q2)


    print("R0 [m] = ", R0)
    print("R1 [m] = ", R1)
    print("R2 [m] = ", R2)

    th0 = angle_from_pos(pos0, theta0=theta0, d=d)  # numpy.arctan( 1.0 / (pos0 / d + 1.0 / numpy.tan(theta0) ))
    th1 = angle_from_pos(pos1, theta0=theta0, d=d)  # numpy.arctan( 1.0 / (pos1 / d + 1.0 / numpy.tan(theta0) ))
    th2 = angle_from_pos(pos2, theta0=theta0, d=d)  # numpy.arctan( 1.0 / (pos2 / d + 1.0 / numpy.tan(theta0) ))

    print(theta0, th0)
    print(theta1, th1)
    print(theta2, th2)

    npoints = 200
    xpos = numpy.linspace(pos2, pos1, npoints)
    thetapos = angle_from_pos(xpos, theta0=theta0, d=d)
    qpos = q0 - xpos / numpy.cos(thetapos)
    Rpos = 2 / ((1/M2 + 1/qpos) * numpy.sin(thetapos))

    from srxraylib.plot.gol import plot
    plot(xpos, Rpos, xtitle="position on M2 [m]", ytitle="R[m]", show=0)

    coeff = numpy.polynomial.polynomial.polyfit(xpos, Rpos, 1)
    print(coeff)

    plot(xpos, Rpos,
         xpos, xpos * coeff[1] + coeff[0], xtitle="position on M2 [m]", ytitle="R[m]")


    plot(xpos, numpy.cumsum(1/Rpos), xtitle="position on M2 [m]", ytitle="height [m]", show=0)


    R00 = coeff[0]
    alpha00 = coeff[1]
    ss = R00 + alpha00 * xpos
    yy = (ss * numpy.log(numpy.abs(ss)) - ss)/ alpha00**2

    plot(xpos, yy, xtitle="position on M2 [m]", ytitle="fit height [m]")

    plot(xpos, 1.0 / numpy.gradient(yy, xpos), xtitle="position on M2 [m]", ytitle="R from slope []")

def integral_result(x, R, c):
    if c == 0:
        return x**2 / (2 * R)  # Special case for c = 0
    else:
        return ((R + c*x)**2) / (2 * R * c) - (R / c**2) * numpy.log(numpy.abs(R + c*x))

if __name__ == "__main__":
    # calc1()


    d = 0.003
    M1 = 29.000
    M2 = 29.700
    L_M2 = 0.8

    thetaN = numpy.arctan(d / (M2 - M1))
    M1M2 = numpy.sqrt(d**2 + (M2 - M1)**2)

    # pos0 = pos(theta0, theta0=theta0,
    #            d=d)  # d / numpy.sin(theta0) * numpy.cos(theta0) -  d_M1M2 * numpy.cos(theta0)
    # pos1 = pos(theta1, theta0=theta0,
    #            d=d)  # d / numpy.sin(theta1) * numpy.cos(theta1) -  d_M1M2 * numpy.cos(theta0)
    # pos2 = pos(theta2, theta0=theta0,
    #            d=d)  # d / numpy.sin(theta2) * numpy.cos(theta2) -  d_M1M2 * numpy.cos(theta0)

    print("d [m] = ", d)
    print("M1 [m] = ", M1)
    print("M2 [m] = ", M2)
    print("M1M2 [m] = ", M1M2)
    print("Nominal thetaN [deg] = ", numpy.degrees(thetaN))

    theta1 = numpy.radians(0.2)
    theta2 = numpy.radians(0.4)
    # theta11 = numpy.arctan( 1.0 / (-L_M2 / d + numpy.cos(thetaN) / numpy.sin(thetaN)) )
    # theta22 = numpy.arctan( 1.0 / ( L_M2 / d + numpy.cos(thetaN) / numpy.sin(thetaN)))
    # theta22 = numpy.arcsin(1.0 / (-L_M2 / 2 + 1.0 / numpy.sin(thetaN)))
    # print("theta1, theta11 [deg] = ", numpy.degrees(theta1), numpy.degrees(theta11),)
    # print("theta2, theta22 [deg] = ", numpy.degrees(theta2), numpy.degrees(theta22),)

    xN = d * numpy.cos(thetaN) / numpy.sin(thetaN)
    x1 = d * numpy.cos(theta1) / numpy.sin(theta1)
    x2 = d * numpy.cos(theta2) / numpy.sin(theta2)
    print("x2, xN, X1, x1-x2 = ", x2, xN, x1, x1-x2)
    print("x2-xN, xN-xN, x1-xN = ", x2-xN, xN-xN, x1-xN)

    xpos = numpy.linspace(x2-xN, x1-xN, 100)
    thetapos = angle_from_pos(xpos, theta0=thetaN, d=d)
    plot(thetapos, numpy.degrees(thetapos), xtitle="position on M2 [m]", ytitle="angle [deg]", show=0)

    q0 = 54.0 - M2
    qpos = q0 - xpos / numpy.cos(thetapos)
    Rpos = 2 / ((1 / M2 + 1 / qpos) * numpy.sin(thetapos))
    RN = 2 / ((1 / M2 + 1 / q0) * numpy.sin(thetaN))
    print("RN = ", RN)

    # plot(xpos, Rpos, xtitle="position on M2 [m]", ytitle="R[m]")

    # linear fot or R
    coeff = numpy.polynomial.polynomial.polyfit(xpos, Rpos, 1)
    print(coeff)
    R00 = coeff[0]
    alpha00 = coeff[1]
    ss = R00 + alpha00 * xpos
    plot(xpos, Rpos, xpos, ss, xtitle="position on M2 [m]", ytitle="R[m]", legend=['calculated', 'linear fit'], show=0)

    y0 = RN * (1 - numpy.sqrt(1 - (xpos / RN)**2 ))


    y0series = xpos**2 / (2 * RN)


    ii = numpy.argmin(numpy.abs(Rpos - R00))
    print('ii: ', ii)
    y = ((R00 + alpha00 * xpos) / alpha00) - (R00 / alpha00) * numpy.log(numpy.abs(R00 + alpha00 * xpos))
    y2 = integral_result(xpos, R00, alpha00)

    coeff2 = numpy.polynomial.polynomial.polyfit(xpos, y2, 1)
    y2_base = coeff2[0] + xpos * coeff2[1]
    y2 = y2 - y2_base
    y2 = y2 - y2[ii]



    plot(xpos, y0,
         xpos, y0series,
         xpos, y - y[ii],
         xpos, y2 - y2,
         legend=['fix exact', 'fix series', 'variable', 'variable 2'], ylog=1)


    # plot(xpos, yy, xtitle="position on M2 [m]", ytitle="fit height [m]")
    # plot(xpos, 1.0 / numpy.gradient(yy, xpos), xtitle="position on M2 [m]", ytitle="R from slope []")




    # print("pos0 [m] = ", pos0)
    # print("pos1 [m] = ", pos1)
    # print("pos2 [m] = ", pos2)
    #
    # q0 = 54.0 - 31.5
    # q1 = q0 - pos1 / numpy.cos(theta1)
    # q2 = q0 - pos2 / numpy.cos(theta2)
    #
    # R0 = 2 / ((1 / M2 + 1 / q0) * numpy.sin(theta0))
    # R1 = 2 / ((1 / M2 + 1 / q1) * numpy.sin(theta1))
    # R2 = 2 / ((1 / M2 + 1 / q2) * numpy.sin(theta2))
    #
    # print("q0 [m] = ", q0)
    # print("q1 [m] = ", q1)
    # print("q2 [m] = ", q2)
    #
    # print("R0 [m] = ", R0)
    # print("R1 [m] = ", R1)
    # print("R2 [m] = ", R2)
    #
    # th0 = angle_from_pos(pos0, theta0=theta0, d=d)  # numpy.arctan( 1.0 / (pos0 / d + 1.0 / numpy.tan(theta0) ))
    # th1 = angle_from_pos(pos1, theta0=theta0, d=d)  # numpy.arctan( 1.0 / (pos1 / d + 1.0 / numpy.tan(theta0) ))
    # th2 = angle_from_pos(pos2, theta0=theta0, d=d)  # numpy.arctan( 1.0 / (pos2 / d + 1.0 / numpy.tan(theta0) ))
    #
    # print(theta0, th0)
    # print(theta1, th1)
    # print(theta2, th2)
    #
    # npoints = 200
    # xpos = numpy.linspace(pos2, pos1, npoints)
    # thetapos = angle_from_pos(xpos, theta0=theta0, d=d)
    # qpos = q0 - xpos / numpy.cos(thetapos)
    # Rpos = 2 / ((1 / M2 + 1 / qpos) * numpy.sin(thetapos))
    #
    # from srxraylib.plot.gol import plot
    # plot(xpos, Rpos, xtitle="position on M2 [m]", ytitle="R[m]", show=0)
    #
    # coeff = numpy.polynomial.polynomial.polyfit(xpos, Rpos, 1)
    # print(coeff)
    #
    # plot(xpos, Rpos,
    #      xpos, xpos * coeff[1] + coeff[0], xtitle="position on M2 [m]", ytitle="R[m]")
    #
    # plot(xpos, numpy.cumsum(1 / Rpos), xtitle="position on M2 [m]", ytitle="height [m]", show=0)
    #
    # R00 = coeff[0]
    # alpha00 = coeff[1]
    # ss = R00 + alpha00 * xpos
    # yy = (ss * numpy.log(numpy.abs(ss)) - ss) / alpha00 ** 2
    #
    # plot(xpos, yy, xtitle="position on M2 [m]", ytitle="fit height [m]")
    #
    # plot(xpos, 1.0 / numpy.gradient(yy, xpos), xtitle="position on M2 [m]", ytitle="R from slope []")