import numpy


def pos(theta, theta0=0.0, d=1.0):
    return d / numpy.sin(theta) * numpy.cos(theta) - d / numpy.sin(theta0) * numpy.cos(theta0)

def angle_from_pos(pos, theta0=0.0, d=1.0):
    return numpy.arctan(1.0 / (pos / d + 1.0 / numpy.tan(theta0)))

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



