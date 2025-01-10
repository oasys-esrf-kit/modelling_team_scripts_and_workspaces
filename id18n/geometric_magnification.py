import numpy
from srxraylib.plot.gol import plot

def calculate(p1=28.0):
    #
    # inputs
    #


    # p1 = 29.0 # horizontal
    q2 = 0.05 # horizontal
    q = 0.1 # vertical
    D = 200.0

    #
    #
    #
    q1 = numpy.linspace(40, 55, 100) - p1
    p2 = D - (p1 + q1 + q2)
    p = D - q # this is for vertical

    return p1, q1, p2, q2, p, q



if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # p1, q1, p2, q2, p, q = calculate(p1=28.0)
    #
    # s_H = 73.7
    # s_V = 14.2
    # a_H = 17.3
    # a_V = 14.3
    #
    #
    # plot(q1 + p1, 1 / (q1 * q2 / p1 / p2),
    #      q1 + p1, numpy.ones_like(q1) * 1 / (q / p),
    #      xtitle="Position of secondary H source [m]",
    #      legend=["H","V"], ytitle="1/|M|", title="DEMAGNIFICATION", show=0)
    #
    # plot(q1 + p1, 1e3 * s_H * (q1 * q2 / p1 / p2),
    #      q1 + p1, numpy.ones_like(q1) * 1e3 * s_V * (q / p),
    #      xtitle="Position of secondary H source [m]", ytitle="sample size [nm]", title="SIZE", legend=["H","V"],  show=0)
    #
    # plot(q1 + p1, a_H / (q1 / p1),
    #      q1 + p1, numpy.ones_like(q1) * a_V,
    #      legend=["H","V"], xtitle="Position of secondary H source [m]", ytitle="divergence at KB [urad]",
    #      title="DIVERGENCE AT MONO OR KB", show=1)


    SCAN = [28.0,29.0,30.0]
    energy = 40

    P1 = []
    Q1 = []
    P2 = []
    Q2 = []
    P = []
    Q = []

    for i in range(len(SCAN)):
        p1, q1, p2, q2, p, q = calculate(p1=SCAN[i])
        P1.append(p1)
        Q1.append(q1)
        P2.append(p2)
        Q2.append(q2)
        P.append(p)
        Q.append(q)

    if energy == 7:
        s_H = 68.555
        s_V = 16.288
        a_H = 18.062
        a_V = 15.817
        darwin = 28.5
    elif energy == 17:
        s_H = 68.074
        s_V = 14.126
        a_H =  16.788
        a_V =  14.345
        darwin = 14.9
    elif energy == 35:
        s_H = 67.899
        s_V = 13.259
        a_H =  15.369
        a_V =  12.655
        darwin = 7.3
    elif energy == 40:
        s_H = 67.873
        s_V = 13.124
        a_H = 14.792
        a_V = 11.948
        darwin = 6.4
    else:
        raise Exception("bad energy")



    plot(Q1[0] + P1[0], 1 / (Q1[0] * Q2[0] / P1[0] / P2[0]),
         Q1[1] + P1[1], 1 / (Q1[1] * Q2[1] / P1[1] / P2[1]),
         Q1[2] + P1[2], 1 / (Q1[2] * Q2[2] / P1[2] / P2[2]),
         Q1[0] + P1[0], numpy.ones_like(Q1[0]) * 1 / (Q[0] / P[0]),
         xtitle="Position of secondary H source [m]",
         legend=["H HF@28m","H HF@29m","H HF@30m","V"], ytitle="1/|M|", title="DEMAGNIFICATION E=%d keV" % energy, show=0)
    plt.grid()

    plot(Q1[0] + P1[0], 1e3 * s_H * (Q1[0] * Q2[0] / P1[0] / P2[0]),
         Q1[1] + P1[1], 1e3 * s_H * (Q1[1] * Q2[1] / P1[1] / P2[1]),
         Q1[2] + P1[2], 1e3 * s_H * (Q1[2] * Q2[2] / P1[2] / P2[2]),
         Q1[0] + P1[0], numpy.ones_like(Q1[0]) * 1e3 * s_V * (Q[0] / P[0]),
         xtitle="Position of secondary H source [m]", ytitle="sample size [nm]", title="SIZE E=%d keV" % energy,
         legend=["H HF@28m","H HF@29m","H HF@30m","V"],  show=0)
    plt.grid()

    plot(Q1[0] + P1[0], a_H / (Q1[0] / P1[0]),
         Q1[1] + P1[1], a_H / (Q1[1] / P1[1]),
         Q1[2] + P1[2], a_H / (Q1[2] / P1[2]),
         Q1[1] + P1[1], numpy.ones_like(Q1[0]) * a_V,
         Q1[1] + P1[1], numpy.ones_like(Q1[0]) * darwin,
         legend=["H HF@28m","H HF@29m","H HF@30m","V","DARWIN %d keV" % energy],
         xtitle="Position of secondary H source [m]", ytitle="divergence at KB [urad]",

         title="DIVERGENCE AT MONO OR KB E=%d keV" % energy, show=0)
    plt.grid()
    plt.show()