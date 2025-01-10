import numpy

energies = [7, 17, 35, 40]

for energy in energies:
    """
    =============== photon energy:  7000 ==========================
    Values are FWHM
    Hsize 10.805 67.690 1.073 68.555
    Vsize 10.805 12.141 1.073 16.288
    Hangle 13.677 9.844 6.499 18.062
    Vangle 13.677 4.568 6.499 15.817
    
    =============== photon energy:  17000 ==========================
    Values are FWHM
    Hsize 6.934 67.690 2.021 68.074
    Vsize 6.934 12.141 2.021 14.126
    Hangle 8.777 9.844 10.387 16.788
    Vangle 8.777 4.568 10.387 14.345
    
    =============== photon energy:  35000 ==========================
    Values are FWHM
    Hsize 4.832 67.690 2.248 67.899
    Vsize 4.832 12.141 2.248 13.259
    Hangle 6.117 9.844 10.093 15.369
    Vangle 6.117 4.568 10.093 12.655
    
    =============== photon energy:  40000 ==========================
    Values are FWHM
    Hsize 4.520 67.690 2.103 67.873
    Vsize 4.520 12.141 2.103 13.124
    Hangle 5.722 9.844 9.441 14.792
    Vangle 5.722 4.568 9.441 11.948
    """

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


    DISTANCES = [28,30,40,45,50,200-0.1, 200-0.05, 200]

    print("\n\npropagated beam size and divergences at E=%d keV" % (energy))
    print("distance[m]  Hsize[um]  Vsize[um]")

    for distance in DISTANCES:
        # distance = 40.0


        print("%.2f %.0f %.0f" % (
            distance,
            numpy.sqrt(s_H**2 + (a_H * distance)**2),
            numpy.sqrt(s_V**2 + (a_V * distance)**2),
            # a_H,
            # a_V
            ))

