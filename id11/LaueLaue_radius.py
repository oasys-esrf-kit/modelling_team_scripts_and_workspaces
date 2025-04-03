import numpy
import scipy.constants as codata


def wavelength_m(energy_keV):
    return codata.c * codata.h / codata.e / (1e3 * energy_keV)  # Wavelength in m

def rowland_radius(energy_keV=1.0, asymmetry_deg=0.0, D_m=10.0):
    """
    Calculate the Rowland radius (R_c) for a bent Si(111) Laue crystal.

    Parameters:
        energy_keV (float): X-ray energy in keV
        asymmetry_deg (float): Asymmetry angle in degrees
        D (float): Detector distance in m

    Returns:
        float: Rowland radius in m
    """
    # Constants
    lattice_constant_A = 5.431  # Silicon lattice constant in Angstroms
    d_111_A = lattice_constant_A / numpy.sqrt(3)  # Si(111) interplanar spacing in Angstroms

    # Convert energy to wavelength
    # wavelength_A = 12.398 / energy_keV  # Wavelength in Angstroms
    # print(codata.c * codata.h / codata.e)
    wavelength_A = 1e10 * wavelength_m(energy_keV) # codata.c * codata.h / codata.e / (1e3 * energy_keV)  # Wavelength in Angstroms

    # Bragg angle calculation
    theta_B_rad = numpy.arcsin(wavelength_A / (2 * d_111_A))  # in radians

    # Apply asymmetry angle correction
    asymmetry_rad = numpy.radians(asymmetry_deg)
    theta_prime_rad = theta_B_rad + asymmetry_rad  # Adjusted angle

    # Rowland radius formula
    Rc_m = D_m / (2 * numpy.sin(theta_prime_rad))
    return Rc_m


if __name__ == "__main__":
    lattice_constant_A = 5.431  # Silicon lattice constant in Angstroms
    d_111_A = lattice_constant_A / numpy.sqrt(3)  # Si(111) interplanar spacing in Angstroms


    F1 = 52.0
    F2 = 8.0

    for Energy_keV in [40.0, 60.0, 80.0, 100, 120, 140]:
        theta_B_rad = numpy.arcsin(1e10 * wavelength_m(Energy_keV) / (2 * d_111_A))  # in radians

        for chi_deg in [-20, -10, 0, 10, 20]:
            # Eq. 1 in https://doi.org/10.1016/S0168-9002(99)00056-X
            chi = numpy.radians(chi_deg)
            Ru = 2 * (numpy.cos(chi + theta_B_rad) / F1 + numpy.cos(chi - theta_B_rad) / F2)**(-1)
            Rl = 2 * (numpy.cos(chi - theta_B_rad) / F1 + numpy.cos(chi + theta_B_rad) / F2)**(-1)
            print("E=%d keV, Bragg angle: %.3f deg, Asymmetry angle = %.1f deg, Ru = %.3f m, , Rl = %.3f m" % \
                  (Energy_keV, numpy.degrees(theta_B_rad), chi_deg, Ru, Rl))
            # print(rowland_radius(energy_keV=Energy_keV, asymmetry_deg=chi_deg, D_m=F2) )



    # Not sure about this...
    energy_range = numpy.linspace(50, 100, 100)  # Energy range from 10 to 100 keV
    D_detector = F2  # m (meter)

    # Compute Rowland radii for energy range
    asymmetry_deg = 10
    Rc_values_10 = [rowland_radius(energy_keV=E, asymmetry_deg=asymmetry_deg, D_m=D_detector) for E in energy_range]

    from srxraylib.plot.gol import plot
    f, a = plot(energy_range, Rc_values_10, legend=f'Asymmetry = {asymmetry_deg}°', color='b',
                xtitle='Energy (keV)', ytitle='Rowland Radius (m)',
                title='Rowland Radius vs. Energy for Si(111) Laue Crystal', grid=1,
                )
