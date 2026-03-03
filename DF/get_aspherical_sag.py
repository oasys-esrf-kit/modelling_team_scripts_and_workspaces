import numpy as np
import matplotlib.pylab as plt


def get_aspheric_sag(r):
    """
    Calculates the sag Z(r) for Edmund Optics #49-663 in METERS
    using values directly from Technical Drawing #376506.
    """
    # Parameters from Drawing (converted to Meters)
    R = 0.32301  # Radius of Curvature (323.01 mm)
    c = 1 / R  # Curvature
    k = 0.0  # Conic Constant

    # Coefficients from Drawing (converted from mm to m)
    # Conversion: F_meters = F_mm * (1000)^3
    F = -2.0800e-8 * (1000 ** 3)  # Result: -20.8
    G = -1.1100e-11 * (1000 ** 5)  # Result: -1.11e3
    H = 0.0

    # Even Asphere Formula
    z_base = (c * r ** 2) / (1 + np.sqrt(1 - (c ** 2 * r ** 2)))
    z_asphere = (F * r ** 4) + (G * r ** 6) + (H * r ** 8)

    return z_base + z_asphere


# Quick verification:
# At r = 12.5mm (edge of the 25mm lens)
print(f"Sag at edge (r=0.0125m): {get_aspheric_sag(0.0125) * 1000:.6f} mm")


# --- Visualization / Test Block ---
radius_max = 0.0125  # 12.5mm semi-aperture for a 25mm lens
r_values = np.linspace(0, radius_max, 100)
z_values = get_aspheric_sag(r_values)

plt.figure(figsize=(8, 5))
plt.plot(r_values * 1000, z_values * 1000, label="S4 Profile (#49-663)")
plt.xlabel("Radial Distance r (mm)")
plt.ylabel("Sag Z (mm)")
plt.title("Back Surface Profile (S4) of Aspherized Achromat")
plt.gca().invert_yaxis()  # Optical sag is usually plotted downward
plt.grid(True, linestyle='--')
plt.legend()
plt.show()

# Print a specific value for verification
# At r = 10mm (0.01m)
test_r = 0.01
print(f"Sag at r=10mm: {get_aspheric_sag(test_r) * 1000:.6f} mm")