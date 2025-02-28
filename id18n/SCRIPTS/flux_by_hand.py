import numpy

#KB mirror lengths
Lv = 60e-3
Lh = 25.07e-3


#
# 7 keV
#

# source sizes and angles
Sh, Sv, Ah, Av = 68.55e-6, 16.30e-6, 18.06e-6, 15.80e-6
theta = 0.0364

# demag V
p_V = 199.900
q_V = 0.1
T_V = Lv * numpy.sin(theta) / (Av * p_V)
print("Transmission in V: ", T_V, 150000 * T_V)


# demag config 2
DM1_config2 = 3.02
DM2_config2 = 3107
p2_config2 = 155.35
# KB angular acceptance = L sin(theta) / p2
A_KB = Lh * numpy.sin(theta) / p2_config2
# Transmission config 2
T_config2 = A_KB / numpy.sqrt( (Ah * DM1_config2)**2 + ((Sh / DM1_config2) / p2_config2)**2 )
print("Config2 H: ", T_config2, 150000 * T_config2)


# demag config 3
DM1_config3 = 1.4
DM2_config3 = 2919
p2_config3 = 145.95
# KB angular acceptance = L sin(theta) / p2
A_KB = Lh * numpy.sin(theta) / p2_config3
# Transmission config 2
T_config3 = A_KB / numpy.sqrt( (Ah * DM1_config3)**2 + ((Sh / DM1_config3) / p2_config3)**2 )
T_slit = 29.19e-6 / (Sh / DM1_config3)
print("Config3 H Slit: ", T_slit, T_slit * 150000)
print("Config3 H KB: ", T_config3, T_config3 * 150000)
T_config3 *= T_slit
print("Config3 H: ", T_config3, 150000 * T_config3)


print("config 2 final", T_V * T_config2, int(150000 * T_V * T_config2))
print("config 3 final", T_V * T_config3, int(150000 * T_V * T_config3))

print("\n\n\n")
print("divergence after HSS config 2: ", Ah * DM1_config2)
print("divergence after HSS config 3: ", Ah * DM1_config3)