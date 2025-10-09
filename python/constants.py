import math
from math import pi


elemPerBar = 5
lumpedMass = 500
desiredFreqNb = 6


rho = 7800  # Density in kg/m^3
nu = 0.3
E = 210 * 10**9  # Young's modulus in Pa
G = E / (2 * (1 + nu))  # Shear modulus in Pa


# Type 1 properties ( type = type of section, length of elem not used)
# A, Iy, Iz, Jx must be defined for each type of section,
# they are used inside the global stiffness and mass matrices
ro1 = 0.06  # Outer radius in m
e1 = 0.005  # Thickness in m
ri1 = ro1 - e1

A1 = pi * (ro1**2 - ri1**2)
Iy1 = (pi / 4.0) * (ro1**4 - ri1**4)
Iz1 = Iy1
Jx1 = (pi / 2.0) * (ro1**4 - ri1**4)

# Type 2 properties ( type = type of section, length of elem not used)
# A, Iy, Iz, Jx must be defined for each type of section,
# they are used inside the global stiffness and mass matrices

ro2 = 0.035  # Outer radius in m
e2 = 0.003  # Thickness in m
ri2 = ro2 - e2

A2 = pi * (ro2**2 - ri2**2)
Iy2 = (pi / 4.0) * (ro2**4 - ri2**4)
Iz2 = Iy2
Jx2 = (pi / 2.0) * (ro2**4 - ri2**4)


# IF more than 2 types of sections are needed, just add them here,
# then go in the file "globalMassStiffMatrices.py" and add an elif clause
# for that type, like what is currently done for type 1 and type 2
