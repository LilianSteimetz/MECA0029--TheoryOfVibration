import math


elemPerBar = 5


l_horizontal = 3.0 / elemPerBar
l_diagonal = math.sqrt(1.5**2 + 1**2) / elemPerBar
l_transverse = 4 / elemPerBar

rho = 7800  # Density in kg/m^3
nu = 0.3
E = 210 * 10**9  # Young's modulus in Pa
ro1 = 0.06  # Outer radius in m
e1 = 0.005  # Thickness in m
ro2 = 0.035  # Outer radius in m
e2 = 0.003  # Thickness in m

lumpedMass = 500/8
