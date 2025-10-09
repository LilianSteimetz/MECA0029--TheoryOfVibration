import numpy as np
from constants import *
from math import pi

G = E / (2 * (1 + nu))

ri1 = ro1 - e1
ri2 = ro2 - e2

A1 = pi * (ro1**2 - ri1**2)
Iy1 = (pi / 4.0) * (ro1**4 - ri1**4)
Iz1 = Iy1
Jx1 = (pi / 2.0) * (ro1**4 - ri1**4)

A2 = pi * (ro2**2 - ri2**2)
Iy2 = (pi / 4.0) * (ro2**4 - ri2**4)
Iz2 = Iy2
Jx2 = (pi / 2.0) * (ro2**4 - ri2**4)


def stiffness_matrix(L, A, Iy, Iz, Jx, E, G):

    K = np.zeros((12, 12))

    K[0, 0] = E*A/L
    K[0, 6] = -E*A/L

    K[1, 1] = 12*E*Iz / L**3
    K[1, 5] = 6*E*Iz / L**2
    K[1, 7] = -12*E*Iz / L**3
    K[1, 11] = 6*E*Iz / L**2

    K[2, 2] = 12*E*Iy / L**3
    K[2, 4] = -6*E*Iy / L**2
    K[2, 8] = -12*E*Iy / L**3
    K[2, 10] = -6*E*Iy / L**2

    K[3, 3] = G*Jx / L
    K[3, 9] = -G*Jx / L

    K[4, 4] = 4*E*Iy / L
    K[4, 8] = 6*E*Iy / L**2
    K[4, 10] = 2*E*Iy / L

    K[5, 5] = 4*E*Iz / L
    K[5, 7] = -6*E*Iz / L**2
    K[5, 11] = 2*E*Iz / L

    K[6, 6] = E*A / L
    K[7, 7] = 12*E*Iz / L**3
    K[7, 11] = -6*E*Iz / L**2
    K[8, 8] = 12*E*Iy / L**3
    K[8, 10] = 6*E*Iy / L**2
    K[9, 9] = G*Jx / L
    K[10, 10] = 4*E*Iy / L
    K[11, 11] = 4*E*Iz / L

    K = K + np.triu(K, 1).T

    return K
