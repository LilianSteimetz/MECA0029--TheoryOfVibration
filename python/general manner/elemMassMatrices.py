import numpy as np
from constants import *
from math import pi


def mass_matrix(L, A, I, rho):
    r_squared = I / A

    M = np.zeros((12, 12))

    M[0, 0] = 1/3
    M[0, 6] = 1/6

    M[1, 1] = 13/35
    M[1, 5] = 11*L/210
    M[1, 7] = 9/70
    M[1, 11] = -13*L/420

    M[2, 2] = 13/35
    M[2, 4] = -11*L/210
    M[2, 8] = 9/70
    M[2, 10] = 13*L/420

    M[3, 3] = r_squared/3
    M[3, 9] = r_squared/6

    M[4, 4] = L**2/105
    M[4, 8] = -13*L/420
    M[4, 10] = -L**2/140

    M[5, 5] = L**2/105
    M[5, 7] = 13*L/420
    M[5, 11] = -L**2/140

    M[6, 6] = 1/3

    M[7, 7] = 13/35
    M[7, 11] = -11*L/210

    M[8, 8] = 13/35
    M[8, 10] = 11*L/210

    M[9, 9] = r_squared/3

    M[10, 10] = L**2/105

    M[11, 11] = L**2/105

    # Symmetrize
    M = M + np.triu(M, 1).T

    return rho * A * L * M
