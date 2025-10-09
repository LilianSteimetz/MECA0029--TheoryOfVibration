from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure
import numpy as np
from math import pi

# compute the mass using a rbm translation

""" Computation of the mass using the rbm in translation """


def compute_total_mass(M):

    nDof = M.shape[0]
    nNodes = nDof // 6

    translat_dofs = []
    for i in range(nNodes):
        translat_dofs.extend([i*6+1])  # ux, uy, uz

    # Extract translational block of the global mass matrix
    M_translat = M[np.ix_(translat_dofs, translat_dofs)]

    # Compute total translational mass
    total_mass = np.sum(M_translat)

    return total_mass


M, K = create_globalMass_and_globalStiffness(constrainedNodes=[])
total_mass = compute_total_mass(M)
print(
    f"Total mass of the structure, using translation rbm : {total_mass:.2f} kg")

""" Computation of the real mass"""

A1 = pi * (ro1**2 - (ro1 - e1)**2)
A2 = pi * (ro2**2 - (ro2 - e2)**2)

L1 = 3.0
L2 = math.sqrt(1.5**2 + 1**2)
L3 = 4.0

m_type1 = A1 * L1 * rho
m_type2 = A1 * L2 * rho
m_type3 = A2 * L2 * rho
m_type4 = A2 * L3 * rho
n_type1 = 18
n_type2 = 4
n_type3 = 16
n_type4 = 9

real_mass = n_type1 * m_type1 + n_type2 * m_type2 + \
    n_type3 * m_type3 + n_type4 * m_type4 + 500
print(
    f"Real mass of the structure : {total_mass:.2f} kg")
