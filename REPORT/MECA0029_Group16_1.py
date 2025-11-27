from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure
import numpy as np
# eigsh is used like eig but computes a given number of eigenvalues
from scipy.sparse.linalg import eigsh
import pandas as pd
from geometry import constrainedNodes

###########################################################################
#                                PART1                                    #
#             Solve The natural frequencies and mode shapes               #
#                                                                         #
###########################################################################


""" Computation of the 6 1st natural frequencies and mode shapes of the structure """
M_global, K_global = create_globalMass_and_globalStiffness()

eigvals, eigvecs = eigsh(K_global, k=desiredFreqNb,
                         M=M_global, sigma=0.0, which='LM')

nat_freqs = np.sqrt(eigvals)/(2*np.pi)

"""Sort natural frequencies and mode shapes"""
idx = np.argsort(nat_freqs)
nat_freqs = nat_freqs[idx]
eigvecs = eigvecs[:, idx]
eigvecs = eigvecs/np.max(np.abs(eigvecs), axis=0)
print("Natural Frequencies (Hz):")
print(nat_freqs)


"""Add the constrained DOFs back (0 displacement)"""
constrainedDOFs = []
for i in range(len(constrainedNodes)):
    node = constrainedNodes[i]
    constrainedDOFs.extend(dofList[node-1, :] - 1)
constrainedDOFs = constrainedDOFs[::-1]

nDOF = np.max(dofList)   # total number of DOFs
nNodes = nodeList.shape[0]
nModes = eigvecs.shape[1]

eigvecs_full = np.zeros((nDOF, nModes))
allDOFs = np.arange(nDOF)
freeDOFs = np.setdiff1d(allDOFs, constrainedDOFs)
eigvecs_full[freeDOFs, :] = eigvecs
eigvecs = eigvecs_full


"""Visualization of a mode shape"""
mode_idx = 3  # mode to visualize


U = np.zeros((nNodes, 3))  # x, y, z displacement per node
for i in range(nNodes):
    U[i, 0] = eigvecs[dofList[i, 0]-1, mode_idx]
    U[i, 1] = eigvecs[dofList[i, 1]-1, mode_idx]
    U[i, 2] = eigvecs[dofList[i, 2]-1, mode_idx]


plot_structure(elemList, nodeList + U*1)


###########################################################################
#                                PART1                                    #
#             computing the total mass with rigid body modes              #
#                                                                         #
###########################################################################


def compute_total_mass(M):

    nDof = M.shape[0]
    nNodes = nDof // 6

    translat_dofs = []
    for i in range(nNodes):
        translat_dofs.extend([i*6+0])  # x
        translat_dofs.extend([i*6+1])  # y
        translat_dofs.extend([i*6+2])  # z

    # Extract translational block of the global mass matrix and sum it
    # is the same thing as doing u^T M u
    M_translat = M[np.ix_(translat_dofs, translat_dofs)]
    total_mass = np.sum(M_translat)

    # because mass is counted 3 times
    return total_mass / 3


M, K = create_globalMass_and_globalStiffness(constrainedNodes=[])
total_mass = compute_total_mass(M)

print(
    f"Total mass of the structure, using translation rbm : {total_mass:.6f} kg")


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
    f"Real mass of the structure : {real_mass:.2f} kg")
