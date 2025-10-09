from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure
import numpy as np
# eigsh is used like eig but computes a given number of eigenvalues
from scipy.sparse.linalg import eigsh

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
constrainedNodes = [1, 6, 12, 17]
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


plot_structure(elemList, nodeList + U*0)
