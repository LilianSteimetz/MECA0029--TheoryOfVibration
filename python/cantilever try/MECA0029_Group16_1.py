import scipy.linalg as la
from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure
import numpy as np
# eigsh is used like eig but computes a given number of eigenvalues
from scipy.sparse.linalg import eigsh


M_global, K_global = create_globalMass_and_globalStiffness()

eigvals, eigvecs = eigsh(K_global, k=6, M=M_global, sigma=0.0, which='LM')

nat_freqs = np.sqrt(eigvals)/(2*np.pi)

# Sort natural frequencies and mode shapes
idx = np.argsort(nat_freqs)
nat_freqs = nat_freqs[idx]
eigvecs = eigvecs[:, idx]
eigvecs = eigvecs/np.max(np.abs(eigvecs), axis=0)
print("Natural Frequencies (Hz):")
print(nat_freqs)


# Add the constrained DOFs back (0 displacement)
constrainedDOFs = []
constrainedNodes = [1]
for i in range(len(constrainedNodes)):
    node = constrainedNodes[i]
    constrainedDOFs.extend(dofList[node-1, :] - 1)
constrainedDOFs = constrainedDOFs[::-1]
"""
for i in range(len(constrainedDOFs)):
    eigvecs = np.insert(
        eigvecs, constrainedDOFs[i], 0, axis=0)
"""
nDOF = np.max(dofList)   # total number of DOFs
nNodes = nodeList.shape[0]
nModes = eigvecs.shape[1]
eigvecs_full = np.zeros((nDOF, nModes))
allDOFs = np.arange(nDOF)
freeDOFs = np.setdiff1d(allDOFs, constrainedDOFs)
eigvecs_full[freeDOFs, :] = eigvecs
eigvecs = eigvecs_full


# Visualization of a mode shape
mode_idx = plot_mode  # mode to visualize


U = np.zeros((nNodes, 3))  # x, y, z displacement per node
for i in range(nNodes):
    U[i, 0] = eigvecs[dofList[i, 0]-1, mode_idx]  # x
    U[i, 1] = eigvecs[dofList[i, 1]-1, mode_idx]  # y
    U[i, 2] = eigvecs[dofList[i, 2]-1, mode_idx]  # z

plot_structure(elemList, nodeList + U*0)

"""
# use your assembled global matrices (before reduction or after? use the reduced ones you feed the solver)
K = K_global.copy()
M = M_global.copy()


def report_matrix(name, A):
    print(f"\n{name}: shape={A.shape}, dtype={A.dtype}")
    print("symmetry error (max abs A - A.T):", np.max(np.abs(A - A.T)))
    vals = la.eigvals(A)
    reals = np.real(vals)
    print("eigs real part: min, max:", reals.min(), reals.max())


report_matrix("K", K)
report_matrix("M", M)

# trace and diag scale
print("\ntrace(K) =", np.trace(K))
print("trace(M) =", np.trace(M))
print("diag(K)[:10] =", np.diag(K)[:10])
print("diag(M)[:10] =", np.diag(M)[:10])

# check definiteness (use small eigenvalue check)
k_eigs = la.eigvalsh(K, subset_by_index=[0, min(9, K.shape[0]-1)])
m_eigs = la.eigvalsh(M, subset_by_index=[0, min(9, M.shape[0]-1)])
print("\nsmallest few eigenvalues K:", k_eigs[:6])
print("smallest few eigenvalues M:", m_eigs[:6])
"""
