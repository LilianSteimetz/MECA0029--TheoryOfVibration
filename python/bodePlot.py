from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure
from geometry import constrainedNodes
from DFT import DFT, IDFT
import numpy as np
# eigsh is used like eig but computes a given number of eigenvalues
from scipy import linalg
import matplotlib.pyplot as plt


def FRF(omega, eigvals, eigvecs, M, nModes):
    # compute the FRF matrix using the modes and natural frequencies already computed
    # also used to compute the mode dispacement method ( if nModes < size of M )

    omega_s = np.sqrt(eigvals)
    idx = np.argsort(omega_s)
    omega_s = omega_s[idx]
    eigvecs = eigvecs[:, idx]
    H = np.zeros((M.shape[0], M.shape[0]), dtype=complex)

    for s in range(nModes):
        mu_s = eigvecs[:, s].T @ M @ eigvecs[:, s]
        H += np.outer(eigvecs[:, s], eigvecs[:, s]) / \
            ((omega_s[s]**2 - omega**2) * mu_s)

    return H


f_dof = dofList[7, 1] - 1  # DOF in Y direction of node 20
f_dof = f_dof - 6 * \
    np.sum(np.where(constrainedNodes[:] < 7))

omegas = np.linspace(0.1, 2 * np.pi * 20, 300)
frf = np.zeros((len(omegas), 1), dtype=complex)
for i, om in enumerate(omegas):
    print(f'Computing FRF at frequency {om/(2*pi):.2f} Hz', end='\r')

    M_global, K_global = create_globalMass_and_globalStiffness()

    eigvals, eigvecs = linalg.eigh(K_global, M_global)

    H = FRF(om, eigvals, eigvecs, M_global, desiredFreqNb)

    frf[i] = H[f_dof, f_dof]

plt.plot(omegas / (2 * np.pi), 20 * np.log10(np.abs(frf)), linewidth=2)
plt.xlabel('Frequency [Hz]')
plt.ylabel('FRF Magnitude [dB]')
plt.show()
