import numpy as np
from scipy.sparse.linalg import eigsh
from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure, vector_to_constrained
from geometry import constrainedNodes
from matplotlib import pyplot as plt

def create_damping_matrix(M, K, damping1=0.01, damping2=0.01):
    # damping ratios
    epsilon1 = damping1
    epsilon2 = damping2
    epsilon = np.array([epsilon1, epsilon2])

    # compute first two natural frequencies (pulsations)
    omegas, _ = eigsh(K, k=2,
                      M=M, sigma=0.0, which='LM')
    omegas = np.sqrt(omegas)
    idx = np.argsort(omegas)
    omegas = omegas[idx]
    omega1, omega2 = omegas[0], omegas[1]
    # solve for alpha and beta, using the modal damping assumption
    # with 2 modal damping ratios given
    Omega_mat = np.array([[omega1, 1/omega1], [omega2, 1/omega2]])
    alpha, beta = 2 * np.linalg.solve(Omega_mat, epsilon)
    print(f"Computed Rayleigh damping coefficients: alpha = {alpha:.6e}, beta = {beta:.6e}")

    C = alpha * K + beta * M

    return C

M_global, K_global = create_globalMass_and_globalStiffness()
damp_mat = create_damping_matrix(M_global, K_global)
eigvals, eigvecs = eigsh(K_global, k=desiredFreqNb,
                         M=M_global, sigma=0.0, which='LM')

nat_freqs = np.sqrt(eigvals)/(2*np.pi)
idx = np.argsort(nat_freqs)
nat_freqs = nat_freqs[idx]
eigvecs = eigvecs[:, idx]
eigvecs = eigvecs/np.max(np.abs(eigvecs), axis=0)
epsilon = np.zeros(desiredFreqNb)

for i in range(6):
    beta_i = eigvecs[:, i].T @ damp_mat @ eigvecs[:, i]
    omega_i = np.sqrt(eigvals[i])
    mu_i =  eigvecs[:, i].T @ M_global@ eigvecs[:, i]
    epsilon[i] = (beta_i) / (2 * mu_i * omega_i)
    
epsilon_2 = np.zeros(desiredFreqNb)
alpha = 3.380138e-04
beta = 2.349904e-01
for i in range(6):
    epsilon_2[i] = ((alpha * np.sqrt(eigvals[i])) + (beta / np.sqrt(eigvals[i]))) / 2
    
plt.plot(np.sqrt(eigvals), epsilon, 'o-')
plt.show()

print(f"First six damping ratios: {epsilon[:]}")
print(f"First six damping ratios: {epsilon_2[:]}")