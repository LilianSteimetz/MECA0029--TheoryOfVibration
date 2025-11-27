import numpy as np
from scipy.sparse.linalg import eigsh


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

    C = alpha * K + beta * M

    return C
