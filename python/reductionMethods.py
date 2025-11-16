import numpy as np
from scipy.sparse.linalg import eigsh


def guyansIrons(K_full, M_full, retained_dofs):
    retained_dofs = np.array(retained_dofs) - 1  # zero indexing

    N_dofs_full = K_full.shape[0]
    all_dofs = np.arange(N_dofs_full)
    # dofs to condense = all dofs - retained dofs
    condensed_dofs = np.setdiff1d(all_dofs, retained_dofs)
    # partition K and M
    K_cr = K_full[np.ix_(condensed_dofs, retained_dofs)]
    K_cc = K_full[np.ix_(condensed_dofs, condensed_dofs)]

    # reduction matrix
    #  x_full (n_dof_full x 1) = R x_reduced (n_retained_dofs x 1) => R ( n_dofs_full x n_retained_dofs)
    R = np.zeros((N_dofs_full, len(retained_dofs)))
    R[np.ix_(retained_dofs, np.arange(len(retained_dofs)))
      ] = np.eye(len(retained_dofs))

    R[np.ix_(condensed_dofs, np.arange(len(retained_dofs)))
      ] = np.linalg.solve(-K_cc, K_cr)
    K_bar = R.T @ K_full @ R
    M_bar = R.T @ M_full @ R

    return K_bar, M_bar, R


def CraigBampton(K_full, M_full, boundary_dofs, n_modes_internal):
    boundary_dofs = np.array(boundary_dofs) - 1  # zero indexing

    N_dofs_full = K_full.shape[0]
    N_boundary_dofs = len(boundary_dofs)

    all_dofs = np.arange(N_dofs_full)
    # dofs to condense = all dofs - boundary dofs
    internal_dofs = np.setdiff1d(all_dofs, boundary_dofs)

    # partition K and M
    K_ib = K_full[np.ix_(internal_dofs, boundary_dofs)]
    K_ii = K_full[np.ix_(internal_dofs, internal_dofs)]

    M_ii = M_full[np.ix_(internal_dofs, internal_dofs)]

    # reduction matrix
    R = np.zeros((N_dofs_full, N_boundary_dofs + n_modes_internal))
    R_11 = np.eye(len(boundary_dofs))
    R_21 = np.linalg.solve(-K_ii, K_ib)
    # compute the n_modes_internal lowest eigenmodes of the internal dofs
    eigvals, eigvecs = eigsh(K_ii, k=n_modes_internal,
                             M=M_ii, sigma=0.0, which='LM')

    R_22 = eigvecs

    R[np.ix_(boundary_dofs, np.arange(N_boundary_dofs))] = R_11
    R[np.ix_(internal_dofs, np.arange(N_boundary_dofs))] = R_21
    R[np.ix_(internal_dofs, np.arange(N_boundary_dofs,
             N_boundary_dofs + n_modes_internal))] = R_22

    K_bar = R.T @ K_full @ R
    M_bar = R.T @ M_full @ R

    return K_bar, M_bar, R
