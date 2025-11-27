import numpy as np
from scipy.linalg import eigh
import time
from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import dofList
from reductionMethods import guyansIrons, CraigBampton
import matplotlib.pyplot as plt
from macMatrix import macMatrix, plotMacMatrix
import pandas as pandas

# compute the natural frequencies and mode shapes to compare the reduction methods with full model


def eigvals_eigvectors_to_sorted_freqs_modes(eigvals, eigvecs):
    nat_freqs = np.sqrt(eigvals)/(2*np.pi)
    idx = np.argsort(nat_freqs)
    nat_freqs = nat_freqs[idx]
    eigvecs = eigvecs[:, idx]
    eigvecs = eigvecs/np.max(np.abs(eigvecs), axis=0)
    return nat_freqs, eigvecs


"""Retained DOFs"""
retained_nodes = [13, 16, 20]
retained_dofs = np.zeros(3*len(retained_nodes), dtype=int)
i = 0
for node in retained_nodes:
    for dof in range(3):
        retained_dofs[i] = dofList[node-1, dof]
        retained_dofs[i] -= 6 * np.sum(constrainedNodes[:] < node)
        i += 1


M_global, K_global = create_globalMass_and_globalStiffness()

""" full model 6 1st frequencies & modes"""
time_start = time.perf_counter()
eigvals, eigvecs = eigh(K_global, M_global, subset_by_index=[
                        0, desiredFreqNb-1])
time_end = time.perf_counter()
time_full = time_end - time_start

nat_freqs_full, eigvecs_full = eigvals_eigvectors_to_sorted_freqs_modes(
    eigvals, eigvecs)


""" Guyan-Irons 6 1st frequencies & modes"""
K_GI, M_GI, R_GI = guyansIrons(K_global, M_global, retained_dofs)
time_start = time.perf_counter()
eigvals, eigvecs = eigh(K_GI, M_GI, subset_by_index=[
    0, desiredFreqNb-1])
time_end = time.perf_counter()
time_GI = time_end - time_start

nat_freqs_GI, eigvecs_GI = eigvals_eigvectors_to_sorted_freqs_modes(
    eigvals, eigvecs)
eigvecs_GI = R_GI @ eigvecs_GI  # map back to full DOFs


errors_GI = np.abs((nat_freqs_GI - nat_freqs_full)/nat_freqs_full)*100
errors_GI = errors_GI.round(3)


""" Craig-Bampton 6 1st frequencies & modes"""
max_n_modes = 50
errors_CB = np.zeros((max_n_modes, 6))
times_reduction = np.zeros(max_n_modes)
times_system = np.zeros(max_n_modes)
nat_freqs_CB = np.zeros((max_n_modes, 6))

for i in range(max_n_modes):
    print(f'Number of retained internal modes: {i+1}')
    n_internal_modes = i + 1

    time_start = time.perf_counter()

    K_CB, M_CB, R_CB = CraigBampton(
        K_global, M_global, retained_dofs, n_internal_modes)

    time_end = time.perf_counter()
    times_reduction[i] = time_end - time_start

    time_start = time.perf_counter()
    eigvals, eigvecs = eigh(K_CB, M_CB, subset_by_index=[
                            0, desiredFreqNb-1])
    time_end = time.perf_counter()
    times_system[i] = time_end - time_start

    nat_freqs_CB[i, :], eigvecs_CB = eigvals_eigvectors_to_sorted_freqs_modes(
        eigvals, eigvecs)
    eigvecs_CB = R_CB @ eigvecs_CB  # map back to full DOFs
    errors_CB[i, :] = np.abs(
        (nat_freqs_CB[i, :] - nat_freqs_full)/nat_freqs_full)*100

df = pandas.DataFrame({
    'Retained internal modes': np.arange(1, max_n_modes + 1),
    'Reduction time (s)': times_reduction,
    'Eigenproblem time (s)': times_system,
    **{f'CB Mode {i+1} error (%)': errors_CB[:, i] for i in range(6)},
    **{f'CB Mode {i+1} freq': nat_freqs_CB[:, i] for i in range(6)},
    **{f'GI Mode {i+1} freq': np.full(max_n_modes, nat_freqs_GI[i]) for i in range(6)},
    **{f'GI Mode {i+1} error (%)': np.full(max_n_modes, errors_GI[i]) for i in range(6)},
    **{f'Full Mode {i+1} freq': np.full(max_n_modes, nat_freqs_full[i]) for i in range(6)},
})

df.to_csv('pt4/CB_errors.csv', index=False)
