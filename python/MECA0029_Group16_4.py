import numpy as np
from scipy.sparse.linalg import eigsh
import time
from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import dofList
from reductionMethods import guyansIrons, CraigBampton
import matplotlib.pyplot as plt
from macMatrix import macMatrix, plotMacMatrix

# compute the natural frequencies and mode shapes to compare the reduction methods with full model


def eigvals_eigvectors_to_sorted_freqs_modes(eigvals, eigvecs):
    nat_freqs = np.sqrt(eigvals)/(2*np.pi)
    idx = np.argsort(nat_freqs)
    nat_freqs = nat_freqs[idx]
    eigvecs = eigvecs[:, idx]
    eigvecs = eigvecs/np.max(np.abs(eigvecs), axis=0)
    return nat_freqs, eigvecs


# compute the retained DOFs, in 1-based indexing, like dofList
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
eigvals, eigvecs = eigsh(K_global, k=desiredFreqNb,
                         M=M_global, sigma=0.0, which='LM')
time_end = time.perf_counter()
time_full = time_end - time_start

nat_freqs_full, eigvecs_full = eigvals_eigvectors_to_sorted_freqs_modes(
    eigvals, eigvecs)
print("Natural Frequencies, full model (Hz):")
print(nat_freqs_full)


""" Guyan-Irons 6 1st frequencies & modes"""
time_start = time.perf_counter()
K_GI, M_GI, R_GI = guyansIrons(K_global, M_global, retained_dofs)
time_end = time.perf_counter()
print("GI reduction time")
print(time_end - time_start)

time_start = time.perf_counter()
eigvals, eigvecs = eigsh(K_GI, k=desiredFreqNb,
                         M=M_GI, sigma=0.0, which='LM')
time_end = time.perf_counter()
time_GI = time_end - time_start

nat_freqs_GI, eigvecs_GI = eigvals_eigvectors_to_sorted_freqs_modes(
    eigvals, eigvecs)
eigvecs_GI = R_GI @ eigvecs_GI  # map back to full DOFs
print("Natural Frequencies, Guyan-Irons (Hz):")
print(nat_freqs_GI)

errors_GI = np.abs((nat_freqs_GI - nat_freqs_full)/nat_freqs_full)*100
errors_GI = errors_GI.round(3)

print("Guyan-Irons relative errors (%):")
print(errors_GI)


""" Craig-Bampton 6 1st frequencies & modes"""
n_internal_modes = 7
time_start = time.perf_counter()
K_CB, M_CB, R_CB = CraigBampton(
    K_global, M_global, retained_dofs, n_internal_modes)
time_end = time.perf_counter()
print("CB reduction time")
print(time_end - time_start)
time_start = time.perf_counter()
eigvals, eigvecs = eigsh(K_CB, k=desiredFreqNb,
                         M=M_CB, sigma=0.0, which='LM')
time_end = time.perf_counter()
time_CB = time_end - time_start
nat_freqs_CB, eigvecs_CB = eigvals_eigvectors_to_sorted_freqs_modes(
    eigvals, eigvecs)
eigvecs_CB = R_CB @ eigvecs_CB  # map back to full DOFs
print("Natural Frequencies, Craig-Bampton (Hz):")
print(nat_freqs_CB)

errors_CB = np.abs((nat_freqs_CB - nat_freqs_full)/nat_freqs_full)*100
errors_CB = errors_CB.round(3)
print("Craig-Bampton relative errors (%):")
print(errors_CB)


print(
    f"Computation times (full, GI, CB): {time_full:.6f}, {time_GI:.6f}, {time_CB:.6f}  s")


""" mac matrix for GI vs Full model"""

mac_matrix = macMatrix(eigvecs_full, eigvecs_GI)   # choose GI or CB

plotMacMatrix(mac_matrix, title="",
              show=True, save_path="presentation/MAC_Full_vs_GI.svg", red='GI')

"""Plot mac matrix for CB vs Full model"""
mac_matrix = macMatrix(eigvecs_full, eigvecs_CB)   # choose GI or CB

plotMacMatrix(mac_matrix, title="",
              show=True, save_path="presentation/MAC_Full_vs_CB.svg", red='CB')
