import numpy as np
from scipy.sparse.linalg import eigsh
import time
from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import dofList, vector_to_constrained
from reductionMethods import guyansIrons, CraigBampton
import matplotlib.pyplot as plt
from macMatrix import macMatrix, plotMacMatrix
from newmark import time_integ, time_integ_slides
from globalDampingMatrix import create_damping_matrix


###########################################################################
#                                PART4                                    #
#                          Reduction methods                              #
#                           Modal Analysis                                #
###########################################################################


def eigvals_eigvectors_to_sorted_freqs_modes(eigvals, eigvecs):
    nat_freqs = np.sqrt(eigvals)/(2*np.pi)
    idx = np.argsort(nat_freqs)
    nat_freqs = nat_freqs[idx]
    eigvecs = eigvecs[:, idx]
    eigvecs = eigvecs/np.max(np.abs(eigvecs), axis=0)
    return nat_freqs, eigvecs


""" 1 : Computation of retained DOFs (1-indexed)"""
retained_nodes = [13, 16, 20]
retained_dofs = np.zeros(3*len(retained_nodes), dtype=int)
i = 0
for node in retained_nodes:
    for dof in range(3):
        retained_dofs[i] = dofList[node-1, dof]
        retained_dofs[i] -= 6 * np.sum(constrainedNodes[:] < node)
        i += 1

""" 2 : Global Matrices computation """
M_global, K_global = create_globalMass_and_globalStiffness()

""" 3 : Full Model modal analysis"""
eigvals, eigvecs = eigsh(K_global, k=desiredFreqNb,
                         M=M_global, sigma=0.0, which='LM')

nat_freqs_full, eigvecs_full = eigvals_eigvectors_to_sorted_freqs_modes(
    eigvals, eigvecs)


""" 4 : Guyan-Irons modal analysis"""
K_GI, M_GI, R_GI = guyansIrons(K_global, M_global, retained_dofs)

eigvals, eigvecs = eigsh(K_GI, k=desiredFreqNb,
                         M=M_GI, sigma=0.0, which='LM')

nat_freqs_GI, eigvecs_GI = eigvals_eigvectors_to_sorted_freqs_modes(
    eigvals, eigvecs)
eigvecs_GI = R_GI @ eigvecs_GI  # map back to full DOFs


""" 5 : Craig-Bampton modal analysis"""
n_internal_modes = 7
K_CB, M_CB, R_CB = CraigBampton(
    K_global, M_global, retained_dofs, n_internal_modes)
eigvals, eigvecs = eigsh(K_CB, k=desiredFreqNb,
                         M=M_CB, sigma=0.0, which='LM')
nat_freqs_CB, eigvecs_CB = eigvals_eigvectors_to_sorted_freqs_modes(
    eigvals, eigvecs)
eigvecs_CB = R_CB @ eigvecs_CB  # map back to full DOFs


""" 6 : prints of the results """
print("Natural Frequencies, full model (Hz):")
print(nat_freqs_full)

print("Natural Frequencies, Guyan-Irons (Hz):")
print(nat_freqs_GI)

print("Natural Frequencies, Craig-Bampton (Hz):")
print(nat_freqs_CB)


errors_GI = np.abs((nat_freqs_GI - nat_freqs_full)/nat_freqs_full)*100
errors_CB = np.abs((nat_freqs_CB - nat_freqs_full)/nat_freqs_full)*100

print(f"Guyan-Irons relative errors (%): {errors_GI.round(3)}")

print(f"Craig-Bampton relative errors (%): {errors_CB.round(3)}")


""" 7 : MAC matrices """

mac_matrix = macMatrix(eigvecs_full, eigvecs_GI)   # full vs gi

plotMacMatrix(mac_matrix, title="",
              show=True, red='GI')

mac_matrix = macMatrix(eigvecs_full, eigvecs_CB)   # full vs CB

plotMacMatrix(mac_matrix, title="",
              show=True, red='CB')


###########################################################################
#                                PART4                                    #
#                          Reduction methods                              #
#                Transient response, time integration                     #
###########################################################################

""" 1 : integration parameters"""
gamma_newmark = 1/2
beta_newmark = 1/4
timeStep = 1e-4  # [s]
integTime = 10  # [s]


""" 2 : transient excitation parameters"""
f_duration = 0.01  # [s]
f_dof_idx = dofList[19, 1] - 1  # DOF in Y direction of node 20
f_dof = f_dof_idx - 6 * \
    np.sum(constrainedNodes[:] < 20)

f_vector = np.zeros(dofList.shape[0] * dofList.shape[1])
f_vector[f_dof_idx] = 7000  # [N]

f_vector = vector_to_constrained(f_vector)

""" 3 : Global Matrices creation """
M_global, K_global = create_globalMass_and_globalStiffness()
C_global = create_damping_matrix(M_global, K_global)


""" 4 : Full model transient response """
q_full, q_dot_full, q_ddot_full, t = time_integ_slides(
    K_global, M_global, C_global, f_vector, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)

""" 5 : Guyan-Irons transient response """
K_GI, M_GI, R_GI = guyansIrons(K_global, M_global, retained_dofs)
C_GI = R_GI.T @ C_global @ R_GI
f_vector_GI = R_GI.T @ f_vector

q_GI, q_dot_GI, q_ddot_GI, t = time_integ_slides(
    K_GI, M_GI, C_GI, f_vector_GI, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)

q_GI, q_dot_GI, q_ddot_GI = R_GI @ q_GI, R_GI @ q_dot_GI, R_GI @ q_ddot_GI


""" 6 : Craig-Bampton transient response """
n_internal_modes = 7
K_CB, M_CB, R_CB = CraigBampton(
    K_global, M_global, retained_dofs, n_internal_modes)
C_CB = R_CB.T @ C_global @ R_CB
f_vector_CB = R_CB.T @ f_vector
q_CB, q_dot_CB, q_ddot_CB, t = time_integ_slides(
    K_CB, M_CB, C_CB, f_vector_CB, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)

q_CB, q_dot_CB, q_ddot_CB = R_CB @ q_CB, R_CB @ q_dot_CB, R_CB @ q_ddot_CB


""" 7 : Plot time responses """

plt.plot(t, q_full[f_dof], label='Full model', linewidth=1)
plt.plot(t, q_GI[f_dof], label='Guyan-Irons', linewidth=1)
plt.plot(t, q_CB[f_dof], label='Craig-Bampton', linewidth=1)
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.xlim(0, integTime)
plt.legend()
plt.tight_layout()
plt.show()


""" 8 : DFT computation """

idx = np.where((t >= 0) & (t <= 5))[0]

t = t[idx]
q_full = q_full[f_dof, idx]              # select columns corresponding to time
q_GI = q_GI[f_dof, idx]
q_CB = q_CB[f_dof, idx]

N = q_full.size
T = t[-1] - t[0]
dt = t[1] - t[0]
fs = 1/dt

Q_full = np.fft.rfft(q_full)
Q_GI = np.fft.rfft(q_GI)
Q_CB = np.fft.rfft(q_CB)

freqs = np.fft.rfftfreq(N, dt)

""" 9 : Plot the DFT"""

# plot pre-processing : cut to avoid the 0 in f = 0, take absolute value
Q_full = np.abs(Q_full[1:])
Q_GI = np.abs(Q_GI[1:])
Q_CB = np.abs(Q_CB[1:])
freqs = freqs[1:]


plt.plot(freqs, Q_full, label="full model")
plt.plot(freqs, Q_GI, label="Guyan-Irons")
plt.plot(freqs, Q_CB, label="Craig-Bampton")
plt.yscale("log")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.xlim(0, 30)
plt.legend
plt.show()
