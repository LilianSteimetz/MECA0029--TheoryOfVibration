from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import dofList, vector_to_constrained
import numpy as np
from reductionMethods import guyansIrons, CraigBampton
from globalDampingMatrix import create_damping_matrix
from newmark import time_integ, time_integ_slides
import time
import matplotlib.pyplot as plt
import pandas as pd


# compute transient response using reduction methods


retained_nodes = [13, 16, 20]
retained_dofs = np.zeros(3*len(retained_nodes), dtype=int)
i = 0
for node in retained_nodes:
    for dof in range(3):
        retained_dofs[i] = dofList[node-1, dof]
        retained_dofs[i] -= 6 * np.sum(constrainedNodes[:] < node)
        i += 1


"""integration parameters"""
gamma_newmark = 1/2
beta_newmark = 1/4


"""transient excitation parameters:"""
f_duration = 0.01  # [s]
f_dof_idx = dofList[19, 1] - 1  # DOF in Y direction of node 20
f_dof = f_dof_idx - 6 * \
    np.sum(constrainedNodes[:] < 20)

f_vector = np.zeros(dofList.shape[0] * dofList.shape[1])
f_vector[f_dof_idx] = 7000  # [N]

f_vector = vector_to_constrained(f_vector)

timeStep = 1e-4  # [s]
integTime = 10  # [s]


M_global, K_global = create_globalMass_and_globalStiffness()
C_global = create_damping_matrix(M_global, K_global)

""" Full model transient response """
time_start = time.perf_counter()
q_full, q_dot_full, q_ddot_full, t = time_integ_slides(
    K_global, M_global, C_global, f_vector, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)
time_end = time.perf_counter()
time_full = time_end - time_start

""" Guyan-Irons transient response """
K_GI, M_GI, R_GI = guyansIrons(K_global, M_global, retained_dofs)
C_GI = R_GI.T @ C_global @ R_GI
f_vector_GI = R_GI.T @ f_vector
time_start = time.perf_counter()
q_GI, q_dot_GI, q_ddot_GI, t = time_integ_slides(
    K_GI, M_GI, C_GI, f_vector_GI, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)
time_end = time.perf_counter()
time_GI = time_end - time_start
q_GI = R_GI @ q_GI
q_dot_GI = R_GI @ q_dot_GI
q_ddot_GI = R_GI @ q_ddot_GI

""" Craig-Bampton transient response """
n_internal_modes = 7
K_CB, M_CB, R_CB = CraigBampton(
    K_global, M_global, retained_dofs, n_internal_modes)
C_CB = R_CB.T @ C_global @ R_CB
f_vector_CB = R_CB.T @ f_vector
time_start = time.perf_counter()
q_CB, q_dot_CB, q_ddot_CB, t = time_integ_slides(
    K_CB, M_CB, C_CB, f_vector_CB, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)
time_end = time.perf_counter()
time_CB = time_end - time_start
q_CB = R_CB @ q_CB
q_dot_CB = R_CB @ q_dot_CB
q_ddot_CB = R_CB @ q_ddot_CB

"""print time results"""
print(
    f"Computation times (full, GI, CB): {time_full:.6f}, {time_GI:.6f}, {time_CB:.6f}  s")


"""plotting results"""
"""
plt.plot(t, q_full[f_dof], label='Full model', linewidth=1)
plt.plot(t, q_GI[f_dof], label='Guyan-Irons', linewidth=1)
plt.plot(t, q_CB[f_dof], label='Craig-Bampton', linewidth=1)
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.xlim(0, integTime)
plt.legend()
plt.tight_layout()
plt.show()
"""


df = pd.DataFrame({
    'time': t,
    'Displacement full': q_full[f_dof, :],
    'Displacement GI': q_GI[f_dof, :],
    'Displacement CB': q_CB[f_dof, :],
})
df.to_csv('pt4/time_integ.csv', index=False)


"""

n_iter = 50
time_full = np.zeros(n_iter)
time_GI = np.zeros(n_iter)
time_CB = np.zeros(n_iter)
for i in range(n_iter):

    time_start = time.perf_counter()
    q_full, q_dot_full, q_ddot_full, t = time_integ_slides(
        K_global, M_global, C_global, f_vector, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)
    time_end = time.perf_counter()
    time_full[i] = time_end - time_start


for i in range(n_iter):
    time_start = time.perf_counter()
    q_CB, q_dot_CB, q_ddot_CB, t = time_integ_slides(
        K_CB, M_CB, C_CB, f_vector_CB, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)
    time_end = time.perf_counter()
    time_CB[i] = time_end - time_start


for i in range(n_iter):
    time_start = time.perf_counter()
    q_GI, q_dot_GI, q_ddot_GI, t = time_integ_slides(
        K_GI, M_GI, C_GI, f_vector_GI, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)
    time_end = time.perf_counter()
    time_GI[i] = time_end - time_start

t_full_avg = np.average(time_full)
t_GI_avg = np.average(time_GI)
t_CB_avg = np.average(time_CB)

std_dev_full = np.std(time_full)
std_dev_GI = np.std(time_GI)
std_dev_CB = np.std(time_CB)


print("average times : full, GI, CB")
print(f"{t_full_avg: .6f}, {t_GI_avg: .6f}, {t_CB_avg: .6f}")
print("std dev : full, GI, CB")
print(f"{std_dev_full:.6f},{std_dev_GI:.6f}, {std_dev_CB:.6f}")


plt.plot(np.arange(0, n_iter), time_full)
plt.show()
"""
