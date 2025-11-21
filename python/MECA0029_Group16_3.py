from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure, vector_to_constrained
from geometry import constrainedNodes
from globalDampingMatrix import create_damping_matrix
from newmark import time_integ, time_integ_slides
from DFT import DFT, IDFT

import numpy as np
# eigsh is used like eig but computes a given number of eigenvalues
import matplotlib.pyplot as plt
import time


# integration parameters
gamma_newmark = 3/4
beta_newmark = min(1.1*(1/4 * (gamma_newmark + 1/2)**2), 1)

# transient excitation parameters:
f_duration = 0.01  # [s]
f_dof_idx = dofList[19, 1] - 1  # DOF in Y direction of node 20
f_dof = f_dof_idx - 6 * \
    np.sum(constrainedNodes[:] < 16)

f_vector = np.zeros(dofList.shape[0] * dofList.shape[1])
f_vector[f_dof_idx] = 7000  # [N]

f_vector = vector_to_constrained(f_vector)

timeStep = 0.00001  # [s]
integTime = 0.02  # [s]


M_global, K_global = create_globalMass_and_globalStiffness()
C_global = create_damping_matrix(M_global, K_global)
"""
q, q_dot, q_ddot, t = time_integ(
    K_global, M_global, C_global, f_vector, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)
# plt.plot(t, q[f_dof_idx], linewidth=2)
q2, q_dot, q_ddot, t = time_integ(
    K_global, M_global, C_global, f_vector, f_duration, timeStep/10, integTime, gamma=gamma_newmark, beta=beta_newmark)
# plt.plot(t, q[f_dof_idx], linewidth=2)
plt.plot(t, q[f_dof_idx] - q2[f_dof_idx][:len(q[f_dof_idx])],
         label='dt=%.4f s' % timeStep, linewidth=2)
"""
q, q_dot, q_ddot, t = time_integ(
    K_global, M_global, C_global, f_vector, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)
plt.plot(t, q[f_dof], linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.xlim(0, integTime)
plt.legend()
plt.tight_layout()
plt.savefig(
    f'plots/pt3_t_{elemPerBar}EPB_{integTime}T_{timeStep}dt.pdf', dpi=300)

# plt.vlines(x=f_duration, ymin=min(q[f_dof]), ymax=max(q[f_dof]), colors= 'r', linestyles = 'dashed', label = 'End of excitation')
plt.close()
"""
dft_q = np.fft.fft(q[f_dof])
f, dft_q = DFT(t, q[f_dof])
plt.plot(f, np.abs(dft_q), linewidth=2)
plt.xlabel('Frequency [Hz]')
plt.ylabel('relative DFT amplitude []')
plt.xlim(0, 30)
plt.show()
"""

y = q[f_dof]                    # shape (N_steps,)
y = y - np.mean(y)                  # remove DC
N = y.size
T = t[-1] - t[0]                    # total duration
dt = t[1] - t[0]                    # sampling period
fs = 1/dt                           # sampling frequency

Y = np.fft.rfft(y)                         # one-sided FFT
freqs = np.fft.rfftfreq(N, dt)             # matching frequency axis

plt.plot(freqs, np.abs(Y))
plt.xlabel('Frequency [Hz]')
plt.vlines(x=[2.573, 6.844, 10.855, 12.916, 13.866, 14.877], ymin=0, ymax=max(np.abs(
    Y)), colors='r', linestyles='dashed', linewidth=0.5, label='Natural Frequencies')
plt.ylabel('|Y|')
plt.xlim(0, 30)
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig(
    f'plots/pt3_f_{elemPerBar}EPB_{integTime}T_{timeStep}dt.pdf', dpi=300)
# plt.show()
plt.close()
