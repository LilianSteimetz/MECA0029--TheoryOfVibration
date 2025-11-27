from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure, vector_to_constrained
from geometry import constrainedNodes
from globalDampingMatrix import create_damping_matrix
from newmark import time_integ, time_integ_slides
from DFT import DFT, IDFT

import numpy as np
import matplotlib.pyplot as plt

###########################################################################
#                                PART3                                    #
#                         Transient response                              #
#                                                                         #
###########################################################################

""" 1 : integration parameters"""
gamma_newmark = 3/4
beta_newmark = min(1.1*(1/4 * (gamma_newmark + 1/2)**2), 1)
timeStep = 0.00001  # [s]
integTime = 0.02  # [s]

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

""" 4 : Time integration """
q, q_dot, q_ddot, t = time_integ_slides(
    K_global, M_global, C_global, f_vector, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)

""" 5 : Plot time response """
plt.plot(t, q[f_dof], linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.xlim(0, integTime)
plt.legend()
plt.tight_layout()
plt.show()

""" 6 : DFT computation"""
y = q[f_dof]                    # shape (N_steps,)
y = y - np.mean(y)                  # remove DC
N = y.size
T = t[-1] - t[0]                    # total duration
dt = t[1] - t[0]                    # sampling period
fs = 1/dt                           # sampling frequency

Y = np.fft.rfft(y)                         # one-sided FFT
freqs = np.fft.rfftfreq(N, dt)             # matching frequency axis

""" 7 : DFT plot """

plt.plot(freqs, np.abs(Y))
plt.xlabel('Frequency [Hz]')
plt.vlines(x=[2.573, 6.844, 10.855, 12.916, 13.866, 14.877], ymin=0, ymax=max(np.abs(
    Y)), colors='r', linestyles='dashed', linewidth=0.5, label='Natural Frequencies')
plt.ylabel('|Y|')
plt.xlim(0, 30)
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.show()
