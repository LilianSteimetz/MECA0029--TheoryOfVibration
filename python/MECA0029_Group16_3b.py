from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure, vector_to_constrained
from geometry import constrainedNodes
from globalDampingMatrix import create_damping_matrix
from newmark import time_integ, time_integ_slides
from DFT import DFT, IDFT

import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd


# integration parameters
gamma_newmark = 0.5
beta_newmark = 1/4

timeStep = 5e-3  # [s]
integTime = 5  # [s]

# transient excitation parameters:
f_duration = 0.01  # [s]
f_dof_idx = dofList[19, 1] - 1  # DOF in Y direction of node 20
f_dof = f_dof_idx - 6 * \
    np.sum(constrainedNodes[:] < 20)

f_vector = np.zeros(dofList.shape[0] * dofList.shape[1])
f_vector[f_dof_idx] = 7000  # [N]

f_vector = vector_to_constrained(f_vector)


M_global, K_global = create_globalMass_and_globalStiffness()
C_global = create_damping_matrix(M_global, K_global)


q, q_dot, q_ddot, t = time_integ_slides(
    K_global, M_global, C_global, f_vector, f_duration, timeStep, integTime, gamma=gamma_newmark, beta=beta_newmark)


df = pd.DataFrame({
    'time': t,
    'Displacement': q[f_dof, :],
    'velocity': q_dot[f_dof, :],
    'Acceleration': q_ddot[f_dof, :],
})
df.to_csv(f'pt3/time_integ_T{integTime}_h{timeStep}.csv', index=False)
