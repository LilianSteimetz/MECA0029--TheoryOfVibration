from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure
from geometry import constrainedNodes
from DFT import DFT, IDFT
import numpy as np
# eigsh is used like eig but computes a given number of eigenvalues
from scipy import linalg
import matplotlib.pyplot as plt

###########################################################################
#                                PART2                                    #
#                        Stationary response                              #
#                                                                         #
###########################################################################


def FRF(omega, eigvals, eigvecs, M, nModes):
    # compute the FRF matrix using the modes and natural frequencies already computed
    # also used to compute the mode dispacement method ( if nModes < size of M )

    omega_s = np.sqrt(eigvals)
    idx = np.argsort(omega_s)
    omega_s = omega_s[idx]
    eigvecs = eigvecs[:, idx]
    H = np.zeros((M.shape[0], M.shape[0]), dtype=complex)

    for s in range(nModes):
        mu_s = eigvecs[:, s].T @ M @ eigvecs[:, s]
        H += np.outer(eigvecs[:, s], eigvecs[:, s]) / \
            ((omega_s[s]**2 - omega**2) * mu_s)

    return H


def FRF_accel(omega, H, s, excitation_DOF):
    # computes the FRF in terms of acceleration,
    # based on the displacement FRF H

    a = -omega**2 * H @ s

    return np.abs(a[excitation_DOF])


def mode_acc(omega, eigvals, eigvecs, M, K, nModes):
    # computes an approx of the FRF matrix using the acceleration mode method

    omega_s = np.sqrt(eigvals)
    idx = np.argsort(omega_s)
    omega_s = omega_s[idx]
    eigvecs = eigvecs[:, idx]
    H = np.zeros((M.shape[0], M.shape[0]), dtype=complex)

    for s in range(nModes):
        mu_s = eigvecs[:, s].T @ M @ eigvecs[:, s]
        H += np.outer(eigvecs[:, s], eigvecs[:, s]) / \
            ((omega_s[s]**2 - omega**2) * mu_s * omega_s[s]**2)

    return np.linalg.inv(K) + (omega**2 * H)


""" 1 : Computation of the natural frequencies and mode shapes of the structure """

M_global, K_global = create_globalMass_and_globalStiffness()
eigvals, eigvecs = linalg.eig(K_global, M_global)


""" 2 : Parameters  """

# Number of mode shapes used in the different (MD and MA) methods
k = 1
# total nDof minus the ones constrained (clamped nodes)
nDOf_reduced = M_global.shape[0]


# Pulsation of the excitation
omega_exc = 2 * np.pi * 2.4

# Construction of s, the vector of amplitude of the excitation (at each DOF)
excitation_DOF = dofList[15, 1] - 1  # DOF in Y direction of node 16
excitation_DOF = excitation_DOF - 6 * \
    np.sum(constrainedNodes[:] < 16)

s = np.zeros(nDOf_reduced)
s[excitation_DOF] = 500  # [N]


""" 3: Computation of the frequency response functions and the time response """

# Transfer functions
H_frf = np.zeros((nDOf_reduced, nDOf_reduced), dtype=complex)
H_frf = FRF(omega_exc, eigvals, eigvecs, M_global, nDOf_reduced)

H_dis = np.zeros((nDOf_reduced, nDOf_reduced), dtype=complex)
H_dis = FRF(omega_exc, eigvals, eigvecs, M_global, k)

H_acc = np.zeros((nDOf_reduced, nDOf_reduced), dtype=complex)
H_acc = mode_acc(omega_exc, eigvals, eigvecs, M_global, K_global, k)

# Spatial response
x_frf = H_frf @ s
x_dis = H_dis @ s
x_acc = H_acc @ s


""" 4 : Computing RMS """

A = FRF_accel(omega_exc, H_frf, s, excitation_DOF)
rms = A / np.sqrt(2)
print(f"RMS value at node 16 in Y-direction: {rms:.4f} m/s²")


""" 5 : omputing the PSD """

fMaxPSD = 50  # [Hz]
omegaMaxPSD = 2*np.pi*fMaxPSD

N = 100 * fMaxPSD
omegaPSD = 2 * np.pi * np.linspace(0, fMaxPSD, N)
dOmega = omegaPSD[1] - omegaPSD[0]
PSDValues = np.zeros(np.shape(omegaPSD)[0])
psdNonZeroIdx = np.where(abs(omegaPSD - omega_exc) < 0.5 * dOmega)
PSDValues[psdNonZeroIdx] = 1/2 * A**2 / dOmega
psdIntegral = np.trapezoid(PSDValues, omegaPSD)
print(f"Value of the integration of PSD: {psdIntegral:.4f} m²/s⁴")
print(f"RMS value squared : {rms**2:.4f} m²/s⁴")


""" 6 : Plots"""

# Plotting the time response at node 16 in the Y-direction for each method
t = np.linspace(0, 1, 1000)

plt.plot(t, np.real(x_frf[excitation_DOF] *
         np.cos(omega_exc * t)), label='Exact solution')
plt.plot(t, np.real(x_dis[excitation_DOF] * np.cos(omega_exc * t)),
         label='Displacement mode approximation')
plt.plot(t, np.real(x_acc[excitation_DOF] * np.cos(omega_exc * t)),
         label='Acceleration mode approximation')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title('Response at node 16 in Y-direction')
plt.legend()
plt.grid(True)
plt.show()

# Plotting the PSD over the frequency
plt.plot(omegaPSD, PSDValues, label='PSD(w)')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('PSD (m²/s⁴)')
plt.title('Power Spectral Density of the acceleration at node 16 in Y-direction')
plt.legend()
plt.grid(True)
plt.show()


""" 7 : Convergence study of the methods"""
# done using the excitation of the project statement,
# the response is at the excitation node too.

# the true amplitude from the frf
NModesConv = 12

# check if it goes from 1 to 5 included
ModesVector = np.arange(1, NModesConv+1, 1, dtype=int)

frfAmplitudes = np.full(NModesConv, x_frf[excitation_DOF])

disMethAmplitudes = np.zeros(NModesConv)
accMethAmplitudes = np.zeros(NModesConv)

for NMode in ModesVector:

    H_dis = np.zeros((nDOf_reduced, nDOf_reduced), dtype=complex)
    H_dis = FRF(omega_exc, eigvals, eigvecs, M_global, NMode)
    x_dis = H_dis @ s
    disMethAmplitudes[NMode-1] = np.real(x_dis[excitation_DOF])

    H_acc = np.zeros((nDOf_reduced, nDOf_reduced), dtype=complex)
    H_acc = mode_acc(omega_exc, eigvals, eigvecs, M_global, K_global, NMode)
    x_acc = H_acc @ s
    accMethAmplitudes[NMode-1] = np.real(x_acc[excitation_DOF])


plt.plot(ModesVector, np.real(frfAmplitudes), label='FRF amplitude')
plt.plot(ModesVector, np.real(disMethAmplitudes),
         label='Mode displacement amplitude')
plt.plot(ModesVector, np.real(accMethAmplitudes),
         label='Mode acceleration amplitude')
plt.title("Amplitude at the excitation DOF")
plt.legend()
plt.xticks(ModesVector)
plt.show()

plt.plot(ModesVector, 100 * np.real(frfAmplitudes - disMethAmplitudes)/np.real(frfAmplitudes),
         label='Mode displacement error')
plt.plot(ModesVector, 100*np.real(frfAmplitudes - accMethAmplitudes)/np.real(frfAmplitudes),
         label='Mode acceleration error')
plt.title("Relative error of MD and MA at the excitation DOF (%)")
plt.yscale('log')
plt.legend()
plt.xticks(ModesVector)
plt.show()


disMethModeContrib = np.zeros(NModesConv)
accMethModeContrib = np.zeros(NModesConv)


for NMode in ModesVector:

    if NMode == 1:
        disMethModeContrib[0] = disMethAmplitudes[0]
        accMethModeContrib[0] = accMethAmplitudes[0]
    else:
        disMethModeContrib[NMode-1] = (disMethAmplitudes[NMode -
                                                         1] - disMethAmplitudes[NMode-2])
        accMethModeContrib[NMode-1] = (accMethAmplitudes[NMode -
                                                         1] - accMethAmplitudes[NMode-2])


plt.plot(ModesVector, np.real(disMethModeContrib)/np.real(np.sum(disMethModeContrib)),
         label='Mode displacement contrib')
plt.plot(ModesVector, np.real(accMethModeContrib)/np.real(np.sum(accMethModeContrib)),
         label='Mode acceleration contrib')
plt.title("Relative contribution of each mode to the response at the excitation DOF")
plt.legend()
plt.yscale('log')
plt.xticks(ModesVector)
plt.show()


""" Plotting the excitation in time"""

t = np.linspace(0, 5, 1000)
x_t = 500 * np.cos(omega_exc * t)

plt.plot(t, x_t, label='Excitation force')
plt.xlabel('Time (s)')
plt.ylabel('Force (N)')
plt.title('Excitation force at node 16 in Y-direction')
plt.legend()
plt.show()

# Plotting the DFT and IDFT functions to check they work well

f, X_f = DFT(t, x_t)
plt.stem(f, np.abs(X_f) / np.sum(X_f), label='DFT of the excitation force')
plt.xlim(0, 10)
plt.show()
