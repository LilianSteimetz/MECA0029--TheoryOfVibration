from constants import *
from globalMassStiffMatrices import create_globalMass_and_globalStiffness
from mesh import elemList, nodeList, dofList, plot_structure
import numpy as np
# eigsh is used like eig but computes a given number of eigenvalues
from scipy import linalg
import matplotlib.pyplot as plt

def FRF(omega, eigvals, eigvecs, M, nDOF_reduced):
    
    omega_s = np.sqrt(eigvals) 
    idx = np.argsort(omega_s)
    omega_s = omega_s[idx]
    eigvecs = eigvecs[:, idx]
    H = np.zeros((M.shape[0], M.shape[0]), dtype=complex)
    
    for s in range(nDOF_reduced): 
        mu_s = eigvecs[:,s].T @ M @ eigvecs[:,s]
        H += np.outer(eigvecs[:,s], eigvecs[:,s]) / ((omega_s[s]**2 - omega**2) * mu_s)
    
    return H

def mode_acc(omega, eigvals, eigvecs, M, K, k):
    
    omega_s = np.sqrt(eigvals) 
    idx = np.argsort(omega_s)
    omega_s = omega_s[idx]
    eigvecs = eigvecs[:, idx]
    H = np.zeros((M.shape[0], M.shape[0]), dtype=complex)
    
    for s in range(k): 
        mu_s = eigvecs[:,s].T @ M @ eigvecs[:,s]
        H += np.outer(eigvecs[:,s], eigvecs[:,s]) / ((omega_s[s]**2 - omega**2) * mu_s * omega_s[s]**2)
    
    return np.linalg.inv(K) + (omega**2 * H)

def FRF_accel(omega, H, s, excitation_DOF):
    
    a = -omega**2 * H @ s
    
    return np.abs(a[excitation_DOF])
    
""" Computation of the natural frequencies and mode shapes of the structure """
M_global, K_global = create_globalMass_and_globalStiffness()
eigvals, eigvecs = linalg.eig(K_global, M_global)

""" Number of mode shapes used in the different methods """
nDOf_reduced = M_global.shape[0]
k = 1 

""" Excitation load parameters """
omega_exc = 2 * np.pi * 2.4  # Excitation frequency of node 16
excitation_DOF = dofList[15, 1] - 1  
excitation_DOF = excitation_DOF - 18 # 3 clamped nodes (before node 16) with 6 DOF's each

s = np.zeros(nDOf_reduced)
s[excitation_DOF] = 500 # Amplitude of the force applied at node 16 in the Y-direction

""" Computation of the frequency response functions and the time response """
H_frf = np.zeros((nDOf_reduced, nDOf_reduced), dtype=complex)
H_frf = FRF(omega_exc, eigvals, eigvecs, M_global, nDOf_reduced)

H_dis = np.zeros((nDOf_reduced, nDOf_reduced), dtype=complex)
H_dis = FRF(omega_exc, eigvals, eigvecs, M_global, k)

H_acc = np.zeros((nDOf_reduced, nDOf_reduced), dtype=complex)
H_acc = mode_acc(omega_exc, eigvals, eigvecs, M_global, K_global, k)

x_frf = H_frf @ s
x_dis = H_dis @ s
x_acc = H_acc @ s

"""Computing PSD and RMS """
A = FRF_accel(omega_exc, H_frf, s, excitation_DOF)

rms  = A / np.sqrt(2)
print(f"RMS value at node 16 in Y-direction: {rms:.4f} m/s²")

Fs = 200
T = 20
N = int(T * Fs)
freq = np.fft.rfftfreq(N, 1/Fs)
dw = 2*np.pi * (freq[1] - freq[0])
psd = np.zeros(len(freq))
idx = np.argmin(np.abs((2*np.pi*freq) - omega_exc))
psd[idx] = A**2 / (2* dw)

psd_value = np.trapezoid(psd, 2*np.pi*freq) # Value of the integration of PSD should be = rms**2 (verification)
print(f"Value of the integration of PSD: {psd_value:.4f} m²/s⁴")
print(f"RMS value squared : {rms**2:.4f} m²/s⁴")



""" Plotting the time response at node 16 in the Y-direction for each method"""
t = np.linspace(0, 1, 1000)

plt.plot(t, np.real(x_frf[excitation_DOF] * np.cos(omega_exc * t)), label='Exact solution')
plt.plot(t, np.real(x_dis[excitation_DOF] * np.cos(omega_exc * t)), label='Displacement mode approximation')
plt.plot(t, np.real(x_acc[excitation_DOF] * np.cos(omega_exc * t)), label='Acceleration mode approximation')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title('Response at node 16 in Y-direction')
plt.legend()
plt.grid(True)
plt.show()

""" Plotting the PSD over the frequency """
plt.plot(2*np.pi*freq, psd, label='PSD(w)')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('PSD (m²/s⁴)')
plt.title('Power Spectral Density of the acceleration at node 16 in Y-direction')
plt.legend()
plt.grid(True)
plt.show()