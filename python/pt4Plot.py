import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("pt4/CB_errors.csv", delimiter=',', skiprows=1)
internal_modes = data[:, 0]
reduction_times = data[:, 1]
eigenproblem_times = data[:, 2]
errors = data[:, 3:9]

plt.plot(internal_modes, reduction_times, marker='o', label='Reduction time')
plt.title('reduction time')
plt.show()

plt.plot(internal_modes, eigenproblem_times,
         marker='o', label='Eigenproblem time')
plt.title('eigenproblem time')
plt.show()

for i in range(6):
    plt.plot(internal_modes, errors[:, i],
             marker='o', markersize=3,  label=f'Mode {i+1}')
plt.yscale('log')
plt.xlabel('Number of retained internal modes')
plt.ylabel('Relative error (%)')
plt.hlines(1, colors='r', linestyles='dashed',
           label='1% error', xmin=1, xmax=internal_modes[-1])
plt.grid(True, ls="--", alpha=0.4)
plt.show()
