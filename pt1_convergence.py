import numpy as np
import matplotlib.pyplot as plt

numbers_py = [i for i in range(1, 11)]
numbers_NX = [1, 2, 3, 4, 5, 10]
file_names_NX = [f"NX/FrequencyGraphs/{n}elemPerBar.csv" for n in numbers_NX]
file_names_py = [f"python/convergencePt1/freq{n}EPB.csv" for n in numbers_py]


def read_data(file_name):
    data = np.loadtxt(file_name, delimiter=',', skiprows=1)
    freq_indices = data[:, 0]
    freq_values = data[:, 1]
    return freq_indices, freq_values


def create_f_array(file_names, numbers):
    # Initialisation des tableaux pour chaque fréquence
    f = np.zeros((6, len(numbers)))
    for i, file_name in enumerate(file_names):
        _, freq_values = read_data(file_name)
        for j in range(6):
            f[j, i] = freq_values[j]
    return f

    # Tracé d'un plot par fréquence


f_py = create_f_array(file_names_py, numbers_py)
f_NX = create_f_array(file_names_NX, numbers_NX)

for j in range(6):

    plt.rcParams.update({
        "font.family": "DejaVu Serif",
        # ou ["DejaVu Serif"], ["STIX"], etc.
        "font.serif": ["Arial"],
        "font.size": 20,
        "axes.labelsize": 20,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
    })

    fig = plt.figure("Figure X", figsize=(9, 6))

    # axes titles
    plt.xlabel('Number of elements per beam (N)')
    plt.ylabel('Frequency (Hz)')
    plt.plot(numbers_py, f_py[j, :], 'o-', label=f'python', linewidth=1)
    plt.plot(numbers_NX, f_NX[j, :], 's--', label=f'NX', linewidth=1)

    plt.locator_params(axis='y', nbins=10)
    plt.xticks([0, 2, 4, 6, 8, 10])

    plt.tick_params(axis='both', which='major',
                    labelsize=16, width=1.5, length=5)
    plt.tick_params(axis='both', which='minor',
                    labelsize=12, width=1, length=3)

    # plt.legend(loc='upper right', frameon=False)

    plt.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.5)
    plt.legend()

    plt.tight_layout()

    fig.savefig(f'pt1_convergence/convergencePlotF{j+1}.pdf', dpi=300)
    plt.close()  # Ferme la figure pour éviter la superposition
