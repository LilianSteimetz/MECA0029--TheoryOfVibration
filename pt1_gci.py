import numpy as np
import matplotlib.pyplot as plt


def read_data(file_name):
    data = np.loadtxt(file_name, delimiter=',', skiprows=1)
    return data[:, 1]  # fréquences


def load_freq_data(file_names, N_array, nmodes):
    f = np.zeros((nmodes, len(N_array)))
    for i, file_name in enumerate(file_names):
        freq_values = read_data(file_name)
        f[:, i] = freq_values[:nmodes]
    return f


def calculate_gci(f, N, nmodes, Fs):
    """Calculates the Grid Convergence Index (GCI)."""
    p = np.zeros((nmodes, len(N)-2))
    GCI = np.zeros((nmodes, len(N)-1))

    for m in range(nmodes):            # boucle modes
        for i in range(len(N)-2):      # boucle maillages
            r21 = N[i+1] / N[i]
            r32 = N[i+2] / N[i+1]
            # formule de p ordre de convergence estimé entre maillages i et i+1
            # Add a small epsilon to prevent division by zero if frequencies are identical
            epsilon = 1e-15
            p[m, i] = np.log(abs((f[m, i+1] - f[m, i] + epsilon) /
                                 (f[m, i+2] - f[m, i+1] + epsilon))) / np.log(r21)

            # GCI entre maillages i et i+1
            GCI[m, i] = Fs * abs(f[m, i] - f[m, i+1]) / (r21**p[m, i] - 1)
        # dernier GCI du mode m
        i = len(N) - 3
        r32 = N[i+2] / N[i+1]
        GCI[m, -1] = Fs * abs(f[m, i+1] - f[m, i+2]) / (r32**p[m, i] - 1)

    return GCI


def calculate_delta_rel(f, N, nmodes):
    """Calculates the relative delta between successive meshes."""
    delta_rel = np.zeros((nmodes, len(N)-1))
    for m in range(nmodes):
        for i in range(len(N)-1):
            # Add a small epsilon to prevent division by zero
            delta_rel[m, i] = abs(f[m, i+1] - f[m, i]) / (f[m, i+1] + 1e-15)
    return delta_rel


def plotSlidesMulti(x, y, marker='o', xLabel=None, yLabel=None, labels=None, yScaleLog=False, savePath=None, hline=None, vline=None, xlim=None):  # multiple lines to plot

    # renaming
    labeling = labels

    fig, ax = plt.subplots()  # Note: figsize is now handled by params

    for i in range(len(y)):
        ax.plot(x, y[i],  marker=marker, label=labeling[i])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if xLabel is not None:
        ax.set_xlabel(xLabel)   # Fontsize is auto-set by params
    if yLabel is not None:
        ax.set_ylabel(yLabel)        # Fontsize is auto-set by params

    if yScaleLog is True:
        ax.set_yscale('log')

    if hline is not None:
        ax.hlines(hline, xmin=np.min(x), xmax=np.max(
            x), color='red', linestyles='--')
    if xlim is not None:
        ax.set_xlim(0, xlim)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
              ncol=3, frameon=False)
    ax.set_xticks([1, 2, 3, 4, 5])

    plt.savefig(savePath + '.svg', transparent=True)


slide_params = {
    # Text Sizes (Large for visibility)
    'font.size': 20,           # General default
    'axes.labelsize': 28,      # x and y labels
    'axes.titlesize': 24,      # Title
    'xtick.labelsize': 20,     # Tick numbers
    'ytick.labelsize': 20,
    'legend.fontsize': 18,     # Legend text

    # Line & Marker Geometries (Thick for projectors)
    'lines.linewidth': 5,    # Thicker data lines
    'lines.markersize': 15,    # Much larger markers (default was too small)
    'lines.markeredgewidth': 0,  # Remove marker outline for cleaner look

    # Structural Geometries
    'axes.linewidth': 2.0,     # Thicker spines (box)
    'xtick.major.width': 2.0,  # Thicker ticks
    'ytick.major.width': 2.0,
    'xtick.major.size': 8.0,   # Longer ticks
    'ytick.major.size': 8.0,

    # Slide Aesthetics
    'font.family': 'sans-serif',  # Sans-serif is more legible on slides than serif
    'figure.autolayout': True,   # Similar to tight_layout
    'figure.figsize': (12, 9),   # 16:9 Aspect Ratio by default
}

# Apply the parameters globally
plt.rcParams.update(slide_params)


# Ensure output directory exists
output_dir = 'pt1_convergence'


nmodes = 6

"""
Fs = 1.25  # Factor of Safety

# --- 1. GCI Calculation ---

N_gci = np.array([1, 2, 4, 10])  # Meshes for GCI calculation

# Define file paths for both Python and NX
file_names_py_gci = [f"python/convergencePt1/freq{n}EPB.csv" for n in N_gci]
file_names_NX_gci = [f"NX/FrequencyGraphs/{n}elemPerBar.csv" for n in N_gci]

# Load data
f_py_gci = load_freq_data(file_names_py_gci, N_gci, nmodes)
f_NX_gci = load_freq_data(file_names_NX_gci, N_gci, nmodes)

# Calculate GCI
GCI_py = calculate_gci(f_py_gci, N_gci, nmodes, Fs)
GCI_NX = calculate_gci(f_NX_gci, N_gci, nmodes, Fs)

# --- Plot 1: GCI for Python ---
plot_x_gci = N_gci[:-1]  # GCI is calculated for N pairs, so [1, 2, 4]

fig = plt.figure(figsize=(9, 6))
for m in range(nmodes):
    plt.plot(plot_x_gci, GCI_py[m, :], 'o-', label=f'Mode {m+1}')

plt.xlabel("Number of elements per beam (N)")
plt.ylabel("GCI")
plt.yscale("log")
plt.xticks(plot_x_gci)
plt.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.8)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(f'pt1_convergence/GCI_Python.pdf', dpi=300)
plt.close()

# --- Plot 2: GCI for NX ---
fig = plt.figure(figsize=(9, 6))
for m in range(nmodes):
    plt.plot(plot_x_gci, GCI_NX[m, :], 's-', label=f'Mode {m+1}')

plt.xlabel("Number of elements per beam (N)")
plt.ylabel("GCI")
plt.yscale("log")
plt.xticks(plot_x_gci)
plt.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.8)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(f'pt1_convergence/GCI_NX.pdf', dpi=300)
plt.close()

"""
# --- 2. Delta Relatif Calculation ---

N_delta = np.array([1, 2, 3, 4, 5, 10])

file_names_py_delta = [
    f"python/convergencePt1/freq{n}EPB.csv" for n in N_delta]
file_names_NX_delta = [
    f"NX/FrequencyGraphs/{n}elemPerBar.csv" for n in N_delta]

# Load data
f_py_delta = load_freq_data(file_names_py_delta, N_delta, nmodes)
f_NX_delta = load_freq_data(file_names_NX_delta, N_delta, nmodes)

# Calculate Delta Relatif
delta_rel_py = calculate_delta_rel(f_py_delta, N_delta, nmodes)
delta_rel_NX = calculate_delta_rel(f_NX_delta, N_delta, nmodes)

# Plot Delta Relatif (One plot for all modes)
plot_x_delta = N_delta[:-1]
labels = [f"Mode {i+1}" for i in range(6)]


plotSlidesMulti(plot_x_delta, delta_rel_py, labels=labels, xLabel="Number of elements per beam",
                yLabel=r"$\delta_{rel}^{N \rightarrow N+1}$", yScaleLog=True, savePath="python/presentation/pt1ConvPython")

"""
# --- Plot 3: Delta Relatif for Python ---
fig = plt.figure(figsize=(9, 6))
for m in range(nmodes):
    plt.plot(plot_x_delta, delta_rel_py[m, :], 'o-', label=f'Mode {m+1}')

plt.xlabel("Number of elements per beam (N)")
plt.ylabel(r"$\delta_{rel}^{N \rightarrow N+1}$")
plt.yscale("log")
plt.ylim(1e-7, 3*1e-2)
plt.xticks(plot_x_delta)
plt.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.8)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(f'pt1_convergence/deltaRelatif_Python.pdf', dpi=300)
plt.close()

# --- Plot 4: Delta Relatif for NX ---
fig = plt.figure(figsize=(9, 6))
for m in range(nmodes):
    plt.plot(plot_x_delta, delta_rel_NX[m, :], 's-', label=f'Mode {m+1}')

plt.xlabel("Number of elements per beam (N)")
plt.ylabel(r"$\delta_{rel}^{N \rightarrow N+1}$")
plt.yscale("log")
plt.ylim(1e-7, 3*1e-2)
plt.xticks(plot_x_delta)
plt.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.8)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(f'pt1_convergence/deltaRelatif_NX.pdf', dpi=300)
plt.close()

# Print delta relatif #
print("Delta Relatif Python:")
print(max(delta_rel_py[:, -2]))

print("\nDelta Relatif NX:")
print(max(delta_rel_NX[:, -2]))

"""
