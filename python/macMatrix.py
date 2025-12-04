import numpy as np
import matplotlib.pyplot as plt


def macMatrix(eigvecsMat1, eigvecsMat2):
    MAC = np.zeros((eigvecsMat1.shape[1], eigvecsMat2.shape[1]))
    for i in range(eigvecsMat1.shape[1]):
        den_i = eigvecsMat1[:, i].T @ eigvecsMat1[:, i]
        for j in range(eigvecsMat2.shape[1]):

            num = (eigvecsMat1[:, i].T @ eigvecsMat2[:, j])**2
            den_j = eigvecsMat2[:, j].T @ eigvecsMat2[:, j]
            den = den_i * den_j

            MAC[i, j] = num/den
    return MAC


slide_params = {
    # Text Sizes (Large for visibility)
    'font.size': 20,           # General default
    'axes.labelsize': 24,      # x and y labels
    'axes.titlesize': 24,      # Title
    'xtick.labelsize': 20,     # Tick numbers
    'ytick.labelsize': 20,
    'legend.fontsize': 18,     # Legend text

    # Line & Marker Geometries (Thick for projectors)
    'lines.linewidth': 3.5,    # Thicker data lines
    'lines.markersize': 10,    # Much larger markers (default was too small)
    'lines.markeredgewidth': 0,  # Remove marker outline for cleaner look

    # Structural Geometries
    'axes.linewidth': 1.0,     # Thicker spines (box)
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


def plotMacMatrix(mac_matrix, title="MAC Matrix", show=True, save_path=None, red=None):

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')
    cax = ax.matshow(mac_matrix, cmap='Greys', vmin=0, vmax=1)

    n_rows, n_cols = mac_matrix.shape

    # Tick positions
    ax.set_xticks(np.arange(n_cols))
    ax.set_yticks(np.arange(n_rows))

    # Tick labels: 1, 2, 3, ...
    ax.set_xticklabels(np.arange(1, n_cols+1))
    ax.set_yticklabels(np.arange(1, n_rows+1))

    # Axis names
    if red == 'GI':
        ax.set_xlabel("Guyan-Irons")
        ax.set_ylabel("Full Model")
    elif red == 'CB':
        ax.set_xlabel("Craig-Bampton")
        ax.set_ylabel("Full model")
    else:
        ax.set_xlabel("Modes, Reduced Model")
        ax.set_ylabel("Modes, Full Model")

    # Colorbar
    cbar = fig.colorbar(cax, shrink=0.6)
    cbar.set_label('Correlation')

    # Annotate values
    for (i, j), value in np.ndenumerate(mac_matrix):

        if abs(value) < 1e-2:
            txt = "0"          # or "" to hide
        else:
            txt = f"{value:.2f}"

        ax.text(
            j, i, txt,
            ha='center', va='center',
            color='white' if value > 0.5 else 'black',
            fontsize=14      # larger text now fits
        )

    ax.set_title(title)
    fig.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()
