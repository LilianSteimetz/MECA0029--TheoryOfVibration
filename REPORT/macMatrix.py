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


def plotMacMatrix(mac_matrix, title="MAC Matrix", show=True, save_path=None, red=None):
    plt.rcParams.update({
        "font.family": "DejaVu Serif",
        # ou ["DejaVu Serif"], ["STIX"], etc.
        "font.serif": ["Times New Roman"],
        "font.size": 16,
        "axes.labelsize": 16,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
    })

    fig, ax = plt.subplots(figsize=(6, 6))
    cax = ax.matshow(mac_matrix, cmap='Blues', vmin=0, vmax=1)

    n_rows, n_cols = mac_matrix.shape

    # Tick positions
    ax.set_xticks(np.arange(n_cols))
    ax.set_yticks(np.arange(n_rows))

    # Tick labels: 1, 2, 3, ...
    ax.set_xticklabels(np.arange(1, n_cols+1))
    ax.set_yticklabels(np.arange(1, n_rows+1))

    # Axis names
    if red == 'GI':
        ax.set_xlabel("Modes, Guyan-Irons Reduced Model")
        ax.set_ylabel("Modes, Full Model")
    elif red == 'CB':
        ax.set_xlabel("Modes, Craig-Bampton Reduced Model")
        ax.set_ylabel("Modes, Full Model")
    else:
        ax.set_xlabel("Modes, Reduced Model")
        ax.set_ylabel("Modes, Full Model")

    # Colorbar
    cbar = fig.colorbar(cax, shrink=0.6)
    cbar.set_label('Correlation')

    # Annotate values
    for (i, j), value in np.ndenumerate(mac_matrix):
        ax.text(j, i, f"{value:.2f}", ha='center', va='center',
                color='white' if value > 0.5 else 'black', fontsize=12)

    ax.set_title(title)
    fig.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()
