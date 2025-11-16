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


def plotMacMatrix(mac_matrix, title="MAC Matrix", show=True, save_path=None):

    fig, ax = plt.subplots(figsize=(6, 6))
    cax = ax.matshow(mac_matrix, cmap='Blues', vmin=0, vmax=1)

    # labels for however many modes you keep
    n_full = mac_matrix.shape[0]
    n_red = mac_matrix.shape[1]

    x_labels = [f"Red {j+1}" for j in range(n_red)]
    y_labels = [f"Full {i+1}" for i in range(n_full)]

    ax.set_xticks(np.arange(n_red))
    ax.set_yticks(np.arange(n_full))

    ax.set_xticklabels(x_labels, rotation=45, ha='left')
    ax.set_yticklabels(y_labels)

    # colorbar
    cbar = fig.colorbar(cax, shrink=0.8)
    cbar.set_label('MAC Value')

    # annotate cells
    for (i, j), value in np.ndenumerate(mac_matrix):
        ax.text(j, i, f"{value:.2f}", ha='center', va='center',
                color='white' if value > 0.5 else 'black')

    ax.set_title("MAC Matrix: Full vs Craig–Bampton Modes")

    fig.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()
