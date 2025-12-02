import numpy as np
import matplotlib.pyplot as plt


slide_params = {
    # Text Sizes (Large for visibility)
    'font.size': 20,           # General default
    'axes.labelsize': 22,      # x and y labels
    'axes.titlesize': 24,      # Title
    'xtick.labelsize': 20,     # Tick numbers
    'ytick.labelsize': 20,
    'legend.fontsize': 18,     # Legend text

    # Line & Marker Geometries (Thick for projectors)
    'lines.linewidth': 3.5,    # Thicker data lines
    'lines.markersize': 10,    # Much larger markers (default was too small)
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


def plotSlides(x, y,  marker='o', xLabel=None, yLabel=None, label=None, yScaleLog=False, savePath=None, hline=None, vline=None):  # one line to plot

    # renaming
    labeling = label

    fig, ax = plt.subplots()  # Note: figsize is now handled by params
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.plot(x, y, marker=marker, label=labeling)

    if xLabel is not None:
        ax.set_xlabel(xLabel)   # Fontsize is auto-set by params
    if yLabel is not None:
        ax.set_ylabel(yLabel)        # Fontsize is auto-set by params

    if yScaleLog is True:
        ax.set_yscale('log')

    if hline is not None:
        ax.hlines(hline, xmin=np.min(x), xmax=np.max(
            x), color='red', linestyles='--')

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
              ncol=4, frameon=False)
    plt.savefig(savePath + '.svg', transparent=True)


def plotSlidesMulti(x, y, marker='o', xLabel=None, yLabel=None, labels=None, yScaleLog=False, savePath=None, hline=None, vline=None):  # multiple lines to plot

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

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
              ncol=4, frameon=False)
    plt.savefig(savePath + '.svg', transparent=True)
