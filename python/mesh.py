import numpy as np
from constants import elemPerBar
import matplotlib.pyplot as plt
from geometry import eL, nL, etL, constrainedNodes


def interpolation(node1, node2, alpha):
    return (1 - alpha) * node1 + alpha * node2


def create_elemList_and_nodeList_and_elemTypeList(elemPerBar, elemList=eL, nodeList=nL, elemTypeList=etL):
    if elemPerBar == 1:
        return elemList, nodeList, elemTypeList
    else:
        elemL = np.zeros(
            (elemPerBar*np.shape(elemList)[0], 2), dtype=int)
        nodeL = nodeList.tolist()
        elemTypeL = []
        arrayPositionIndex = 0
        nodeCounter = elemList.max() + 1

        for i in range(len(elemList)):
            node1 = elemList[i, 0]
            node2 = elemList[i, 1]
            coord1 = nodeList[node1 - 1, :]
            coord2 = nodeList[node2 - 1, :]
            elemType = elemTypeList[i]

            for j in range(elemPerBar):
                alpha = (j+1) / elemPerBar
                if j < elemPerBar - 1:
                    newNode = nodeCounter
                    nodeCounter += 1
                    nodeL.append(interpolation(coord1, coord2, alpha))
                else:
                    newNode = node2

                elemL[arrayPositionIndex, :] = [node1, newNode]
                elemTypeL.append(elemType)
                arrayPositionIndex += 1

                node1 = newNode
    return np.array(elemL), np.array(nodeL), np.array(elemTypeL, dtype=int)


def create_dofList(nodeList=nL):
    dofList = np.zeros((np.shape(nodeList)[0], 6), dtype=int)

    for i in range(np.shape(nodeList)[0]):
        dofList[i, :] = [i*6 + 1, i*6+2, i*6+3, i*6 + 4, i*6+5, i*6+6]

    return dofList


def create_locel(dofList, elemList=eL):
    locel = np.zeros((np.shape(elemList)[0], 12), dtype=int)

    for i in range(np.shape(elemList)[0]):
        node1 = elemList[i, 0] - 1
        node2 = elemList[i, 1] - 1

        locel[i, :6] = dofList[node1, :]
        locel[i, 6:] = dofList[node2, :]

    return locel


def set_axes_equal(ax):
    # Equal aspect ratio for 3D plots
    x_lim = ax.get_xlim3d()
    y_lim = ax.get_ylim3d()
    z_lim = ax.get_zlim3d()

    x_range = x_lim[1] - x_lim[0]
    y_range = y_lim[1] - y_lim[0]
    z_range = z_lim[1] - z_lim[0]

    max_range = max([x_range, y_range, z_range]) / 2

    mid_x = 0.5 * (x_lim[0] + x_lim[1])
    mid_y = 0.5 * (y_lim[0] + y_lim[1])
    mid_z = 0.5 * (z_lim[0] + z_lim[1])


slide_params = {
    # Text Sizes (Large for visibility)
    'font.size': 24,           # General default
    'axes.labelsize': 20,      # x and y labels
    'axes.titlesize': 24,      # Title
    'xtick.labelsize': 20,     # Tick numbers
    'ytick.labelsize': 20,
    'legend.fontsize': 18,     # Legend text


    # Structural Geometries
    'axes.linewidth': 2.0,     # Thicker spines (box)
    'xtick.major.width': 2.0,  # Thicker ticks
    'ytick.major.width': 2.0,
    'xtick.major.size': 8.0,   # Longer ticks
    'ytick.major.size': 8.0,

    # Slide Aesthetics
    'font.family': 'sans-serif',  # Sans-serif is more legible on slides than serif
    'figure.autolayout': True,   # Similar to tight_layout
}

# Apply the parameters globally
plt.rcParams.update(slide_params)


def plot_structure(elemList, nodeList, view=None):

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Plot nodes
    ax.scatter(nodeList[:, 0], nodeList[:, 1],
               nodeList[:, 2], c='r', s=10, label='Nodes')

    # Plot elements
    for elem in elemList:
        x = [nodeList[elem[0]-1, 0], nodeList[elem[1]-1, 0]]
        y = [nodeList[elem[0]-1, 1], nodeList[elem[1]-1, 1]]
        z = [nodeList[elem[0]-1, 2], nodeList[elem[1]-1, 2]]
        ax.plot(x, y, z, 'b')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    plt.tight_layout()

    if view == 'XY':
        ax.view_init(elev=90, azim=-90)  # XY
        plt.savefig('plots/mode_6XY_py.png', bbox_inches='tight', pad_inches=0)
    elif view == 'XZ':
        ax.view_init(elev=0, azim=90)  # plan XZ vu de face
        plt.savefig('plots/mode_6XZ_py.png', bbox_inches='tight')

    elif view is None:
        plt.show()


def plot_structureBis(elemList, nodeList, view=None):

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Plot nodes
    ax.scatter(nodeList[:, 0], nodeList[:, 1], nodeList[:, 2],
               c='b', s=20)

    # Plot elements
    for elem in elemList:
        x = [nodeList[elem[0]-1, 0], nodeList[elem[1]-1, 0]]
        y = [nodeList[elem[0]-1, 1], nodeList[elem[1]-1, 1]]
        z = [nodeList[elem[0]-1, 2], nodeList[elem[1]-1, 2]]
        ax.plot(x, y, z, 'grey', lw=2)

    # No grid, clean background
    ax.grid(False)
    ax.set_facecolor((1, 1, 1, 0))      # transparent axis background
    fig.patch.set_facecolor((1, 1, 1, 0))  # transparent figure background

    # Keep axes visible but clean
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_xlim3d(0, 15)
    ax.set_ylim3d(0, 4)
    ax.set_zlim3d(0, 1)
    ax.set_xticks([0, 15])
    ax.set_yticks([0, 4])
    ax.set_zticks([0, 1])

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.tight_layout()

    # Views
    if view == 'XY':
        ax.view_init(elev=90, azim=-90)
        plt.savefig('plots/mode_6XY_py.png', bbox_inches='tight',
                    pad_inches=0, dpi=300, transparent=True)

    elif view == 'XZ':
        ax.view_init(elev=0, azim=90)
        plt.savefig('plots/mode_6XZ_py.png', bbox_inches='tight',
                    pad_inches=0, dpi=300, transparent=True)

    elif view is None:
        plt.savefig("presentation/pt1_model.svg",
                    transparent=True)
        plt.show()


# 1) create elemList, elemTypeList and nodeList with create_elemList_and_nodeList_and_elemTypeList
# 2) create dofList with create_dofList
# 3) create locel with create_locel

elemList, nodeList, elemTypeList = create_elemList_and_nodeList_and_elemTypeList(
    elemPerBar)
dofList = create_dofList(nodeList)
locel = create_locel(dofList, elemList)


# plot_structure(elemList, nodeList)


def vector_to_constrained(vec):
    constrainedDOFs = []
    for i in range(len(constrainedNodes)):
        node = constrainedNodes[i]
        constrainedDOFs.extend(dofList[node-1, :] - 1)
    freeDOFs = np.array([i for i in range(vec.shape[0])
                        if i not in constrainedDOFs])

    vec = vec[freeDOFs].reshape(-1)
    return vec
