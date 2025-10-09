import numpy as np
from constants import elemPerBar
import matplotlib.pyplot as plt
from geometry import eL, nL, etL


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


def plot_structure(elemList, nodeList):

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Plot nodes
    ax.scatter(nodeList[:, 0], nodeList[:, 1],
               nodeList[:, 2], c='r', s=5, label='Nodes')

    # Plot elements
    for elem in elemList:
        x = [nodeList[elem[0]-1, 0], nodeList[elem[1]-1, 0]]
        y = [nodeList[elem[0]-1, 1], nodeList[elem[1]-1, 1]]
        z = [nodeList[elem[0]-1, 2], nodeList[elem[1]-1, 2]]
        ax.plot(x, y, z, 'b')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Frame Structure')
    ax.legend()
    plt.show()


# 1) create elemList, elemTypeList and nodeList with create_elemList_and_nodeList_and_elemTypeList
# 2) create dofList with create_dofList
# 3) create locel with create_locel

elemList, nodeList, elemTypeList = create_elemList_and_nodeList_and_elemTypeList(
    elemPerBar)
dofList = create_dofList(nodeList)
locel = create_locel(dofList, elemList)


# plot_structure(elemList, nodeList)
