import numpy as np

eL = np.array([[1, 2],  # elemList with initial geometry of 1 elem per bar
               [2, 3],
               [3, 4],
               [4, 5],
               [5, 6],
               [12, 13],
               [13, 14],
               [14, 15],
               [15, 16],
               [16, 17],  # Elem10
               [1, 7],
               [11, 6],
               [12, 18],
               [22, 17],
               [7, 8],
               [8, 9],
               [9, 10],
               [10, 11],
               [18, 19],
               [19, 20],  # Elem20
               [20, 21],
               [21, 22],
               [7, 2],
               [2, 8],
               [8, 3],
               [3, 9],
               [9, 4],
               [4, 10],
               [10, 5],
               [5, 11],  # Elem30
               [18, 13],
               [13, 19],
               [19, 14],
               [14, 20],
               [20, 15],
               [15, 21],
               [21, 16],
               [16, 22],
               [7, 18],
               [2, 13],  # Elem 40
               [8, 19],
               [3, 14],
               [9, 20],
               [4, 15],
               [10, 21],
               [5, 16],
               [11, 22]])

nL = np.array([[0, 4, 1],   # nodeList with initial geometry of 1 elem per bar
              [3, 4, 1],
              [6, 4, 1],
              [9, 4, 1],
              [12, 4, 1],
              [15, 4, 1],
              [1.5, 4, 0],
              [4.5, 4, 0],
              [7.5, 4, 0],
              [10/5, 4, 0],
              [13.5, 4, 0],
              [0, 0, 1],
              [3, 0, 1],
              [6, 0, 1],
              [9, 0, 1],
              [12, 0, 1],
              [15, 0, 1],
              [1.5, 0, 0],
              [4.5, 0, 0],
              [7.5, 0, 0],
              [10.5, 0, 0],
              [13.5, 0, 0]])


def interpolation(node1, node2, alpha):
    return (1 - alpha) * node1 + alpha * node2


def create_elemList_and_nodeList(elemPerBar, elemList=eL, nodeList=nL):
    if elemPerBar == 1:
        return elemList, nodeList
    else:
        elemL = np.zeros(
            (elemPerBar*np.shape(elemList)[0], 2), dtype=int)
        nodeL = nodeList.tolist()
        arrayPositionIndex = 0
        nodeCounter = elemList.max() + 1

        for i in range(len(elemList)):
            node1 = elemList[i, 0]
            node2 = elemList[i, 1]
            coord1 = nodeList[node1 - 1, :]
            coord2 = nodeList[node2 - 1, :]

            for j in range(elemPerBar):
                alpha = (j+1) / elemPerBar
                if j < elemPerBar - 1:
                    newNode = nodeCounter
                    nodeCounter += 1
                    nodeL.append(interpolation(coord1, coord2, alpha))
                else:
                    newNode = node2

                elemL[arrayPositionIndex, :] = [node1, newNode]
                arrayPositionIndex += 1

                node1 = newNode
    return np.array(elemL), np.array(nodeL)


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


# 1) create elemList and nodeList with create_elemList_and_nodeList
# 2) create dofList with create_dofList
# 3) create locel with create_locel


elemPerBar = 3
elemList, nodeList = create_elemList_and_nodeList(elemPerBar)
dofList = create_dofList(nodeList)
locel = create_locel(dofList, elemList)
