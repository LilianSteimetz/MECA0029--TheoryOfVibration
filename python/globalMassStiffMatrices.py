from constants import *
from elemMassMatrices import M_frame_hor, M_frame_diag, M_support_diag, M_support_trans
from elemStiffnessMatrices import K_frame_hor, K_frame_diag, K_support_diag, K_support_trans
from mesh import elemList, nodeList, elemTypeList, dofList, locel
from geometry import constrainedNodes, lumpedNodes
import numpy as np


def create_rotation_matrix_12x12(elemList, nodeList, l, i):
    # local system of coord
    e_zRef = np.array([0, 0, 1])
    e_xLoc = 1/l * np.array([nodeList[elemList[i, 1]-1, 0] - nodeList[elemList[i, 0]-1, 0],
                             nodeList[elemList[i, 1]-1, 1] -
                             nodeList[elemList[i, 0]-1, 1],
                             nodeList[elemList[i, 1]-1, 2] - nodeList[elemList[i, 0]-1, 2]])

    if abs(np.dot(e_xLoc, e_zRef)) > 0.99:
        # choose another reference (global x) when nearly parallel
        e_ref = np.array([1.0, 1.0, 1.0])
    else:
        e_ref = e_zRef
    e_xLoc /= np.linalg.norm(e_xLoc)
    e_yLoc = np.cross(e_zRef, e_xLoc)
    e_yLoc = 1/np.linalg.norm(e_yLoc) * e_yLoc
    e_zLoc = np.cross(e_xLoc, e_yLoc)

    # 3x3 rotation matrix
    R_3x3 = np.zeros((3, 3))
    R_3x3[:, 0] = e_xLoc
    R_3x3[:, 1] = e_yLoc
    R_3x3[:, 2] = e_zLoc

    # 12x12 rotation matrix
    R_6x6 = np.zeros((6, 6))
    R_6x6[:3, :3] = R_3x3
    R_6x6[3:, 3:] = R_3x3
    R_12x12 = np.zeros((12, 12))
    R_12x12[:6, :6] = R_6x6
    R_12x12[6:, 6:] = R_6x6

    return R_12x12


def addLumpedMasses(M, lumpedNodes, lumpedMass):
    for i in range(len(lumpedNodes)):
        node = lumpedNodes[i]
        for i in range(3):
            dofIdx = 6 * (node - 1) + i
            M[dofIdx, dofIdx] += lumpedMass
    return M


def constraintMatReduction(K, M, constrainedDOFs):
    freeDOFs = np.array([i for i in range(K.shape[0])
                        if i not in constrainedDOFs])

    K_reduced = K[np.ix_(freeDOFs, freeDOFs)]
    M_reduced = M[np.ix_(freeDOFs, freeDOFs)]

    return K_reduced, M_reduced


def create_globalMass_and_globalStiffness(elemList=elemList, elemTypeList=elemTypeList, dofList=dofList, locel=locel, constrainedNodes=constrainedNodes, lumpedNodes=lumpedNodes, lumpedMass=lumpedMass):
    nDOF = np.max(dofList)
    M_global = np.zeros((nDOF, nDOF))
    K_global = np.zeros((nDOF, nDOF))
    # reference axis for the local axes construction,
    # taken as e_z global because no vertical element

    for i in range(np.shape(elemList)[0]):
        elemType = elemTypeList[i]
        if elemType == 1:
            M_elem = M_frame_hor
            K_elem = K_frame_hor
            l = l_horizontal
        elif elemType == 2:
            M_elem = M_frame_diag
            K_elem = K_frame_diag
            l = l_diagonal
        elif elemType == 3:
            M_elem = M_support_diag
            K_elem = K_support_diag
            l = l_diagonal
        elif elemType == 4:
            M_elem = M_support_trans
            K_elem = K_support_trans
            l = l_transverse

        # Rotation matrix for the current element
        R_12x12 = create_rotation_matrix_12x12(
            elemList, nodeList, l, i)

        dofIndices = locel[i, :] - 1
        M_es = R_12x12.T @ M_elem @ R_12x12
        K_es = R_12x12.T @ K_elem @ R_12x12

        """
        for a in range(12):
            for b in range(12):
                M_global[dofIndices[a], dofIndices[b]] += M_es[a, b]
                K_global[dofIndices[a], dofIndices[b]] += K_es[a, b]
        """
        # More efficient way to do the same as above
        M_global[np.ix_(dofIndices, dofIndices)] += M_es
        K_global[np.ix_(dofIndices, dofIndices)] += K_es

    # Add lumped masses
    M_global = addLumpedMasses(M_global, lumpedNodes, lumpedMass)

    # Reduce matrices by applying clamped boundary conditions
    constrainedDOFs = []
    for i in range(len(constrainedNodes)):
        node = constrainedNodes[i]
        constrainedDOFs.extend(dofList[node-1, :] - 1)

    K_global, M_global = constraintMatReduction(
        K_global, M_global, constrainedDOFs)

    return M_global, K_global
