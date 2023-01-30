import numpy as np

import functions


class Node:
    BC = bool()

    def __init__(self, x, y):
        self.x = x
        self.y = y


class Element:
    ID = []  # up to 4!
    vectorP = [.0, .0, .0, .0]

    def __init__(self, ID):
        self.ID = ID


class GlobalData:
    simTime = int()
    simStepTime = int()
    conductivity = int()
    alpha = int()
    tot = int()
    initialTemp = int()
    density = int()
    specificHeat = int()


class Grid:
    nodesNumber = int()
    elementsNumber = int()
    nodes = []
    elements = []


class SC:
    intPt2_ksi = [-1 / np.sqrt(3), 1 / np.sqrt(3), -1 / np.sqrt(3), 1 / np.sqrt(3)]  # for n = 2
    intPt2_eta = [-1 / np.sqrt(3), -1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)]
    ptWeight2 = [1, 1]  # weights

    intPt3_ksi = [-np.sqrt(3 / 5), 0, np.sqrt(3 / 5), -np.sqrt(3 / 5), 0,
              np.sqrt(3 / 5), -np.sqrt(3 / 5), 0, np.sqrt(3 / 5)]  # points for n = 3
    intPt3_eta = [-np.sqrt(3 / 5), -np.sqrt(3 / 5), -np.sqrt(3 / 5), 0, 0,
                  0, np.sqrt(3 / 5), np.sqrt(3 / 5), np.sqrt(3 / 5)]  # points for n = 3

    ptWeight3 = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]  # weights

    intPt4_ksi = [-np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)), -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
                  np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)), np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)),
                  -np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)), -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
                  np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)), np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)),
                  -np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)), -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
                  np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)), np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)),
                  -np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)), -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
                  np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)), np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5))]

    intPt4_eta = [-np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)), -np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)),
                  -np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)), -np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)),
                  -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)), -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
                  -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)), -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
                  np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)), np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
                  np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)), np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
                  np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)), np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)),
                  np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)), np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5))]

    ptWeight4 = [(18 - np.sqrt(30)) / 36, (18 + np.sqrt(30)) / 36, (18 + np.sqrt(30)) / 36, (18 - np.sqrt(30)) / 36]

class El4:
    ksi = []
    eta = []
    N = []

    @staticmethod
    def fill(number):
        El4.ksi.clear()
        El4.eta.clear()
        El4.N.clear()

        if number == 2:
            for i in range(number ** 2):
                El4.ksi.append([])
                El4.eta.append([])
                for j in range(4):
                    El4.ksi[i].append(functions.Nksi(j, SC.intPt2_eta[i]))
                    El4.eta[i].append(functions.Neta(j, SC.intPt2_ksi[i]))
            for i in range(number ** 2):
                El4.N.append(functions.n1(SC.intPt2_ksi[i], SC.intPt2_eta[i]))
                El4.N.append(functions.n2(SC.intPt2_ksi[i], SC.intPt2_eta[i]))
                El4.N.append(functions.n3(SC.intPt2_ksi[i], SC.intPt2_eta[i]))
                El4.N.append(functions.n4(SC.intPt2_ksi[i], SC.intPt2_eta[i]))
            El4.N = np.array(El4.N).reshape(4, 4).tolist()
        if number == 3:
            for i in range(number ** 2):
                El4.ksi.append([])
                El4.eta.append([])
                for j in range(4):
                    El4.ksi[i].append(functions.Nksi(j, SC.intPt3_eta[i]))
                    El4.eta[i].append(functions.Neta(j, SC.intPt3_ksi[i]))
            for i in range(number ** 2):
                El4.N.append(functions.n1(SC.intPt3_ksi[i], SC.intPt3_eta[i]))
                El4.N.append(functions.n2(SC.intPt3_ksi[i], SC.intPt3_eta[i]))
                El4.N.append(functions.n3(SC.intPt3_ksi[i], SC.intPt3_eta[i]))
                El4.N.append(functions.n4(SC.intPt3_ksi[i], SC.intPt3_eta[i]))
            El4.N = np.array(El4.N).reshape(9, 4).tolist()

        if number == 4:
            for i in range(number ** 2):
                El4.ksi.append([])
                El4.eta.append([])
                for j in range(4):
                    El4.ksi[i].append(functions.Nksi(j, SC.intPt4_eta[i]))
                    El4.eta[i].append(functions.Neta(j, SC.intPt4_ksi[i]))
            for i in range(number ** 2):
                El4.N.append(functions.n1(SC.intPt4_ksi[i], SC.intPt4_eta[i]))
                El4.N.append(functions.n2(SC.intPt4_ksi[i], SC.intPt4_eta[i]))
                El4.N.append(functions.n3(SC.intPt4_ksi[i], SC.intPt4_eta[i]))
                El4.N.append(functions.n4(SC.intPt4_ksi[i], SC.intPt4_eta[i]))
            El4.N = np.array(El4.N).reshape(16, 4).tolist()


class MatrixH:
    H = []
    C = []
    matC = []
    tr = []
    dxdksi = []
    dxdeta = []
    dydksi = []
    dydeta = []
    revJacobian = []
    matJacobian = []
    nx = []
    ny = []
    GH = []
    GC = []
    GP = []
    T = []
    matH = []

    globalP = []

    # coordinates = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    coordinates = [[0.0, 0.0], [0.025, 0.0], [0.025, 0.025], [0.0, 0.025]]


class Side:
    intPt2_ksi = [-1 / np.sqrt(3), 1 / np.sqrt(3), 1, 1, -1 / np.sqrt(3), 1 / np.sqrt(3), -1, -1]
    intPt2_eta = [-1, -1, -1 / np.sqrt(3), 1 / np.sqrt(3), 1, 1, 1 / np.sqrt(3), -1 / np.sqrt(3)]

    intPt3_ksi = [-np.sqrt(3 / 5), 0, np.sqrt(3 / 5), 1, 1, 1, -np.sqrt(3 / 5), 0, np.sqrt(3 / 5), -1, -1, -1]
    intPt3_eta = [-1, -1, -1, -np.sqrt(3 / 5), 0, np.sqrt(3 / 5), 1, 1, 1, -np.sqrt(3 / 5), 0, np.sqrt(3 / 5)]

    intPt4_ksi = [-np.sqrt((3.0 / 7.0) + (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  -np.sqrt((3.0 / 7.0) - (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  np.sqrt((3.0 / 7.0) - (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  np.sqrt((3.0 / 7.0) + (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  1, 1, 1, 1, np.sqrt((3.0 / 7.0) + (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  np.sqrt((3.0 / 7.0) - (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  -np.sqrt((3.0 / 7.0) - (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  -np.sqrt((3.0 / 7.0) + (2.0 / 7.0) * np.sqrt(6.0 / 5.0)), -1, -1, -1, -1]

    intPt4_eta = [-1, -1, -1, -1, -np.sqrt((3.0 / 7.0) + (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  -np.sqrt((3.0 / 7.0) - (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  np.sqrt((3.0 / 7.0) - (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  np.sqrt((3.0 / 7.0) + (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  1, 1, 1, 1, np.sqrt((3.0 / 7.0) + (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  np.sqrt((3.0 / 7.0) - (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  -np.sqrt((3.0 / 7.0) - (2.0 / 7.0) * np.sqrt(6.0 / 5.0)),
                  -np.sqrt((3.0 / 7.0) + (2.0 / 7.0) * np.sqrt(6.0 / 5.0))]

    resultsA = []
    resultsB = []
    resultsC = []
    resultsD = []
    Hbc = []
    vecP = []
    bc = []

    @staticmethod
    def fill(number):
        Side.resultsA.clear()
        Side.resultsB.clear()
        Side.resultsC.clear()
        Side.resultsD.clear()

        if number == 2:
            for i in range(0, 8, 2):
                Side.resultsA.append(functions.n1(Side.intPt2_ksi[i], Side.intPt2_eta[i]))
                Side.resultsA.append(functions.n2(Side.intPt2_ksi[i], Side.intPt2_eta[i]))
                Side.resultsA.append(functions.n3(Side.intPt2_ksi[i], Side.intPt2_eta[i]))
                Side.resultsA.append(functions.n4(Side.intPt2_ksi[i], Side.intPt2_eta[i]))
            for i in range(1, 8, 2):
                Side.resultsB.append(functions.n1(Side.intPt2_ksi[i], Side.intPt2_eta[i]))
                Side.resultsB.append(functions.n2(Side.intPt2_ksi[i], Side.intPt2_eta[i]))
                Side.resultsB.append(functions.n3(Side.intPt2_ksi[i], Side.intPt2_eta[i]))
                Side.resultsB.append(functions.n4(Side.intPt2_ksi[i], Side.intPt2_eta[i]))
            Side.resultsA = np.array(Side.resultsA).reshape(4, 4).tolist()
            Side.resultsB = np.array(Side.resultsB).reshape(4, 4).tolist()

        if number == 3:
            for i in range(0, 12, 3):
                Side.resultsA.append(functions.n1(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
                Side.resultsA.append(functions.n2(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
                Side.resultsA.append(functions.n3(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
                Side.resultsA.append(functions.n4(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
            for i in range(1, 12, 3):
                Side.resultsB.append(functions.n1(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
                Side.resultsB.append(functions.n2(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
                Side.resultsB.append(functions.n3(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
                Side.resultsB.append(functions.n4(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
            for i in range(2, 12, 3):
                Side.resultsC.append(functions.n1(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
                Side.resultsC.append(functions.n2(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
                Side.resultsC.append(functions.n3(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
                Side.resultsC.append(functions.n4(Side.intPt3_ksi[i], Side.intPt3_eta[i]))
            Side.resultsA = np.array(Side.resultsA).reshape(4, 4).tolist()
            Side.resultsB = np.array(Side.resultsB).reshape(4, 4).tolist()
            Side.resultsC = np.array(Side.resultsC).reshape(4, 4).tolist()

        if number == 4:
            for i in range(0, 16, 4):
                Side.resultsA.append(functions.n1(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsA.append(functions.n2(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsA.append(functions.n3(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsA.append(functions.n4(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
            for i in range(1, 16, 4):
                Side.resultsB.append(functions.n1(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsB.append(functions.n2(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsB.append(functions.n3(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsB.append(functions.n4(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
            for i in range(2, 16, 4):
                Side.resultsC.append(functions.n1(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsC.append(functions.n2(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsC.append(functions.n3(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsC.append(functions.n4(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
            for i in range(3, 16, 4):
                Side.resultsD.append(functions.n1(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsD.append(functions.n2(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsD.append(functions.n3(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
                Side.resultsD.append(functions.n4(Side.intPt4_ksi[i], Side.intPt4_eta[i]))
            Side.resultsA = np.array(Side.resultsA).reshape(4, 4).tolist()
            Side.resultsB = np.array(Side.resultsB).reshape(4, 4).tolist()
            Side.resultsC = np.array(Side.resultsC).reshape(4, 4).tolist()
            Side.resultsD = np.array(Side.resultsD).reshape(4, 4).tolist()
            # print(np.array(Side.resultsA))
            # print(np.array(Side.resultsB))
            # print(np.array(Side.resultsC))
            # print(np.array(Side.resultsD))

    @staticmethod
    def multiplication(number, node1, node2):
        tempArr = []
        Side.vecP.clear()
        Side.Hbc.clear()
        tempArr.clear()

        detJ = np.sqrt((node2.y - node1.y) * (node2.y - node1.y) + (node2.x - node1.x) * (node2.x - node1.x)) / 2

        #
        # Matrix Hbc
        #

        if number == 2:
            for k in range(number ** 2):
                tempArr.append([])
                for i in range(4):
                    for j in range(4):
                        temp = Side.resultsA[k][j] * Side.resultsA[k][i] + Side.resultsB[k][j] * Side.resultsB[k][i]
                        temp = temp * (SC.ptWeight2[k % 2] * SC.ptWeight2[k % 2])
                        tempArr[k].append(GlobalData.alpha * temp * detJ)
            for i in range(number ** 2):
                tempArr2 = np.array(tempArr[i]).reshape(4, 4).tolist()
                Side.Hbc.append(tempArr2)
        # print("\n", np.array(Side.Hbc))

        if number == 3:
            for k in range(4):
                tempArr.append([])
                for i in range(4):
                    for j in range(4):
                        temp = SC.ptWeight3[0] * Side.resultsA[k][j] * Side.resultsA[k][i] + \
                               SC.ptWeight3[1] * Side.resultsB[k][j] * Side.resultsB[k][i] + \
                               SC.ptWeight3[2] * Side.resultsC[k][j] * Side.resultsC[k][i]

                        tempArr[k].append(GlobalData.alpha * temp * detJ)
            for i in range(4):
                tempArr2 = np.array(tempArr[i]).reshape(4, 4).tolist()
                Side.Hbc.append(tempArr2)
        # print(np.array(Side.Hbc))

        if number == 4:
            for k in range(4):
                tempArr.append([])
                for i in range(4):
                    for j in range(4):
                        temp = SC.ptWeight4[0] * Side.resultsA[k][j] * Side.resultsA[k][i] + \
                               SC.ptWeight4[1] * Side.resultsB[k][j] * Side.resultsB[k][i] + \
                               SC.ptWeight4[2] * Side.resultsC[k][j] * Side.resultsC[k][i] + \
                               SC.ptWeight4[3] * Side.resultsD[k][j] * Side.resultsD[k][i]

                        tempArr[k].append(GlobalData.alpha * temp * detJ)
            for i in range(4):
                tempArr2 = np.array(tempArr[i]).reshape(4, 4).tolist()
                Side.Hbc.append(tempArr2)
        # print(np.array(Side.Hbc))

        #
        # Vector P
        #

        if number == 2:
            for i in range(number ** 2):
                Side.vecP.append([])
                for j in range(4):
                    Side.vecP[i].append(
                        GlobalData.alpha * (SC.ptWeight2[0] * (Side.resultsA[i][j] * GlobalData.tot) +
                                            SC.ptWeight2[1] * (Side.resultsB[i][j] * GlobalData.tot)) * detJ)
        if number == 3:
            for i in range(4):
                Side.vecP.append([])
                for j in range(4):
                    Side.vecP[i].append(
                        GlobalData.alpha * (SC.ptWeight3[0] * (Side.resultsA[i][j] * GlobalData.tot) +
                                            SC.ptWeight3[1] * (Side.resultsB[i][j] * GlobalData.tot) +
                                            SC.ptWeight3[2] * (Side.resultsC[i][j] * GlobalData.tot)) * detJ)
        if number == 4:
            for i in range(4):
                Side.vecP.append([])
                for j in range(4):
                    Side.vecP[i].append(
                        GlobalData.alpha * (SC.ptWeight4[0] * (Side.resultsA[i][j] * GlobalData.tot) +
                                            SC.ptWeight4[1] * (Side.resultsB[i][j] * GlobalData.tot) +
                                            SC.ptWeight4[2] * (Side.resultsC[i][j] * GlobalData.tot) +
                                            SC.ptWeight4[3] * (Side.resultsD[i][j] * GlobalData.tot)) * detJ)
            # print(np.array(Side.vecP))
        return Side.Hbc
