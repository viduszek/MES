# Required libraries
import numpy as np

# Required .py files
import classes


#
# Funkcje kształtu
#

def n1(ksi, eta):
    return 0.25 * (1.0 - ksi) * (1.0 - eta)


def n2(ksi, eta):
    return 0.25 * (1.0 + ksi) * (1.0 - eta)


def n3(ksi, eta):
    return 0.25 * (1.0 + ksi) * (1.0 + eta)


def n4(ksi, eta):
    return 0.25 * (1.0 - ksi) * (1.0 + eta)


#
# Pochodne funkcji kształtu 
#

def Nksi(number, x):
    if number == 0:
        return -0.25 * (1 - x)
    if number == 1:
        return 0.25 * (1 - x)
    if number == 2:
        return 0.25 * (x + 1)
    if number == 3:
        return -0.25 * (x + 1)


def Neta(number, x):
    if number == 0:
        return -0.25 * (1 - x)
    if number == 1:
        return -0.25 * (1 + x)
    if number == 2:
        return 0.25 * (1 + x)
    if number == 3:
        return 0.25 * (1 - x)

#
# Funkcje do obliczeń
#


def f(x):
    return x + 2


def f1(x):
    return 2 * x ** 2 + 3 * x - 8


def f2(ksi, eta):
    return -5 * ksi ** 2 * eta + 2 * ksi * eta ** 2 + 10


def f3(x):
    return 3 * x ** 2 - 6 * x + 1


def gauss(n, matH, T):
    temp = 0.0
    temp2 = 0.0

    for j in range(0, n-1, 1):
        for i in range(j+1, n, 1):
            temp = matH[i][j] / matH[j][j]
            for k in range(0, n+1, 1):
                matH[i][k] -= matH[j][k] * temp

    for i in range(n-1, -1, -1):
        temp2 = 0
        for j in range(i+1, n, 1):
            temp2 += matH[i][j] * T[j]
        T[i] = (matH[i][n] - temp2) / matH[i][i]


#
# Tworzenie macierzy A: dxdksi, dxdeta, dydksi, dydeta
#

def createMatrix(number):  #

    # print(np.array(classes.El4.ksi))
    # print(np.array(classes.El4.eta))

    matrix = []
    matrix.clear()
    if (number == 2) | (number == 3) | (number == 4):
        for i in range(number ** 2):
            xksi, xeta, yksi, yeta = 0.0, 0.0, 0.0, 0.0
            for j in range(4):
                yeta = yeta + (classes.El4.eta[i][j] * classes.MatrixH.coordinates[j][1])
                yksi = yksi + (classes.El4.ksi[i][j] * classes.MatrixH.coordinates[j][1]) * (-1)
                xeta = xeta + (classes.El4.eta[i][j] * classes.MatrixH.coordinates[j][0]) * (-1)
                xksi = xksi + (classes.El4.ksi[i][j] * classes.MatrixH.coordinates[j][0])
            classes.MatrixH.dxdksi[i] = xksi
            classes.MatrixH.dydksi[i] = yksi
            classes.MatrixH.dxdeta[i] = xeta
            classes.MatrixH.dydeta[i] = yeta
    matrix.append(classes.MatrixH.dxdksi)
    matrix.append(classes.MatrixH.dydksi)
    matrix.append(classes.MatrixH.dxdeta)
    matrix.append(classes.MatrixH.dydeta)
    return matrix


#
# Wyznaczanie jakobianu; revJacobian = 1/detJ, matJacobian = detJ
#

def jacobian(number):
    classes.MatrixH.revJacobian.clear()
    classes.MatrixH.matJacobian.clear()
    if (number == 2) | (number == 3) | (number == 4):
        for i in range(number ** 2):
            temp = (classes.MatrixH.dydeta[i] * classes.MatrixH.dxdksi[i]) - \
                   (classes.MatrixH.dydksi[i] * classes.MatrixH.dxdeta[i])
            classes.MatrixH.revJacobian.append(1 / temp)
            classes.MatrixH.matJacobian.append(temp)
    else:
        print("Error! Number must be in range [2, 4]. Function: ", jacobian.__name__)
    # print(np.array(classes.MatrixH.matJacobian))
    # print(*['{:.20f}'.format(l) for l in classes.MatrixH.revJacobian])

    # print(np.array(classes.MatrixH.revJacobian))
    return classes.MatrixH.revJacobian, classes.MatrixH.matJacobian


#
# Mnożenie (jakobian * macierz A (dxdksi...))
#

def translate(number):  # 80 0 0 80 // jakobiany
    tr = []
    classes.MatrixH.tr.clear()
    # print(classes.MatrixH.matJacobian)

    if (number == 2) | (number == 3) | (number == 4):
        for i in range(number**2):
            xksi, xeta, yksi, yeta = 0.0, 0.0, 0.0, 0.0
            classes.MatrixH.tr.append([])
            for j in range(4):
                yeta = yeta + (classes.El4.eta[i][j] * classes.MatrixH.coordinates[j][1])
                yksi = yksi + (classes.El4.ksi[i][j] * classes.MatrixH.coordinates[j][1]) * (-1)
                xeta = xeta + (classes.El4.eta[i][j] * classes.MatrixH.coordinates[j][0]) * (-1)
                xksi = xksi + (classes.El4.ksi[i][j] * classes.MatrixH.coordinates[j][0])
            classes.MatrixH.tr[i].append(yeta * classes.MatrixH.revJacobian[i])
            classes.MatrixH.tr[i].append(yksi * classes.MatrixH.revJacobian[i])
            classes.MatrixH.tr[i].append(xeta * classes.MatrixH.revJacobian[i])
            classes.MatrixH.tr[i].append(xksi * classes.MatrixH.revJacobian[i])
    else:
        print("Error! Number must be in range [2, 4]. Function: ", translate.__name__)

    # print(np.array(classes.MatrixH.tr))
    return tr

#
# Przekształcanie jakobianu, macierze dn/dx i dn/dy
#


def transform(number):  # dla punktow calkowania dndx, dndy
    classes.MatrixH.nx.clear()
    classes.MatrixH.ny.clear()
    if (number == 2) | (number == 3) | (number == 4):
        for i in range(number**2):
            classes.MatrixH.nx.append([])
            classes.MatrixH.ny.append([])
            for j in range(4):
                classes.MatrixH.nx[i].append(classes.MatrixH.tr[i][0] * classes.El4.ksi[i][j] +
                                             classes.MatrixH.tr[i][1] * classes.El4.eta[i][j])
                classes.MatrixH.ny[i].append(classes.MatrixH.tr[i][2] * classes.El4.ksi[i][j] +
                                             classes.MatrixH.tr[i][3] * classes.El4.eta[i][j])
    else:
        print("Error! Number must be in range [2, 4]. Function: ", transform.__name__)
    # print("\n", np.array(classes.MatrixH.nx), "\n\n", np.array(classes.MatrixH.ny))
    # print(np.array(classes.MatrixH.nx))
    return classes.MatrixH.nx, classes.MatrixH.ny


#
# Wyznaczanie macierzy H dla poszczególnych punktów całkowania, jest ich bardzo dużo, bo dużo raz się powtarza.
#

def calcMatrixH(number):  # calka kt((nx*nx-1)+(ny*ny-1))
    tempArr = []
    tempArr3 = []
    tempArr3.clear()
    tempArr.clear()
    classes.MatrixH.H.clear()
    classes.MatrixH.C.clear()
    if (number == 2) | (number == 3) | (number == 4):
        for k in range(number**2):
            tempArr.append([])
            tempArr3.append([])
            for i in range(4):
                for j in range(4):
                    temp = (classes.MatrixH.nx[k][j] * classes.MatrixH.nx[k][i] + classes.MatrixH.ny[k][j] *
                            classes.MatrixH.ny[k][i]) * classes.MatrixH.matJacobian[k] * classes.GlobalData.conductivity
                    tempArr[k].append(temp)

                    temp = classes.GlobalData.specificHeat * classes.GlobalData.density * classes.MatrixH.matJacobian[k]
                    temp = temp * classes.El4.N[k][j] * classes.El4.N[k][i]
                    tempArr3[k].append(temp)
        for i in range(number ** 2):
            tempArr2 = np.array(tempArr[i]).reshape(4, 4)
            classes.MatrixH.H.append(tempArr2)
            tempArr4 = np.array(tempArr3[i]).reshape(4, 4)
            classes.MatrixH.C.append(tempArr4)

    # print(np.array(classes.GlobalData.conductivity))
    return classes.MatrixH.H

#
# Ostateczna macierz H dla danego elementu skończonego
#


def H(number):
    matH = []
    matH.clear()
    classes.MatrixH.matC.clear()
    tempp = 0
    if number == 2:
        for j in range(4):
            temp = 0
            temp2 = 0
            for k in range(number**2):
                temp = temp + classes.MatrixH.H[k][j] * classes.SC.ptWeight2[k % 2] * classes.SC.ptWeight2[k % 2]
                temp2 = temp2 + classes.MatrixH.C[k][j] * classes.SC.ptWeight2[k % 2] * classes.SC.ptWeight2[k % 2]
            # print(temp)
            matH.append(temp)
            classes.MatrixH.matC.append(temp2)
    elif number == 3:
        for j in range(4):
            temp = 0
            temp2 = 0
            for k in range(number**2):
                # print(classes.MatrixH.H)
                if k > 2:
                    tempp = 1
                if k > 5:
                    tempp = 2
                temp = temp + classes.MatrixH.H[k][j] * classes.SC.ptWeight3[tempp] * classes.SC.ptWeight3[k % 3]
                temp2 = temp2 + classes.MatrixH.C[k][j] * classes.SC.ptWeight3[tempp] * classes.SC.ptWeight3[k % 3]
            # print(temp)
            matH.append(temp)
            classes.MatrixH.matC.append(temp2)
    elif number == 4:
        for j in range(4):
            temp = 0
            temp2 = 0
            for k in range(number**2):
                # print(classes.MatrixH.H)
                if k > 3:
                    tempp = 1
                if k > 7:
                    tempp = 2
                if k > 11:
                    tempp = 3
                temp = temp + classes.MatrixH.H[k][j] * classes.SC.ptWeight4[tempp] * classes.SC.ptWeight4[k % 4]
                temp2 = temp2 + classes.MatrixH.C[k][j] * classes.SC.ptWeight4[tempp] * classes.SC.ptWeight4[k % 4]
            # print(temp)
            matH.append(temp)
            classes.MatrixH.matC.append(temp2)
    # print(np.array(matH))
    return matH


#
# Funkcja "mes": lokalne macierze H, agregacja
#

def mes(number):

    for i in range(classes.Grid.nodesNumber):
        classes.MatrixH.GH.append([])
        for j in range(classes.Grid.nodesNumber):
            classes.MatrixH.GH[i].append(0)

    for i in range(classes.Grid.nodesNumber):
        classes.MatrixH.matH.append([])
        for j in range(classes.Grid.nodesNumber+1):
            classes.MatrixH.matH[i].append(0)

    for i in range(classes.Grid.nodesNumber):
        classes.MatrixH.GC.append([])
        for j in range(classes.Grid.nodesNumber):
            classes.MatrixH.GC[i].append(0)

    for i in range(classes.Grid.elementsNumber):
        nodesArr = []
        for j in range(4):
            temp = int(classes.Grid.elements[i].ID[j]) - 1
            classes.MatrixH.coordinates[j][0] = classes.Grid.nodes[temp].x
            classes.MatrixH.coordinates[j][1] = classes.Grid.nodes[temp].y
            nodesArr.append(classes.Grid.nodes[temp])

        # print(classes.MatrixH.coordinates)
        pom = [[.0, .0, .0, .0], [.0, .0, .0, .0], [.0, .0, .0, .0], [.0, .0, .0, .0]]

        classes.El4.fill(number)  # number # fills ksi and eta arrays
        createMatrix(number)
        jacobian(number)
        translate(number)
        transform(number)
        calcMatrixH(number)

        if nodesArr[0].BC and nodesArr[1].BC:
            node1 = nodesArr[0]
            node2 = nodesArr[1]
            classes.Side.fill(number)
            classes.Side.multiplication(number, node1, node2)
            pom = np.add(pom, classes.Side.Hbc[0])
            classes.Grid.elements[i].vectorP = np.add(classes.Grid.elements[i].vectorP, (classes.Side.vecP[0]))
        if nodesArr[1].BC and nodesArr[2].BC:
            node1 = nodesArr[1]
            node2 = nodesArr[2]
            classes.Side.fill(number)
            classes.Side.multiplication(number, node1, node2)
            pom = np.add(pom, classes.Side.Hbc[1])
            classes.Grid.elements[i].vectorP = np.add(classes.Grid.elements[i].vectorP, (classes.Side.vecP[1]))
        if nodesArr[2].BC and nodesArr[3].BC:
            node1 = nodesArr[2]
            node2 = nodesArr[3]
            classes.Side.fill(number)
            classes.Side.multiplication(number, node1, node2)
            pom = np.add(pom, classes.Side.Hbc[2])
            classes.Grid.elements[i].vectorP = np.add(classes.Grid.elements[i].vectorP, (classes.Side.vecP[2]))
        if nodesArr[3].BC and nodesArr[0].BC:
            node1 = nodesArr[3]
            node2 = nodesArr[0]
            classes.Side.fill(number)
            classes.Side.multiplication(number, node1, node2)
            pom = np.add(pom, classes.Side.Hbc[3])
            classes.Grid.elements[i].vectorP = np.add(classes.Grid.elements[i].vectorP, (classes.Side.vecP[3]))

        temp = np.add(H(number), pom)
        # temp = pom
        # print(np.array(temp))
        # temp = H(number) + pom

#
# Agregacja
#
        # print(classes.Grid.elements[i].vectorP)
        for j in range(4):
            classes.MatrixH.globalP[int(classes.Grid.elements[i].ID[j]) - 1] += classes.Grid.elements[i].vectorP[j]
            for k in range(4):
                classes.MatrixH.GH[int(classes.Grid.elements[i].ID[j]) - 1][int(classes.Grid.elements[i].ID[k]) - 1] \
                    += temp[j][k]
                classes.MatrixH.GC[int(classes.Grid.elements[i].ID[j]) - 1][int(classes.Grid.elements[i].ID[k]) - 1] \
                    += classes.MatrixH.matC[j][k]

    # for i in range(classes.Grid.nodesNumber):
        # print(*['{:.2f}'.format(i) for i in classes.MatrixH.GC[i]])
        # print(*['{:.2f}'.format(i) for i in classes.MatrixH.GH[i]])
        # pass

    t_size = classes.Grid.nodesNumber + 1

    for i in range(classes.Grid.nodesNumber):
        for j in range(classes.Grid.nodesNumber):
            classes.MatrixH.GH[i][j] = classes.MatrixH.GH[i][j] + \
                                       (classes.MatrixH.GC[i][j] / classes.GlobalData.simStepTime)

    # for l in range(classes.Grid.nodesNumber):
    #     print(*['{:.2f}'.format(l) for l in classes.MatrixH.GH[l]])

    for i in range(classes.GlobalData.simStepTime, classes.GlobalData.simTime+1, classes.GlobalData.simStepTime):
        for j in range(classes.Grid.nodesNumber):
            classes.MatrixH.GP[j] = classes.MatrixH.globalP[j]
            for k in range(classes.Grid.nodesNumber):
                classes.MatrixH.GP[j] = classes.MatrixH.GP[j] + \
                                        classes.MatrixH.GC[j][k] / classes.GlobalData.simStepTime \
                                        * classes.MatrixH.T[k]

        for j in range(classes.Grid.nodesNumber):
            for k in range(classes.Grid.nodesNumber):
                classes.MatrixH.matH[j][k] = classes.MatrixH.GH[j][k]
            classes.MatrixH.matH[j][classes.Grid.nodesNumber] = classes.MatrixH.GP[j]

        gauss(classes.Grid.nodesNumber, classes.MatrixH.matH, classes.MatrixH.T)

        # print(np.array(classes.MatrixH.T))

        t0 = classes.MatrixH.T[0]  # min
        t1 = classes.MatrixH.T[1]  # max

        for j in range(classes.Grid.nodesNumber):
            if classes.MatrixH.T[j] > t1:
                t1 = classes.MatrixH.T[j]
            if classes.MatrixH.T[j] < t0:
                t0 = classes.MatrixH.T[j]

        print("Time: ", i, "\tTmin: ", t0, "\tTmax: ", t1)
        # print("Tmin: ", t0)
        # print("Tmax: ", t1)
        # print("===")

    # print("\n", classes.MatrixH.globalP)
    # print(np.array(classes.MatrixH.matC))
    # print(np.array(classes.El4.N))

