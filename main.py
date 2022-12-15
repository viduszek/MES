# =============================== Libraries/Packages == #
import numpy as np  #
#
#
# =============================== Classes ============= #


class Node:
    x = float()
    y = float()
    BC = bool()

    def __init__(self, x, y):
        self.x = x;
        self.y = y;


class Element:
    ID = []  # up to 4!
    # ID = [4]

    def __init__(self, ID):
        self.ID = ID


class GlobalData:  # there's actually no need to declare those variables, but you can do it tho
    simTime = int()
    simStepTime = int()
    conductivity = int()
    alpha = int()
    tot = int()
    initialTemp = int()
    density = int()
    specificHeat = int()

    # def __init__(self):
    # print("class: GlobalData initialized.")


class Grid:
    nodesNumber = int()
    elementsNumber = int()
    nodes = []
    elements = []

    # def __init__(self):
    # print("class: Grid initialized.")


class SC:
    intPt2 = [-1 / np.sqrt(3), 1 / np.sqrt(3)]  # for n = 2
    ptWeight2 = [1, 1]  # weights

    intPt3 = [np.sqrt(3 / 5), np.sqrt(3 / 5), np.sqrt(3 / 5)]  # points for n = 3
    ptWeight3 = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]  # weights

    intPt4 = [-np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)), -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
              np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)), np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5))]
    ptWeight4 = [(18 - np.sqrt(30)) / 36, (18 + np.sqrt(30)) / 36, (18 + np.sqrt(30)) / 36, (18 - np.sqrt(30)) / 36]

    # def __init__(self):
    # print("class: SC initialized.")


class El4:
    ksi = []
    eta = []

    ksi2 = [1, -1, -1, 1]  # additional weights for the integral points / ksi version
    ksi3 = [-1, -1, -1, 0, 0, 0, 1, 1, 1]
    ksi4 = [1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1]

    eta2 = [1, 1, 1, 1]  # additional weights for the integral points / eta version
    eta3 = [-1, 0, 1, -1, 0, 1, -1, 0, 1]
    eta4 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    @staticmethod
    def fill(number):
        if number == 2:
            for i in range(number ** 2):
                El4.ksi.append([])
                El4.eta.append([])
                for j in range(4):
                    El4.ksi[i].append(Nksi(j, SC.intPt2[i % number] * El4.ksi2[i]))
                    El4.eta[i].append(Neta(j, SC.intPt2[i % number] * El4.eta2[i]))
        elif number == 3:
            for i in range(number ** 2):
                El4.ksi.append([])
                El4.eta.append([])
                for j in range(4):
                    El4.ksi[i].append(Nksi(j, SC.intPt3[i % number] * El4.ksi3[i]))
                    El4.eta[i].append(Neta(j, SC.intPt3[i % number] * El4.eta3[i]))
            # print(np.array(El4.ksi))
            # print(np.array(El4.eta))
        elif number == 4:
            for i in range(number ** 2):
                El4.ksi.append([])
                El4.eta.append([])
                for j in range(4):
                    El4.ksi[i].append(Nksi(j, SC.intPt4[i % number] * El4.ksi4[i]))
                    El4.eta[i].append(Neta(j, SC.intPt4[i % number] * El4.eta4[i]))


class MatrixH:
    H = []
    dxdksi = []
    dxdeta = []
    dydksi = []
    dydeta = []
    revJacobian = []
    matJacobian = []
    nx = []
    ny = []
    GH = []

    # coordinates = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    coordinates = [[0.0, 0.0], [0.025, 0.0], [0.025, 0.025], [0.0, 0.025]]


def createMatrix(number):  #
    matrix = []
    if number == 2:
        for i in range(number ** 2):
            xksi, xeta, yksi, yeta = 0.0, 0.0, 0.0, 0.0
            for j in range(4):
                xksi = xksi + (El4.ksi[i][j] * MatrixH.coordinates[j][0])
                xeta = xeta + (El4.eta[i][j] * MatrixH.coordinates[j][0])
                yksi = yksi + (El4.ksi[i][j] * MatrixH.coordinates[j][1])
                yeta = yeta + (El4.eta[i][j] * MatrixH.coordinates[j][1])
            MatrixH.dxdksi.append(xksi)
            MatrixH.dxdeta.append(xeta)
            MatrixH.dydksi.append(yksi)
            MatrixH.dydeta.append(yeta)
        matrix = list(zip(MatrixH.dxdksi, MatrixH.dxdeta, MatrixH.dydksi, MatrixH.dydeta))
    elif number == 3:
        for i in range(number ** 2):
            xksi, xeta, yksi, yeta = 0.0, 0.0, 0.0, 0.0
            for j in range(4):
                xksi = xksi + (El4.ksi[i][j] * MatrixH.coordinates[j][0])
                xeta = xeta + (El4.eta[i][j] * MatrixH.coordinates[j][0])
                yksi = yksi + (El4.ksi[i][j] * MatrixH.coordinates[j][1])
                yeta = yeta + (El4.eta[i][j] * MatrixH.coordinates[j][1])
            MatrixH.dxdksi.append(xksi)
            MatrixH.dxdeta.append(xeta)
            MatrixH.dydksi.append(yksi)
            MatrixH.dydeta.append(yeta)
        matrix = list(zip(MatrixH.dxdksi, MatrixH.dxdeta, MatrixH.dydksi, MatrixH.dydeta))

    # print(np.array(matrix))
    return matrix


def jacobian(number):
    if number == 2:
        for i in range(number ** 2):
            temp = (MatrixH.dxdksi[i] * MatrixH.dydeta[i]) - (MatrixH.dxdeta[i] * MatrixH.dydksi[i])
            MatrixH.revJacobian.append(1 / temp)
            MatrixH.matJacobian.append(temp)
    elif number == 3:
        for i in range(number ** 2):
            temp = (MatrixH.dxdksi[i] * MatrixH.dydeta[i]) - (MatrixH.dxdeta[i] * MatrixH.dydksi[i])
            MatrixH.revJacobian.append(1 / temp)
            MatrixH.matJacobian.append(temp)
    # print(MatrixH.revJacobian)
    return MatrixH.revJacobian, MatrixH.matJacobian


def translate(number):  # 80 0 0 80
    xksi, xeta, yksi, yeta = 0.0, 0.0, 0.0, 0.0
    tr = []
    for i in range(number**2):
        tr.append([])
        for j in range(4):
            xksi = xksi + (El4.ksi[i][j] * MatrixH.coordinates[j][0])
            xeta = xeta + (El4.eta[i][j] * MatrixH.coordinates[j][0])
            yksi = yksi + (El4.ksi[i][j] * MatrixH.coordinates[j][1])
            yeta = yeta + (El4.eta[i][j] * MatrixH.coordinates[j][1])
        tr[i].append(xksi * MatrixH.revJacobian[i])
        tr[i].append(xeta * MatrixH.revJacobian[i])
        tr[i].append(yksi * MatrixH.revJacobian[i])
        tr[i].append(yeta * MatrixH.revJacobian[i])

    # print(tr)
    return tr


def transform(number):  # dla punktow calkowania dndx, dndy
    if number == 2:
        for i in range(number**2):
            MatrixH.nx.append([])
            MatrixH.ny.append([])
            for j in range(4):
                MatrixH.nx[i].append(translate(2)[0][0] * El4.ksi[i][j] + translate(2)[0][1] * El4.eta[i][j])
                MatrixH.ny[i].append(translate(2)[0][2] * El4.ksi[i][j] + translate(2)[0][3] * El4.eta[i][j])
    if number == 3:
        for i in range(number**2):
            MatrixH.nx.append([])
            MatrixH.ny.append([])
            for j in range(4):
                MatrixH.nx[i].append(translate(3)[0][0] * El4.ksi[i][j] + translate(3)[0][1] * El4.eta[i][j])
                MatrixH.ny[i].append(translate(3)[0][2] * El4.ksi[i][j] + translate(3)[0][3] * El4.eta[i][j])


    # print(MatrixH.nx), "\n\n", np.array(MatrixH.ny))
    # print(np.array(MatrixH.nx))
    return MatrixH.nx, MatrixH.ny


def calcMatrixH(kt, number):  # calka kt((nx*nx-1)+(ny*ny-1))
    tempArr = []
    temp2 = 0
    # print(np.array(MatrixH.nx), "\n\n", np.array(MatrixH.ny))
    if number == 2:
        for k in range(number**2):
            tempArr.append([])
            for i in range(number**2):
                # tempArr.append([])
                for j in range(4):
                    temp = (MatrixH.nx[k][j] * MatrixH.nx[k][i] + MatrixH.ny[k][j] * MatrixH.ny[k][i]) * \
                           MatrixH.matJacobian[i]
                    temp = temp * kt
                    tempArr[k].append(temp)
    # MatrixH.H.insert(-1, tempArr)
    # print(np.array(tempArr))
    # print(np.array(MatrixH.H))
        for i in range(number ** 2):
            tempArr2 = np.array(tempArr[i]).reshape(4, 4)
            MatrixH.H.append(tempArr2)
    elif number == 3:
        for k in range(number**2):
            tempArr.append([])
            for i in range(4):
                for j in range(4):
                    temp = (MatrixH.nx[k][j] * MatrixH.nx[k][i] + MatrixH.ny[k][j] * MatrixH.ny[k][i]) * \
                           MatrixH.matJacobian[k]
                    temp = temp * kt
                    tempArr[k].append(temp)
        for i in range(number ** 2):
            tempArr2 = np.array(tempArr[i]).reshape(4, 4)
            MatrixH.H.append(tempArr2)

    # print(np.array(MatrixH.H))
    return MatrixH.H


def H(number):
    matH = []
    tempp = 0
    if number == 2:
        for j in range(4):
            temp = 0
            for k in range(number**2):
                temp = temp + MatrixH.H[k][j] * SC.ptWeight2[k % 2] * SC.ptWeight2[k % 2]
            # print(temp)
            matH.append(temp)
    elif number == 3:
        for j in range(4):
            temp = 0
            for k in range(number**2):
                # print(MatrixH.H)
                if k > 2:
                    tempp = 1
                if k > 5:
                    tempp = 2
                temp = temp + MatrixH.H[k][j] * SC.ptWeight3[tempp] * SC.ptWeight3[k % 3]
            # print(temp)
            matH.append(temp)
    print(np.array(matH))
    return matH

# =============================== Input Functions
def getGlobalData(input_file):
    with open(input_file, 'r') as f:
        GlobalData.simTime = int(f.readline().split()[-1])
        GlobalData.simStepTime = int(f.readline().split()[-1])
        GlobalData.conductivity = int(f.readline().split()[-1])
        GlobalData.alpha = int(f.readline().split()[-1])
        GlobalData.tot = int(f.readline().split()[-1])
        GlobalData.initialTemp = int(f.readline().split()[-1])
        GlobalData.density = int(f.readline().split()[-1])
        GlobalData.specificHeat = int(f.readline().split()[-1])

    # print(GlobalData.simTime, GlobalData.simStepTime, GlobalData.conductivity, GlobalData.alpha,
    #      GlobalData.tot, GlobalData.initialTemp, GlobalData.density, GlobalData.specificHeat)  # control print
    return GlobalData.simTime, GlobalData.simStepTime, GlobalData.conductivity, GlobalData.alpha, \
           GlobalData.tot, GlobalData.initialTemp, GlobalData.density, GlobalData.specificHeat


def getGridData(input_file):
    with open(input_file) as f:
        for i in range(8):  # skip first data sets, which are imported already
            f.readline()
        Grid.nodesNumber = int(f.readline().split()[-1])
        Grid.elementsNumber = int(f.readline().split()[-1])

        f.readline()

        for i in range(Grid.nodesNumber):
            temp = f.readline().replace(",", "").strip().split()
            temp.pop(0)  # removing line number
            temp = [float(temp) for temp in temp]  # let's float them
            Grid.nodes.append(Node(temp[0], temp[1]))
        # print(Grid.nodes[0].x)  # control print

        f.readline()

        for i in range(Grid.elementsNumber):
            temp = f.readline().replace(",", "").strip().split()
            temp.pop(0)
            temp = [float(temp) for temp in temp]
            Grid.elements.append(Element(temp))
        # print(Grid.elements[0].ID[2])  # control print

        f.readline()

        temp = f.readline().replace(",", "").strip().split()
        temp = [int(temp) for temp in temp]
        Node.BC = [False for i in range(Grid.nodesNumber)]

        for i in range(len(temp)):
            num = temp[i]
            Node.BC[num - 1] = True

        # print(Node.BC)  # control print

        return Grid.nodes, Grid.elements, Node.BC


# =============================== Solving Functions
def gaussianQuadrature(number, dimension):
    result = 0.0

    if dimension == 1:
        if number == 2:
            for i in range(number):
                result = result + f1(c.intPt2[i]) * c.ptWeight2[i]
        elif number == 3:
            for i in range(number):
                result = result + f1(c.intPt3[i]) * c.ptWeight3[i]
    elif dimension == 2:
        if number == 2:
            for i in range(number):
                for j in range(number):
                    result = result + f2(c.intPt2[i], c.intPt2[j]) * c.ptWeight2[i] * c.ptWeight2[j]
        elif number == 3:
            for i in range(number):
                for j in range(number):
                    result = result + f2(c.intPt3[i], c.intPt3[j]) * c.ptWeight3[i] * c.ptWeight3[j]
    return result


def N1(x):
    return (1 - x) / 2


def N2(x):
    return (1 + x) / 2


def integral(number, x1, x2):
    result = 0

    detJ = (x2 - x1) / 2

    PC = [(N1(SC.intPt2[0]) * x1 + N2(SC.intPt2[0]) * x2), (N1(SC.intPt2[1]) * x1 + N2(SC.intPt2[1]) * x2)]

    if number == 2:
        for i in range(number):
            result = result + f(PC[i]) * SC.ptWeight2[i]
    # elif number == 3:
    #    for i in range(number):
    #        result = result + f(PC[i]) * SC.ptWeight3[i]

    return result * detJ


def Nksi(number, x):
    if number == 0:
        return 0.25 * (x - 1)
    if number == 1:
        return -0.25 * (x - 1)
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


# =============================== To be solved
def f(x):
    return x + 2


def f1(x):
    return 2 * x ** 2 + 3 * x - 8


def f2(ksi, eta):
    return -5 * ksi ** 2 * eta + 2 * ksi * eta ** 2 + 10


def f3(x):
    return 3 * x ** 2 - 6 * x + 1


def mes(number):
    for i in range(Grid.nodesNumber):
        MatrixH.GH.append([])
        for j in range(Grid.nodesNumber):
            MatrixH.GH[i].append(0)

    for i in range(Grid.elementsNumber):
        for j in range(4):
            temp = int(Grid.elements[i].ID[j]) - 1
            MatrixH.coordinates[j][0] = Grid.nodes[temp].x
            MatrixH.coordinates[j][1] = Grid.nodes[temp].y
        # print(MatrixH.coordinates)

        g.fill(number)  # number # fills ksi and eta arrays
        createMatrix(number)
        jacobian(number)
        translate(number)
        transform(number)
        calcMatrixH(GlobalData.conductivity, number)
        temp = H(number)

        for j in range(4):
            for k in range(4):
                MatrixH.GH[int(Grid.elements[i].ID[j]) - 1][int(Grid.elements[i].ID[k]) - 1] += temp[j][k]

    print(np.array(MatrixH.GH))



# =============================== Main
# Choose input file.

inFile = '1'

# inFile = input("""Choose input file
#            1. Test1_4_4
#            2. Test2_4_4_MixGrid
#            3. Test3_31_31_kwadrat
# File index: """)

if inFile == '1':
    filename = "Test1_4_4.txt"
elif inFile == '2':
    filename = "Test2_4_4_MixGrid.txt"
elif inFile == '3':
    filename = "Test3_31_31_kwadrat.txt"
else:
    print("Error.")
    filename = ""
    exit()

# Class initialization
a = GlobalData()
b = Grid()
c = SC()
# d = Node()
# e = Element()
g = El4()
h = MatrixH()

# Import data from the file.
getGlobalData(filename)
getGridData(filename)

# Printing (tests, we call functions here)
# print("Kwadratura Gaussa: \t", gaussianQuadrature(3,  1))   # number, dimension
# print("Calka: \t\t\t\t", integral(2, 4, 12))                # number, x1, x2
# g.fill(2)  # number # fills ksi and eta arrays
# createMatrix(2)
# jacobian(2)
# translate(2)
# transform(2)
# calcMatrixH(30, 2)
# H(2)
mes(3)

# print(np.matrix(g.ksi), "\n\n", np.matrix(g.eta))