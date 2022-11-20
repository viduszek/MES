# =============================== Libraries/Packages
import numpy as np


# =============================== Classes
class Node:
    x = float()
    y = float()
    BC = bool()


class Element:
    ID = []  # up to 4!
    # ID = [4]


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
    ptWeight2 = [1, 1]                          # weights

    intPt3 = [-np.sqrt(3/5), 0, np.sqrt(3/5)]  # points for n = 3
    ptWeight3 = [5.0/9.0, 8.0/9.0, 5.0/9.0]    # weights

    intPt4 = []
    ptWeight4 = []

    # def __init__(self):
        # print("class: SC initialized.")


class El4:
    ksi = [] 
    eta = []

    ksi2 = [1, -1, -1, 1]  # additional weights for the integral points / ksi version
    # ksi3 = []
    # ksi4 = []

    eta2 = [1, 1, 1, 1]  # additional weights for the integral points / eta version
    # eta3 = []
    # eta4 = []

    @staticmethod
    def fill(number):
        if number == 2:
            for i in range(number**2):
                El4.ksi.append([])
                El4.eta.append([])
                for j in range(4):
                    El4.ksi[i].append(Nksi(j, SC.intPt2[i % number] * El4.ksi2[i]))
                    El4.eta[i].append(Neta(j, SC.intPt2[i % number] * El4.eta2[i]))
        # elif number == 3:
        #     print("tu jeszcze nic nie ma.")
        # elif number == 4:
        #     print("tu tez jeszcze nic nie ma")


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
        nodesNumber = int(f.readline().split()[-1])
        elementsNumber = int(f.readline().split()[-1])

        f.readline()

        for i in range(nodesNumber):
            temp = f.readline().replace(",", "").strip().split()
            temp.pop(0)  # removing line number
            temp = [float(temp) for temp in temp]  # let's make them float
            Grid.nodes.append(temp)
        # print(Grid.nodes)  # control print

        f.readline()

        for i in range(elementsNumber):
            temp = f.readline().replace(",", "").strip().split()
            temp.pop(0)
            temp = [float(temp) for temp in temp]
            Grid.elements.append(temp)
        # print(Grid.elements)  # control print

        f.readline()

        temp = f.readline().replace(",", "").strip().split()
        temp = [int(temp) for temp in temp]
        Node.BC = [False for i in range(nodesNumber)]

        for i in range(len(temp)):
            num = temp[i]
            Node.BC[num-1] = True

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
    return (1-x)/2


def N2(x):
    return (1+x)/2


def integral(number, x1, x2):
    result = 0

    detJ = (x2 - x1)/2

    PC = [(N1(SC.intPt2[0])*x1+N2(SC.intPt2[0])*x2), (N1(SC.intPt2[1])*x1+N2(SC.intPt2[1])*x2)]

    if number == 2:
        for i in range(number):
            result = result + f(PC[i]) * SC.ptWeight2[i]
    # elif number == 3:
    #    for i in range(number):
    #        result = result + f(PC[i]) * SC.ptWeight3[i]

    return result * detJ


def Nksi(number, x):
    if number == 0:
        return 0.25 * (x-1)
    if number == 1:
        return -0.25 * (x-1)
    if number == 2:
        return 0.25 * (x+1)
    if number == 3:
        return -0.25 * (x+1)
# 0.25(1-k)(1-n), 0.25(1+k)(1-n), 0.25(1+k)(1+n), 0.25(1-k)(1+n)


def Neta(number, x):
    if number == 0:
        return -0.25 * (1-x)
    if number == 1:
        return -0.25 * (1+x)
    if number == 2:
        return 0.25 * (1+x)
    if number == 3:
        return 0.25 * (1-x)


# =============================== To be solved
def f(x):
    return x + 2


def f1(x):
    return 2 * x**2 + 3 * x - 8


def f2(ksi, eta):
    return -5 * ksi**2 * eta + 2 * ksi * eta**2 + 10


def f3(x):
    return 3 * x**2 - 6 * x + 1


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
d = Node()
e = Element()
g = El4()

# Import data from the file.
getGlobalData(filename)
getGridData(filename)

# Printing (tests, we call functions here)
print("Kwadratura Gaussa: \t", gaussianQuadrature(3,  1))   # number, dimension
print("Calka: \t\t\t\t", integral(2, 4, 12))                # number, x1, x2
g.fill(2)                                                   # number # fills ksi and eta arrays

print(np.matrix(El4.ksi), "\n\n", np.matrix(El4.eta))
