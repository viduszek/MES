# =============================== Libraries/Packages
import numpy as np
import scipy as sp

# =============================== Classes


class Node:
    x = float()
    y = float()
    BC = bool()

    def __init__(self):
        print("class: Node initialized.")


class Element:
    ID = []  # up to 4!
    # ID = [4]

    def __init__(self):
        print("class: Element initialized.")


class GlobalData:  # there's actually no need to declare those variables, but you can do it tho
    simTime = int()
    simStepTime = int()
    conductivity = int()
    alpha = int()
    tot = int()
    initialTemp = int()
    density = int()
    specificHeat = int()

    def __init__(self):
        print("class: GlobalData initialized.")


class Grid:
    nodesNumber = int()
    elementsNumber = int()
    nodes = []
    elements = []

    def __init__(self):
        print("class: Grid initialized.")


class SC(object):
    intPt2 = [-1 / np.sqrt(3), 1 / np.sqrt(3)]  # integral points
    ptWeight2 = [1, 1]

    intPt3 = [-np.sqrt(0.6), 0, np.sqrt(0.6)]  # integral points (#3)
    ptWeight3 = [0.555555, 0.888888, 0.555555]

    def __init__(self):
        print("class: SC initialized.")


class Element4(object):
    ksi = []
    eta = []

    def __init__(self):
        print("class: Element4 initizalized")


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

    print(GlobalData.simTime, GlobalData.simStepTime, GlobalData.conductivity, GlobalData.alpha,
          GlobalData.tot, GlobalData.initialTemp, GlobalData.density, GlobalData.specificHeat)
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
        print(Grid.nodes)

        f.readline()

        for i in range(elementsNumber):
            temp = f.readline().replace(",", "").strip().split()
            temp.pop(0)
            temp = [float(temp) for temp in temp]
            Grid.elements.append(temp)
        print(Grid.elements)

        return Grid.nodes, Grid.elements


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


# =============================== To be solved

def f1(x):
    return 2 * x**2 + 3 * x - 8


def f2(ksi, eta):
    return -5 * ksi**2 * eta + 2 * ksi * eta**2 + 10


def f3(x):
    return 3 * x**2 - 6 * x + 1


# =============================== Main

# Choose input file.
filename = "Test1_4_4.txt"
# filename = "Test2_4_4_MixGrid.txt"
# filename = "Test3_31_31_kwadrat.txt"

# Class initialization
a = GlobalData()
b = Grid()
c = SC()
d = Node()
e = Element()
f = Element4()

# Import data from the file.
getGlobalData(filename)
getGridData(filename)

print(gaussianQuadrature(3,  1))
