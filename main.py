import numpy as np


class Node(object):
    x = float()
    y = float()
    BC = bool()


class Element(object):
    ID = [4]


class GlobalData(object):
    global simTime
    global simStepTime
    global conductivity
    global alpha
    global tot
    global initialTemperature
    global density
    global specificHeat


    def getdata(filename):
        with open(filename, 'r') as f:
            simTime = f.readline().split()[-1]
            simStepTime = f.readline().split()[-1]
            conductivity = f.readline().split()[-1]
            alpha = f.readline().split()[-1]
            tot = f.readline().split()[-1]
            initialTemperature = f.readline().split()[-1]
            density = f.readline().split()[-1]
            specificHeat = f.readline().split()[-1]
            print(simTime, simStepTime, conductivity, alpha, tot, initialTemperature, density, specificHeat)
            return simTime, simStepTime, conductivity, alpha, tot, initialTemperature, density, specificHeat


class Grid(object):
    nodesNumber = int()
    elementsNumber = int()
    nodes = []
    elements = []

    def getdata(filename):
        with open(filename) as f:
            for i in range(8):          # skip first data sets, which are imported already
                trash = f.readline()
            nodesNumber = int(f.readline().split()[-1])
            elementsNumber = int(f.readline().split()[-1])

            trash = f.readline()

            for i in range(nodesNumber):
                temp = f.readline().replace(",", "").strip().split()
                temp.pop(0)  # removing line number
                Grid.nodes.append(temp)
                np.array(Grid.nodes, dtype=float)
            print(Grid.nodes)

            trash = f.readline()

            for i in range(elementsNumber):
                temp = f.readline().replace(",", "").strip().split()
                temp.pop(0)  # removing line number
                Grid.elements.append(temp)
                np.array(Grid.elements, dtype=float)
            print(Grid.elements)


class SC(object):
    intPt2 = [-1/np.sqrt(3), 1/np.sqrt(3)]  # integral points
    ptWeight2 = [1, 1]

    intPt3 = [-np.sqrt(0.6), 0, np.sqrt(0.6)]  # integral points (#3)
    ptWeight3 = [0.555555, 0.888888, 0.555555]


class Element4(object):
    ksi = []
    eta = []


# Choose input file:
filename = "Test1_4_4.txt"
# filename = "Test2_4_4_MixGrid.txt"
# filename = "Test3_31_31_kwadrat.txt"

# Import data from the file.
GlobalData.getdata(filename)
Grid.getdata(filename)
