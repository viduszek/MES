import numpy as np


class Node(object):
    x = float()
    y = float()
    WB = bool()


class Element(object):
    ID = [4]


class GlobalData(object):

    def getdata(filename):
        with open(filename, 'r') as f:
            simTime = int(f.readline().split()[-1])
            simStepTime = int(f.readline().split()[-1])
            conductivity = int(f.readline().split()[-1])
            alpha = int(f.readline().split()[-1])
            tot = int(f.readline().split()[-1])
            initialTemperature = int(f.readline().split()[-1])
            density = int(f.readline().split()[-1])
            specificHeat = int(f.readline().split()[-1])
            # print(simTime, simStepTime, conductivity, alpha, tot, initialTemperature, density, specificHeat)


class Grid(object):
    nodesNumber = int()
    elementsNumber = int()

    def getdata(filename):
        with open(filename, 'r') as f:
            for i in range(8):          # skip first data sets
                trash = f.readline()
            nodesNumber = int(f.readline().split()[-1])
            elementsNumber = int(f.readline().split()[-1])

            trash = f.readline()

            nodes = []
            for i in range(nodesNumber):
                temp = f.readline().replace(",", "").strip().split()
                nodes.append(temp)
            print(nodes)

            trash = f.readline()

            elements = []
            for i in range(elementsNumber):
                temp = f.readline().replace(",", "").strip().split()
                elements.append(temp)
            # print(elements)


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

#
