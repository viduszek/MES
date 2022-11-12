import numpy


class Node(object):
    x = float(0)
    y = float(0)
    temp = float(0)
    WB = bool(False)


class Element(object):
    ID = []


class GlobalData(object):
    simTime = int(0)
    simStepTime = int(0)
    conductivity = int(0)
    alpha = int(0)
    tot = int(0)
    initialTemperature = int(0)
    density = int(0)
    specificHeat = int(0)


class Grid(object):
    nodesNumber = int(0)
    elementsNumber = int(0)
