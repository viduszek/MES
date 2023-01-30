import classes


def getGlobalData(input_file):
    with open(input_file, 'r') as f:
        classes.GlobalData.simTime = int(f.readline().split()[-1])
        classes.GlobalData.simStepTime = int(f.readline().split()[-1])
        classes.GlobalData.conductivity = int(f.readline().split()[-1])
        classes.GlobalData.alpha = int(f.readline().split()[-1])
        classes.GlobalData.tot = int(f.readline().split()[-1])
        classes.GlobalData.initialTemp = int(f.readline().split()[-1])
        classes.GlobalData.density = int(f.readline().split()[-1])
        classes.GlobalData.specificHeat = int(f.readline().split()[-1])

    # print(GlobalData.simTime, GlobalData.simStepTime, GlobalData.conductivity, GlobalData.alpha,
    #      GlobalData.tot, GlobalData.initialTemp, GlobalData.density, GlobalData.specificHeat)  # control print
    return classes.GlobalData.simTime, classes.GlobalData.simStepTime, classes.GlobalData.conductivity, \
           classes.GlobalData.alpha, classes.GlobalData.tot, classes.GlobalData.initialTemp, \
           classes.GlobalData.density, classes.GlobalData.specificHeat


def getGridData(input_file):
    with open(input_file) as f:
        for i in range(8):  # skip first data sets, which are imported already
            f.readline()
        classes.Grid.nodesNumber = int(f.readline().split()[-1])
        classes.Grid.elementsNumber = int(f.readline().split()[-1])

        f.readline()

        for i in range(classes.Grid.nodesNumber):
            temp = f.readline().replace(",", "").strip().split()
            temp.pop(0)  # removing line number
            temp = [float(temp) for temp in temp]  # let's float them
            classes.Grid.nodes.append(classes.Node(temp[0], temp[1]))
        # print(Grid.nodes[0].x)  # control print

        f.readline()

        for i in range(classes.Grid.elementsNumber):
            temp = f.readline().replace(",", "").strip().split()
            temp.pop(0)
            temp = [float(temp) for temp in temp]
            classes.Grid.elements.append(classes.Element(temp))
        # print(Grid.elements[0].ID[2])  # control print

        f.readline()

        temp = f.readline().replace(",", "").strip().split()
        temp = [int(temp) for temp in temp]
        # classes.Node[i].BC = [False for i in range(classes.Grid.nodesNumber)]

        for i in range(classes.Grid.nodesNumber):
            classes.Grid.nodes[i].BC = False

        for i in range(len(temp)):
            num = temp[i]
            classes.Grid.nodes[num-1].BC = True

        # print(classes.Node.BC)  # control print

        return classes.Grid.nodes, classes.Grid.elements, classes.Node.BC