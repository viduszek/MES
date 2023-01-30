import numpy as np

import classes
import functions
import data

# Choose input file and number of points:
inFile = '1'
number = 2


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

# Import data from the file.
data.getGlobalData(filename)
data.getGridData(filename)

classes.MatrixH.dxdksi = [0.0]*number*number
classes.MatrixH.dxdeta = [0.0]*number*number
classes.MatrixH.dydksi = [0.0]*number*number
classes.MatrixH.dydeta = [0.0]*number*number

for i in range(classes.Grid.nodesNumber):
    classes.MatrixH.globalP.append(0.0)
    classes.MatrixH.GP.append(0.0)
    classes.MatrixH.T.append(classes.GlobalData.initialTemp)

functions.mes(number)

# print(np.array(classes.Side.resultsA))
# print(np.array(classes.Side.resultsB))
# print(np.array(classes.Side.bc))

# print(np.array(classes.El4.ksi))
# print(np.array(classes.El4.eta))
