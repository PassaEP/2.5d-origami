import numpy as np
import math

def parseNPM2SVG(matrix):
    outList = []
    firstRow = matrix[0]
    secondRow = matrix[1]
    for i in range(0, 3):
        outList.append(str(float(firstRow[0,i])))
        outList.append(str(float(secondRow[0,i])))

    return outList
