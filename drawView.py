from pysvg.structure import G
from pysvg.builders import ShapeBuilder, TransformBuilder, StyleBuilder
from pysvg.text import *
from pysvg.shape import *
from strgen import strPathStart, strPathPoint, pGram
import math

FRAMELEN = 200

SPACING = 30
SPACING_AXIS = 80
AXIS_LENGTH = 20
PRISM_SPACING = 250
SCALING_FACTOR = 2
CIRCLE_RADIUS = 1.125*SCALING_FACTOR
ROW_SPACE = 3.375*SCALING_FACTOR
LABEL_SHIFT = 15


class Slice:
    def __init__(self, dimA, dimB, corner, color):
        self.dimA = float(dimA)
        self.dimB = float(dimB)
        self.gridCoord = corner
        self.color = str(color)

    def calcPosition(self):
        self.gridCoord[0] = self.gridCoord[0] + (FRAMELEN - SCALING_FACTOR*self.dimA)/2                               # finetuned later
        self.gridCoord[1] = self.gridCoord[1] + (FRAMELEN - SCALING_FACTOR*self.dimB)/2


    def draw(self):
        rectangleStyle = StyleBuilder()
        rectangleStyle.setStroke('black')
        rectangleStyle.setStrokeWidth(2)
        rectangleStyle.setFilling(self.color)
        self.view = Rect(self.gridCoord[0], self.gridCoord[1], self.dimA*SCALING_FACTOR, self.dimB*SCALING_FACTOR)
        self.view.set_style(rectangleStyle.getStyle())
        return self.view

    def genLabels(self):
        labels = G()
        labelStyle = StyleBuilder()
        labelStyle.setTextAnchor('middle')
        labelStyle.setFontSize(10)
        lineStyle = StyleBuilder()
        lineStyle.setStroke('#a9a9a9')
        dimLabelA = str(self.dimA) + ' nm'
        dimLabelB = str(self.dimB) + ' nm'
        # horizontal label
        cornerHL = self.view.getTopLeft()
        cornerHR = self.view.getTopRight()

        hLabelStart = [cornerHL[0], cornerHL[1] + 5]
        hLabelEnd = [cornerHR[0], cornerHR[1] + 5]

        hLine = Line(hLabelStart[0], hLabelStart[1], hLabelEnd[0], hLabelEnd[1])
        hLabel = Text(dimLabelA, (hLabelStart[0] + hLabelEnd[0])/2, hLabelEnd[1] + 20)
        hLine.set_style(lineStyle.getStyle())
        hLabel.set_style(labelStyle.getStyle())

        cornerVL = self.view.getTopLeft()
        cornerVR = self.view.getBottomLeft()

        vLabelStart = [cornerVL[0] - 5, cornerVL[1]]
        vLabelEnd = [cornerVR[0] - 5, cornerVR[1]]

        vLine = Line(vLabelStart[0], vLabelStart[1], vLabelEnd[0], vLabelEnd[1])
        vLabel = Text(dimLabelB, vLabelEnd[0] - 30, (vLabelStart[1] + vLabelEnd[1])/2)
        vLine.set_style(lineStyle.getStyle())
        vLabel.set_style(labelStyle.getStyle())
        labels.addElement(hLine)
        labels.addElement(hLabel)
        labels.addElement(vLine)
        labels.addElement(vLabel)

        #center = [self.gridCoord[0] + self.dimA*SCALING_FACTOR/2, self.gridCoord[1] + self.dimB*SCALING_FACTOR + LABEL_SHIFT]
        #dimLabel = str(self.dimA) + 'nm x ' +str(self.dimB) + 'nm'
        #labels = Text(dimLabel, center[0], center[1])
        #labels.set_style(labelStyle.getStyle())
        return labels


class Prism:
    def __init__(self, x, y, z, YZgridCoord):
        self.x = float(x)*SCALING_FACTOR
        self.y = float(y)*SCALING_FACTOR
        self.z = float(z)*SCALING_FACTOR
        print(str(self.z) + ' prism')
        self.startCoord = [YZgridCoord[0], YZgridCoord[1] + self.y + SPACING]
        self.out = G()
    def draw(self):
        faceStyle = StyleBuilder()
        faceStyle.setFilling('blue')
        faceStyle.setStroke('black')
        faceStyle.setFillOpacity(0.5)
        topStyle = StyleBuilder()
        topStyle.setFilling('green')
        topStyle.setStroke('black')
        topStyle.setStrokeWidth(2)
        topStyle.setFillOpacity(0.5)
        sideStyle = StyleBuilder()
        sideStyle.setFilling('red')
        sideStyle.setStroke('black')
        sideStyle.setStrokeWidth(2)
        sideStyle.setFillOpacity(0.5)
        frontFace = Rect(self.startCoord[0], self.startCoord[1], self.x, self.z)
        backFace = Rect(self.startCoord[0] + self.y/math.sqrt(2), self.startCoord[1] - self.y/math.sqrt(2), self.x, self.z)
        # draw parallelogram, grab corners
        cornerBack = backFace.getEdgePoints()
        cornerFront = frontFace.getEdgePoints()
        # Top face string
        pointListT = []
        pointListT.append(strPathStart(cornerBack[0]))
        pointListT.append(strPathPoint(cornerBack[1]))
        pointListT.append(strPathPoint(cornerFront[1]))
        pointListT.append(strPathPoint(cornerFront[0]))
        topFaceString = pGram(pointListT)
        # right side face
        pointListR = [strPathStart(cornerBack[1])]
        pointListR.append(strPathPoint(cornerBack[2]))
        pointListR.append(strPathPoint(cornerFront[2]))
        pointListR.append(strPathPoint(cornerFront[1]))
        # bottom face
        pointListB = [strPathStart(cornerBack[2])]
        pointListB.append(strPathPoint(cornerBack[3]))
        pointListB.append(strPathPoint(cornerFront[3]))
        pointListB.append(strPathPoint(cornerFront[2]))
        # left
        pointListL = [strPathStart(cornerBack[0])]
        pointListL.append(strPathPoint(cornerBack[3]))
        pointListL.append(strPathPoint(cornerFront[3]))
        pointListL.append(strPathPoint(cornerFront[0]))
        # make faces
        topFace = Path(pGram(pointListT))
        rightFace = Path(pGram(pointListR))
        leftFace = Path(pGram(pointListL))
        bottomFace = Path(pGram(pointListB))

        topFace.set_style(topStyle.getStyle())
        bottomFace.set_style(topStyle.getStyle())

        rightFace.set_style(sideStyle.getStyle())
        leftFace.set_style(sideStyle.getStyle())

        frontFace.set_style(faceStyle.getStyle())
        backFace.set_style(faceStyle.getStyle())


        self.out.addElement(backFace)
        self.out.addElement(frontFace)

        self.out.addElement(rightFace)
        self.out.addElement(leftFace)

        self.out.addElement(bottomFace)
        self.out.addElement(topFace)

        return self.out

class Honeycomb:
    def __init__(self, row, col, XZgridCoord, x,z):
        self.col = int(col)
        self.row = int(row)
        self.x = float(x)*SCALING_FACTOR#4*CIRCLE_RADIUS*((self.col-2)/2) + 2*CIRCLE_RADIUS
        self.z = float(z)*SCALING_FACTOR
        self.gridCoord = [XZgridCoord[0] + (FRAMELEN - self.x)/2,XZgridCoord[1]+ (FRAMELEN - self.z)/2]
        self.hCSlice = Rect(self.gridCoord[0],self.gridCoord[1] , self.x, self.z)
    def genHoneycomb(self):
        outComb = G()
        circleStyle = StyleBuilder()
        #circleStyle.setStrokeWidth(0.5)
        circleStyle.setStroke('orange')
        circleStyle.setFilling('#edd239')
        CORNER = self.hCSlice.getBottomRight()

        circleY = []
        circleY.append(CORNER[1] + (ROW_SPACE - CIRCLE_RADIUS))
        for g in range(0, self.row):
            circleX = CORNER[0] - CIRCLE_RADIUS
            if (g != 0 and g % 2 != 0):
                circleY.append(circleY[g-1] + 2*CIRCLE_RADIUS)
            elif (g!= 0 and g % 2 == 0):
                circleY.append(circleY[g-1] + 2*(ROW_SPACE - CIRCLE_RADIUS))

            #initCircle = Circle(circleX, circleY[g], CIRCLE_RADIUS)
            #outComb.addElement(initCircle)
            for i in range(0, self.col):
                newShift = TransformBuilder()
                xShift = -i*CIRCLE_RADIUS*math.sqrt(3)
                yShift = 0
                if i % 2 != 0 and g % 2 != 0:
                    yShift = CIRCLE_RADIUS
                elif i % 2 != 0 and g % 2 == 0:
                    yShift = -1*CIRCLE_RADIUS
                newShift.setTranslation(str(xShift) + ' ' + str(yShift))
                a = Circle(circleX, circleY[g], CIRCLE_RADIUS)
                a.set_transform(newShift.getTransform())
                outComb.addElement(a)

        outComb.set_style(circleStyle.getStyle())
        return outComb
    def genLabels(self):

        labels = G()
        labelStyle = StyleBuilder()
        labelStyle.setTextAnchor('middle')
        labelStyle.setFontSize(10)
        lineStyle = StyleBuilder()
        lineStyle.setStroke('#a9a9a9')
        dimLabelA = str(self.x/SCALING_FACTOR) + ' nm'
        dimLabelB = str(self.z/SCALING_FACTOR) + ' nm'
        # horizontal label
        cornerHL = self.hCSlice.getTopLeft()
        cornerHR = self.hCSlice.getTopRight()

        hLabelStart = [cornerHL[0], cornerHL[1] + 5]
        hLabelEnd = [cornerHR[0], cornerHR[1] + 5]

        hLine = Line(hLabelStart[0], hLabelStart[1], hLabelEnd[0], hLabelEnd[1])
        hLabel = Text(dimLabelA, (hLabelStart[0] + hLabelEnd[0])/2, hLabelEnd[1] + 20)
        hLine.set_style(lineStyle.getStyle())
        hLabel.set_style(labelStyle.getStyle())

        cornerVL = self.hCSlice.getTopLeft()
        cornerVR = self.hCSlice.getBottomLeft()

        vLabelStart = [cornerVL[0] - 5, cornerVL[1]]
        vLabelEnd = [cornerVR[0] - 5, cornerVR[1]]

        vLine = Line(vLabelStart[0], vLabelStart[1], vLabelEnd[0], vLabelEnd[1])
        vLabel = Text(dimLabelB, vLabelEnd[0] - 30, (vLabelStart[1] + vLabelEnd[1])/2)
        vLine.set_style(lineStyle.getStyle())
        vLabel.set_style(labelStyle.getStyle())
        labels.addElement(hLine)
        labels.addElement(hLabel)
        labels.addElement(vLine)
        labels.addElement(vLabel)
        return labels





    def genHoneyRect(self):
        rectangleStyle = StyleBuilder()
        rectangleStyle.setFilling('blue')
        rectangleStyle.setStroke('black') # temp
        self.hCSlice.set_style(rectangleStyle.getStyle())
        return(self.hCSlice)
