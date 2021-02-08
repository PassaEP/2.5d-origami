
SCALING_FACTOR = 2
FRAMELEN = 200
AXIS_LENGTH = 35
from pysvg.builders import ShapeBuilder
from pysvg.text import *
from pysvg.structure import Svg
from pysvg.shape import *
from drawView import *
from page import *
import numpy as np
import math
from strgen import *
from handleMatrix import *

CIRCLE_RADIUS = 1.125*SCALING_FACTOR
ROW_SPACE = 3.375*SCALING_FACTOR

#check the bracket problem
class View:
    def __init__(self, x, y, z, frame, angle):
        self.x = SCALING_FACTOR*float(x)
        self.y = SCALING_FACTOR*float(y)
        self.z = SCALING_FACTOR*float(z)
        self.angle = angle

        self.frame = frame
        self.out = G()
        # figure out the aspect ratio stuff
        self.syms = {'red': 'x', 'blue': 'y', 'green': 'z'}
        self.dimDict = {'x': self.x / SCALING_FACTOR, 'y': self.y / SCALING_FACTOR, 'z': self.z/SCALING_FACTOR}
        dims = [self.x, self.y, self.z]
        areas = [self. x * self.z, self.x * self.y, self.y*self.z]
        dimPairs = [[self.x, self.z], [self.x, self.y], [self.y,self.z]]
        colors = ['blue', 'green', 'red']
        largestFace = max(areas)
        index = areas.index(largestFace)
        self.topDims = dimPairs[index]
        dims.remove(self.topDims[0])
        dims.remove(self.topDims[1])
        self.height = dims[0]
        self.topColor = colors[index]
        self.sideColors = colors[areas.index(self.height*self.topDims[0])]
        self.frontColors = colors[areas.index(self.height*self.topDims[1])]
        # initializing transformation matrices

        cosAngle = math.cos(math.radians(self.angle))
        sinAngle = math.sin(math.radians(self.angle))
        tanAngle = math.tan(math.radians(self.angle))

        m_scaleY = np.matrix([[1, 0, 0],[0, cosAngle, 0],[0, 0, 1]])

        m_shearX = np.matrix([[1, tanAngle, 0],[0, 1, 0],[0, 0, 1]])

        m_rotate = np.matrix([[cosAngle, sinAngle, 0],[-sinAngle, cosAngle, 0],[0, 0, 1]])

        self.transformMatrix = m_rotate * (m_shearX * m_scaleY)
        self.transformMatrix[0, 1] = float(self.transformMatrix[0, 1]) * -1
        self.transformMatrix[1, 0] = float(self.transformMatrix[1, 0]) * -1

        print(self.transformMatrix)

        self.invtransform = np.linalg.inv(self.transformMatrix)
        self.topDownTransform = TransformBuilder()
        stringMatrix = parseNPM2SVG(self.transformMatrix)
        self.topDownTransform.setMatrix(stringMatrix[0],stringMatrix[1],stringMatrix[2],stringMatrix[3],stringMatrix[4],stringMatrix[5])



    def CalcStart(self):
        transformedMidpoint = np.array([[self.frame[0] + FRAMELEN/2], [self.frame[1] + FRAMELEN/2], [1]])
        untransformedMidpoint = self.invtransform * transformedMidpoint

        # getting bottom left to start initializing rectangle
        self.startCoords = [float(untransformedMidpoint[0]) - 0.5*self.topDims[0], float(untransformedMidpoint[1]) - 0.5*self.topDims[1]]





    def makeFaces(self):
        #wantedCoords = np.array([[self.start[0]],[self.start[1]],[1]])
        #self.startCoords = self.invtransform * wantedCoords
        topFace = Rect(float(self.startCoords[0]), float(self.startCoords[1]), self.topDims[0], self.topDims[1])

        vertTransform = np.matrix([[1, 0, 0], [0.577,1,0], [0,0,1]])

        oldYLength = np.matrix([[0],[self.height],[1]])
        # transform Y to get coords for bottom face parts
        #transformedYLengthV = vertTransform * oldYLength

        #transformedYLength = float(transformedYLengthV[1])
        points = topFace.getEdgePoints()

        self.transformedTop = []
        self.transformedBottom = []
        # get untransformed points of the top face
        i = 0
        for stuff in points:
            vector = np.matrix([[stuff[0]], [stuff[1]], [1]])
            self.transformedTop.append(self.transformMatrix * vector)
            self.transformedBottom.append([float(self.transformedTop[i][0]),float(self.transformedTop[i][1]) + self.height])
            i+=1

    def genAxisLabels(self):
        outLabel = G()
        leftArmStart = [self.frame[0] + 0.8 * FRAMELEN, self.frame[1] + 0.1 * FRAMELEN]
        # setting style
        leftStyle = StyleBuilder()
        leftStyle.setStroke(self.sideColors)
        leftStyle.setStrokeWidth(0.5)
        downStyle = StyleBuilder()
        downStyle.setStroke(self.topColor)
        downStyle.setStrokeWidth(0.5)
        rightStyle = StyleBuilder()
        rightStyle.setStroke(self.frontColors)
        rightStyle.setStrokeWidth(0.5)

        labelStyle = StyleBuilder()
        labelStyle.setTextAnchor('middle')
        labelStyle.setFontSize(10)

        # building axis
        leftArmEnd = [leftArmStart[0] +  AXIS_LENGTH * math.cos(math.radians(self.angle)), leftArmStart[1] + AXIS_LENGTH * math.sin(math.radians(self.angle))]
        downArmEnd =  [leftArmEnd[0], leftArmEnd[1] + AXIS_LENGTH]
        rightArmEnd = [leftArmEnd[0] + AXIS_LENGTH * math.cos(math.radians(self.angle)), leftArmEnd[1] - AXIS_LENGTH * math.sin(math.radians(self.angle))]

        # drawing lines and labels
        leftAxis = Line(leftArmStart[0], leftArmStart[1], leftArmEnd[0], leftArmEnd[1])
        leftT = Text(self.syms[self.sideColors], leftArmStart[0], leftArmStart[1] - 5)
        downAxis = Line(leftArmEnd[0], leftArmEnd[1], downArmEnd[0],downArmEnd[1])
        downT = Text(self.syms[self.topColor], downArmEnd[0], downArmEnd[1] + 6)
        rightAxis = Line(leftArmEnd[0], leftArmEnd[1], rightArmEnd[0], rightArmEnd[1])
        rightT = Text(self.syms[self.frontColors], rightArmEnd[0], rightArmEnd[1] - 5)

        leftAxis.set_style(leftStyle.getStyle())
        rightAxis.set_style(rightStyle.getStyle())
        downAxis.set_style(downStyle.getStyle())
        leftT.set_style(labelStyle.getStyle())
        rightT.set_style(labelStyle.getStyle())
        downT.set_style(labelStyle.getStyle())

        outLabel.addElement(leftAxis)
        outLabel.addElement(leftT)
        outLabel.addElement(rightAxis)
        outLabel.addElement(rightT)
        outLabel.addElement(downAxis)
        outLabel.addElement(downT)
        return outLabel

    def callSideTransform(self):
        a = math.cos(math.radians(self.angle))
        b = math.tan(math.radians(self.angle))
        mScaleX = np.matrix([[a, 0, 0], [0, 1, 0], [0, 0, 1]])
        mShearY = np.matrix([[1, 0, 0], [b, 1, 0], [0, 0, 1]])
        newMatrix = mShearY * mScaleX
        newMatrix[0, 1] = float(newMatrix[0, 1]) * -1
        newMatrix[1, 0] = float(newMatrix[1, 0]) * -1
        return newMatrix


    def callFrontTransform(self):
        a = math.cos(math.radians(self.angle))
        b = math.tan(math.radians(-self.angle))
        mScaleX = np.matrix([[a, 0, 0], [0, 1, 0], [0, 0, 1]])
        mShearY = np.matrix([[1, 0, 0], [b, 1, 0], [0, 0, 1]])
        newMatrix = mShearY * mScaleX
        newMatrix[0, 1] = float(newMatrix[0, 1]) * -1
        newMatrix[1, 0] = float(newMatrix[1, 0]) * -1
        return newMatrix


    def helixSlice(self, row, col):
        # figure out if XZ slice is left side, right side, or top panel
        # get the coordinate of where the bottom right corner of the transformed rectangle would be
        # one of the points on top or bottom face
        self.row = row
        self.col = col

        helices = G()

        colors = ['blue', 'green', 'red']

        if self.topColor == colors[0]:
            helixTransform = self.transformMatrix
            helixCorner = np.matrix([[self.startCoords[0] + self.x], [self.startCoords[1]], [1]])
        elif self.sideColors == colors[0]:
            helixTransform = self.callFrontTransform()

            helixCorner = np.matrix([[self.transformedTop[2][0]], [self.transformedTop[2][1]], [1]])
        elif self.frontColors == colors[0]:
            helixTransform = self.callSideTransform()

            helixCorner = np.matrix([[self.transformedTop[3][0]], [self.transformedTop[3][1]], [1]])

        # apply inv transform of that coordinate
        helixinvTransform = np.linalg.inv(helixTransform)

        untransformedBLCorner = helixinvTransform * helixCorner

        stringMatrix = parseNPM2SVG(helixTransform)


        # set the svg matrix

        svgHelixTransform = TransformBuilder()

        svgHelixTransform.setMatrix(stringMatrix[0], stringMatrix[1], stringMatrix[2], stringMatrix[3], stringMatrix[4], stringMatrix[5])


        # initialize the honeycomb lattice with xz coords

        circleStyle = StyleBuilder()
        circleStyle.setStrokeWidth(0.5)
        circleStyle.setStroke('orange')
        circleStyle.setFilling('#edd239')
        CORNER = [float(untransformedBLCorner[0]), float(untransformedBLCorner[1])]


        circleY = []
        circleY.append(CORNER[1] + (ROW_SPACE - CIRCLE_RADIUS))
        for g in range(0, self.row):
            circleX = CORNER[0] - CIRCLE_RADIUS
            if (g != 0 and g % 2 != 0):
                circleY.append(circleY[g-1] + 2*CIRCLE_RADIUS)
            elif (g!= 0 and g % 2 == 0):
                circleY.append(circleY[g-1] + 2*(ROW_SPACE - CIRCLE_RADIUS))

            #initCircle = Circle(circleX, circleY[g], CIRCLE_RADIUS)
            #helices.addElement(initCircle)
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
                helices.addElement(a)

        helices.set_style(circleStyle.getStyle())
        # apply the transform
        helices.set_transform(svgHelixTransform.getTransform())
        return helices


    def drawPerpLabels(self):
        output = G()
        vertLabelCoords = [float(self.transformedTop[2][0] + self.transformedBottom[2][0])/2 + 15, float(self.transformedTop[2][1] + self.transformedBottom[2][1])/2]

        topMiddleCoords = self.transformedTop[2]

        # find the inverse of the top middle corner 2

        untransformedMiddle = self.invtransform * np.array([[topMiddleCoords[0]], [topMiddleCoords[1]], [1]])

        # compute the untransformed coordinates of the right and left labels

        rightLabelCoords = [untransformedMiddle[0] + self.topDims[0], untransformedMiddle[1] + self.height + 15]
        leftLabelCoords = [untransformedMiddle[0] - self.topDims[1], untransformedMiddle[1] + self.height + 15]

        print(rightLabelCoords)
        print(leftLabelCoords)
        # let's generate the style, same for every dim

        perpLabelStyle = StyleBuilder()
        #perpLabelStyle.setTextAnchor('middle')
        perpLabelStyle.setFontSize(7)

        # let's get the transforms
        rightNPTransform = self.callFrontTransform()
        rightStringMatrix = parseNPM2SVG(rightNPTransform)

        leftNPTransform = self.callSideTransform()
        leftStringMatrix = parseNPM2SVG(leftNPTransform)

        # put into SVG format

        rightTransform = TransformBuilder()
        rightTransform.setMatrix(rightStringMatrix[0], rightStringMatrix[1], rightStringMatrix[2], rightStringMatrix[3], rightStringMatrix[4], rightStringMatrix[5])
        leftTransform = TransformBuilder()
        leftTransform.setMatrix(leftStringMatrix[0], leftStringMatrix[1], leftStringMatrix[2], leftStringMatrix[3], leftStringMatrix[4], leftStringMatrix[5])



        # make text objects

        vertLabel = Text(str(self.height / SCALING_FACTOR) + ' nm', vertLabelCoords[0], vertLabelCoords[1])
        rightLabel = Text(str(self.topDims[0] / SCALING_FACTOR) + ' nm', float(rightLabelCoords[0]), float(rightLabelCoords[1]))
        leftLabel = Text(str(self.topDims[1] / SCALING_FACTOR) + ' nm', float(leftLabelCoords[0]), float(leftLabelCoords[1]))

        # apply style

        vertLabel.set_style(perpLabelStyle.getStyle())
        rightLabel.set_style(perpLabelStyle.getStyle())
        leftLabel.set_style(perpLabelStyle.getStyle())

        # apply transforms

        rightLabel.set_transform(rightTransform.getTransform())
        leftLabel.set_transform(leftTransform.getTransform())

        # add to group

        output.addElement(vertLabel)
        output.addElement(rightLabel)
        output.addElement(leftLabel)

        return output


    def Render(self):
        topStyle = StyleBuilder()
        topStyle.setFilling(self.topColor)
        topStyle.setFillOpacity(.75)
        topFace = Rect(float(self.startCoords[0]), float(self.startCoords[1]), self.topDims[0],self.topDims[1])
        topFace.set_transform(self.topDownTransform.getTransform())
        topFace.set_style(topStyle.getStyle())
        self.out.addElement(drawFrontFacePair(self.transformedTop, self.transformedBottom,self.frontColors))
        self.out.addElement(drawSideFacePair(self.transformedTop, self.transformedBottom,self.sideColors))
        self.out.addElement(topFace)
        return self.out
