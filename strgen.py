from pysvg.structure import G
from pysvg.shape import Path
from pysvg.builders import StyleBuilder
def strPathStart(coord):
    pathCoord = 'M' + str(float(coord[0])) + ' ' + str(float(coord[1])) + ' '
    return pathCoord

def strPathPoint(coord):
    pathCoord = 'L' + str(float(coord[0])) + ' ' + str(float(coord[1])) + ' '
    return pathCoord

def pGram(pList):
    outString = pList[0] + pList[1] + pList[2] + pList[3] + 'Z'
    return outString

def drawSideFacePair(topCorners, bottomCorners, color):
    pointList = []
    pointList2 = []
    output = G()
    XYStyle1 = StyleBuilder()
    #XYStyle1.setStroke('black')
    #XYStyle1.setStrokeWidth(1.5)
    XYStyle1.setFilling(color)
    #XYStyle1.setFillOpacity(0.8)
    XYStyle2 = StyleBuilder()
    XYStyle2.setStroke(color)
    XYStyle2.setStrokeWidth(0.5)
    XYStyle2.setFilling(color)
    pointList.append(strPathStart(topCorners[2]))
    pointList.append(strPathPoint(topCorners[3]))
    pointList.append(strPathPoint(bottomCorners[3]))
    pointList.append(strPathPoint(bottomCorners[2]))
    svgCode1 = pGram(pointList)
    pointList2.append(strPathStart(topCorners[0]))
    pointList2.append(strPathPoint(topCorners[1]))
    pointList2.append(strPathPoint(bottomCorners[1]))
    pointList2.append(strPathPoint(bottomCorners[0]))
    svgCode2 = pGram(pointList2)
    rightFace = Path(svgCode1)
    rightFace.set_style(XYStyle1.getStyle())
    leftFace = Path(svgCode2)
    leftFace.set_style(XYStyle2.getStyle())
    output.addElement(leftFace)
    output.addElement(rightFace)

    return output








def drawFrontFacePair(topCorners, bottomCorners, color):
    pointList = []
    pointList2 = []
    output = G()
    YZStyle1 = StyleBuilder()
    #YZStyle1.setStroke('black')
    #YZStyle1.setStrokeWidth(1.5)
    YZStyle1.setFilling(color)
    YZStyle1.setFillOpacity(0.75)
    YZStyle2 = StyleBuilder()
    YZStyle2.setStroke('black')
    YZStyle2.setStrokeWidth(0.5)
    YZStyle2.setFillOpacity(0.75)
    YZStyle2.setFilling('#ffffff')
    pointList.append(strPathStart(topCorners[1]))
    pointList.append(strPathPoint(topCorners[2]))
    pointList.append(strPathPoint(bottomCorners[2]))
    pointList.append(strPathPoint(bottomCorners[1]))
    svgCode1 = pGram(pointList)
    pointList2.append(strPathStart(topCorners[0]))
    pointList2.append(strPathPoint(topCorners[3]))
    pointList2.append(strPathPoint(bottomCorners[3]))
    pointList2.append(strPathPoint(bottomCorners[0]))
    svgCode2 = pGram(pointList2)
    frontFace = Path(svgCode1)
    frontFace.set_style(YZStyle1.getStyle())
    backFace = Path(svgCode2)
    backFace.set_style(YZStyle2.getStyle())
    output.addElement(backFace)
    output.addElement(frontFace)

    return output
