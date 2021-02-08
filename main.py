import csv
from page import RowFrame
from perpsview import View
from drawView import *
from pysvg.structure import *
import math
import argparse
from blockmaker import cnFileMake

parser = argparse.ArgumentParser(description = 'SVG Doc of Origami Schematics')
parser.add_argument('-n', '--csvFileName', action='store', metavar='', required=True, help='Name of CSV File of Dimensions in directory' )

parser.add_argument('-t', '--targetFileName', action='store', metavar='', required=True, help='SVG Target File Name')

parser.add_argument('-p', '--perpsViewCount', type=int, metavar='', required=True, help='Number of Perspective Views')

parser.add_argument('-a','--angleList', nargs='+', required=True, help='list of angles for perspective view')

parser.add_argument('-sh', dest='showHelix', action='store_true')

parser.set_defaults(showHelix=False)

args = parser.parse_args()



def genSVGFile(csvFile, svgFile, perpsNum, angles, showHelix):
    framesPerRow = 3 + perpsNum
    newDoc = Svg()

    with open(csvFile, 'r', newline='') as csvFile:
        specreader = csv.reader(csvFile)
        next(specreader)
        count = 0
        for row in specreader:
            frames = RowFrame(count, framesPerRow)
            newDoc.addElement(frames.makerow())
            origamiRows = int(row[2])
            origamiCols = int(row[4])
            x, y, z = row[6], row[3], row[8]
            id = int(row[0])

            # some pre processing

            origamiLength = float(row[8])

            origamiDesign = {'id': id, 'rows': origamiRows, 'cols': origamiCols, 'length': origamiLength}

            minX, minY, minZ, maxX, maxY, maxZ, handleEnds = cnFileMake(origamiDesign)

            minXY = [minX, minY]

            maxXY = [maxX, maxY]

            minXZ = [minX, minZ]

            maxXZ = [maxX, maxZ]

            minZY = [minZ, minY]

            maxZY = [maxZ, maxY]

            xySlice = Honeycomb(origamiRows, origamiCols, frames.rowCoords[0], x, y)
            newDoc.addElement(xySlice.genHoneyRect())
            newDoc.addElement(xySlice.genHoneycomb())
            newDoc.addElement(xySlice.genLabels())

            ## these slices take args
            ## dimA, dimB, corner, color, planeString,
            ## minDesignCorner, maxDesignCorner, handleEnds

            xzSlice = Slice(x, z, frames.rowCoords[1], 'green', 'xz', minXZ, maxXZ, handleEnds)
            xzSlice.calcPosition()
            newDoc.addElement(xzSlice.draw())
            newDoc.addElement(xzSlice.drawHandles())
            newDoc.addElement(xzSlice.genLabels())
            zySlice = Slice(z, y, frames.rowCoords[2], 'red', 'zy', minZY, maxZY, handleEnds)
            zySlice.calcPosition()
            newDoc.addElement(zySlice.draw())
            newDoc.addElement(zySlice.drawHandles())
            newDoc.addElement(zySlice.genLabels())
            # iterate over the angles
            for i in range(0, len(angles)):
                isoView = View(x,y,z, frames.rowCoords[3 + i], float(angles[i]))
                isoView.CalcStart()
                isoView.makeFaces()
                newDoc.addElement(isoView.Render())
                #newDoc.addElement(isoView.genAxisLabels())
                #newDoc.addElement(isoView.drawPerpLabels())
                if showHelix:
                    newDoc.addElement(isoView.helixSlice(origamiRows, origamiCols))

            count += 1

    newDoc.save(svgFile)
    return 0


if __name__ == '__main__':
    genSVGFile(args.csvFileName, args.targetFileName, args.perpsViewCount, args.angleList, args.showHelix)
