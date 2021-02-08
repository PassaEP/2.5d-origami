from pysvg.structure import G
from pysvg.builders import ShapeBuilder, TransformBuilder, StyleBuilder
from pysvg.text import *
from pysvg.shape import *
import math

FRAMELEN = 200

class RowFrame:
    def __init__(self, i, fpr):
        self.i = i
        self.rowCoords = []
        self.fpr = fpr
    def makerow(self):
        frameStyle = StyleBuilder()
        frameStyle.setStrokeWidth(1.5)
        frameStyle.setStroke('black')
        frameStyle.setFilling('white')
        out = G()
        yCoord = 10 + (FRAMELEN + 20)*self.i
        for z in range(0, self.fpr):
            newFrame = Rect(10 + (FRAMELEN + 20)*z, yCoord, FRAMELEN,FRAMELEN)
            newFrame.set_style(frameStyle.getStyle())
            out.addElement(newFrame)
            self.rowCoords.append([10 + (FRAMELEN + 20)*z, yCoord])

        return out
