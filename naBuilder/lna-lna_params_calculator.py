__author__ = 'Chaoren'

import pandas as pd
from na_params_funcs import *
from database.pdbparser import *
import numpy as np

pdbfilename = "lna.pdb"
csvfilename = "lna.csv"

layerindex = ("layer1", "layer2", "layer3", "layer4", "layer5", "layer6", "layer7")
# layerPoints has 6 elements, each element has three ints, ints are the atom index of C8' C8' and an optional one
layerPoints = ((165, 134, 163), (188, 112, 186), (204, 100, 210), (234, 67, 232), (258, 48, 256), (281, 23, 279), (305, 4, 303))
layerInfo = zip(layerindex, layerPoints)


pdb2csv(pdbfilename, csvfilename)
pdbframe = pd.read_csv(csvfilename)
unitVectorPlane = unitVectorCalculator(pdbframe, layerInfo)
planeRise = riseCalculator(pdbframe, layerInfo, unitVectorPlane)
angle = angle(pdbframe, layerInfo)
slide = slide(pdbframe, layerInfo, unitVectorPlane)

# print unitVectorPlane
print planeRise
print angle
print slide
# 1 represents two parallel plane
# print [np.dot(i,j) for (i,j) in zip(unitVectorPlanes[1:], unitVectorPlanes[:-1])]

