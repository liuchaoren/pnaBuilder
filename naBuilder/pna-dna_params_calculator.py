__author__ = 'Chaoren'

import pandas as pd
from na_params_funcs import *
from database.pdbparser import *
import numpy as np

pdbfilename = "1pdt-noh.pdb"
csvfilename = "1pdt.csv"

layerindex = ("layer1", "layer2", "layer3", "layer4", "layer5", "layer6", "layer7", "layer8")
# layerPoints has 6 elements, each element has three ints, ints are the atom index of C8' C8' and an optional one
layerPoints = ((7, 299, 13), (29, 280, 35), (50, 259, 56), (69, 240, 75), (90, 220, 96), (110, 201, 116), (131, 183, 137), (153, 162, 159))
layerInfo = zip(layerindex, layerPoints)


# pdb2csv(pdbfilename, csvfilename)
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

