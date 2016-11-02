__author__ = 'Chaoren'

import pandas as pd
from na_params_funcs import *
from database.pdbparser import *
import numpy as np

pdbfilename = "1pup.pdb"
csvfilename = "1pup.csv"

layerindex = ("layer1", "layer2", "layer3", "layer4", "layer5", "layer6")
# layerPoints has 6 elements, each element has three ints, ints are the atom index of C8' C8' and an optional one
layerPoints = ((128, 106, 145), (146, 88, 161), (167, 63, 188), (191, 39, 211), (216, 18, 233), (234, 0, 249))
layerInfo = zip(layerindex, layerPoints)


# pdb2csv(pdbfilename, csvfilename)
pdbframe = pd.read_csv(csvfilename)
unitVectorPlane = unitVectorCalculator(pdbframe, layerInfo)
planeRise = riseCalculator(pdbframe, layerInfo, unitVectorPlane)
angle = angle(pdbframe, layerInfo)
slide = slide(pdbframe, layerInfo, unitVectorPlane)

print unitVectorPlane
print planeRise
print angle
print slide
# 1 represents two parallel plane
# print [np.dot(i,j) for (i,j) in zip(unitVectorPlanes[1:], unitVectorPlanes[:-1])]

