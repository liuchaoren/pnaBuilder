__author__ = 'Chaoren'

import sys
import numpy as np
import pandas as pd
from database.pdbparser import *
from database.geoOperations import *
import copy

# from lna_lna_params import *
from pna_params import *


pdb2csv(atPDBfilename, atCSVfilename)
pdb2csv(taPDBfilename, taCSVfilename)
pdb2csv(gcPDBfilename, gcCSVfilename)
pdb2csv(cgPDBfilename, cgCSVfilename)


# import database
atframe = pd.read_csv(atCSVfilename, dtype = dtypes)
taframe = pd.read_csv(taCSVfilename, dtype = dtypes)
gcframe = pd.read_csv(gcCSVfilename, dtype = dtypes)
cgframe = pd.read_csv(cgCSVfilename, dtype = dtypes)

frameLookupTable = {}
frameLookupTable['A'] = atframe
frameLookupTable['T'] = taframe
frameLookupTable['G'] = gcframe
frameLookupTable['C'] = cgframe

class bpsSeq():
    # contains the information of a sequence
    def __init__(self, s):
        self.name = s
        self.frame = frameLookupTable[s].copy()
        twoPointsatTopLayer = np.array(self.frame.ix[twoPointsLookupTable[s], ['x', 'y', 'z']])
        self.topBasePoint = twoPointsatTopLayer[0]
        self.bottomBasePoint = self.topBasePoint
        self.topVector = twoPointsatTopLayer[1] - twoPointsatTopLayer[0]
        self.bottomVector = self.topVector
        self.topz = self.topBasePoint[-1]
        self.bottomz = self.topz


    def __add__(self, other):
        resultbpsSeq = copy.copy(self)
        othercopy = copy.copy(other)
        othercopy.moveFrameby(np.array([0., 0., self.topz - othercopy.bottomz + rise])) # move other in z direction
        # print othercopy.frame
        angleBetweenLayers = angleBetweenvectors(resultbpsSeq.topVector, othercopy.bottomVector)
        z = np.array([0., 0., 1.])
        if np.cross(self.topVector, othercopy.bottomVector)[-1] < 0:
            angleBetweenLayers = -1 * angleBetweenLayers
        othercopy.rotateFrameby(z, twistAngle - angleBetweenLayers) # rotation other around z

        rotation_matrix = rotationAround(z, slideAngle) # calculate the slides
        # print resultbpsSeq.frame
        # print othercopy.frame
        vectorAfterRotation = np.dot(self.topVector, np.transpose(rotation_matrix))
        slide = slideLen * vectorAfterRotation / np.linalg.norm(vectorAfterRotation)
        slidexy = -1 * othercopy.bottomBasePoint + self.topBasePoint + slide
        slidexy[-1] = 0.
        othercopy.moveFrameby(slidexy) # slide
        othercopy.frame.ix[:, ['x', 'y', 'z']] = np.round(np.array(othercopy.frame.ix[:, ['x', 'y', 'z']]), decimals=3)  # correct the xyz float value decimals
        bpsSeq.correctFrameIndexResid(resultbpsSeq, othercopy)
        resultbpsSeq.frame = resultbpsSeq.frame.append(othercopy.frame, verify_integrity=True)
        resultbpsSeq.frame = resultbpsSeq.frame.sort()
        resultbpsSeq.topBasePoint = othercopy.topBasePoint
        resultbpsSeq.topVector = othercopy.topVector
        resultbpsSeq.name = self.name + othercopy.name
        resultbpsSeq.topz = self.topz + (othercopy.topz - othercopy.bottomz) + rise
        return resultbpsSeq


    def moveFrameby(self, u):
        # move frame by u
        oldxyz = np.array(self.frame.ix[:, ['x', 'y', 'z']])
        newxyz = oldxyz + u
        self.frame.ix[:, ['x', 'y', 'z']] = newxyz
        self.topz = self.topz + u[-1]
        self.bottomz = self.bottomz + u[-1]
        self.topBasePoint = self.topBasePoint + u
        self.bottomBasePoint = self.topBasePoint + u

    def rotateFrameby(self, u, theta):
        # rotate around z axis
        rotation_matrix = rotationAround(u, theta)
        framexyz_old = np.array(self.frame.ix[:, ['x', 'y', 'z']])
        self.frame.ix[:, ['x', 'y', 'z']] = np.round(np.dot(framexyz_old, np.transpose(rotation_matrix)), decimals=3)
        self.topVector = np.dot(self.topVector, np.transpose(rotation_matrix))
        self.bottomVector = np.dot(self.bottomVector, np.transpose(rotation_matrix))
        self.topBasePoint = np.dot(self.topBasePoint, np.transpose(rotation_matrix))
        self.bottomBasePoint = np.dot(self.bottomBasePoint, np.transpose(rotation_matrix))


    def centerize(self):
        # put the center the molecule at (0, 0, 0)
        coors = self.frame.ix[:, ['x', 'y', 'z']]
        center = np.round(np.mean(coors, axis=0), decimals=3)
        coors = coors - center
        self.frame.ix[:, ['x', 'y', 'z']] = coors


    @staticmethod
    def correctFrameIndexResid(bpsSeq1, bpsSeq2):
        seq1 = bpsSeq1.name
        seq2 = bpsSeq2.name
        index1 = bpsSeq1.frame.index
        index2 = bpsSeq2.frame.index
        resid1 = bpsSeq1.frame['resSeq']
        resid2 = bpsSeq2.frame['resSeq']
        atomNum2 = len(bpsSeq2.frame)
        resNum1 = len(seq1)
        resNum2 = len(seq2)
        insertionPoint = np.sum(np.array([baseLenLookupTable[s] for s in seq1]))
        index1 = index1[:insertionPoint].append(index1[insertionPoint:] + atomNum2)
        index2 = index2 + insertionPoint
        resid1 = resid1[:insertionPoint].append(resid1[insertionPoint:] + resNum2 * 2)
        resid2 = resid2 + resNum1
        bpsSeq1.frame['resSeq'] = resid1
        bpsSeq2.frame['resSeq'] = resid2
        bpsSeq1.frame.index = index1
        bpsSeq2.frame.index = index2
        bpsSeq1.frame['atom_serial_number'] = bpsSeq1.frame.index + 1
        bpsSeq2.frame['atom_serial_number'] = bpsSeq2.frame.index + 1


for eachBP in sequence:
    if "sequenceFrame" in globals():
        sequenceFrame = sequenceFrame + bpsSeq(eachBP)
    else:
        sequenceFrame = bpsSeq(eachBP)

sequenceFrame.centerize()  #



print sequenceFrame.frame
# print bpsSeq("A").topBasePoint
# output file
csv2pdb(sequenceFrame.frame, sequenceFilename)
