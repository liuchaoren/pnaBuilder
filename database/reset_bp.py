__author__ = 'Chaoren'

import pandas as pd
from geoOperations import  *
from pdbparser import *
import numpy as np

# constants
# pdbfilename = "gc_original.pdb"
# csvfilename = "pna-gc.csv"
# newpdbfilename1 = "pna-gc.pdb"
# newpdbfilename2 = "pna-cg.pdb"
# planeAtomIndexes = (0, 17, 18)     # index from 0
#
# pdbfilename = "at_original.pdb"
# csvfilename = "at_original.csv"
# newpdbfilename1 = "pna-ta.pdb"
# newpdbfilename2 = "at_original.pdb"
# planeAtomIndexes = (0, 15, 20)     # index from 0

#
# pdbfilename = "pna-dna-gc.pdb"
# csvfilename = "pna-dna-gc.csv"
# newpdbfilename1 = "pna-dna-gc-1.pdb"
# newpdbfilename2 = "pna-dna-gc-2.pdb"
# planeAtomIndexes = (10, 22, 16)     # index from 0



pdbfilename = "lna-lna-at.pdb"
csvfilename = "pna-dna-at.csv"
newpdbfilename1 = "lna-lna-at-1.pdb"
newpdbfilename2 = "lna-lna-at-2.pdb"
planeAtomIndexes = (15, 30, 36)     # index from 0



# generate a csv file from a pdb file
pdb2csv(pdbfilename, csvfilename)

pdbframe = pd.read_csv(csvfilename)
bpPlane = pdbframe.ix[planeAtomIndexes, ['x', 'y', 'z']]
bpPlaneVector = planeUnitVector(bpPlane)

z = np.array([0, 0 ,1])
u = np.cross(bpPlaneVector, z)
unorm = np.linalg.norm(u)
u = u / unorm

theta = angleBetweenvectors(bpPlaneVector, z)

# rotation_matrix 1 and 2 are two matrixs: rotate basepair such that it is perpendicular to z and -z
rotation_matrix1 = rotationAround(u, theta)
rotation_matrix2 = rotationAround(u, theta - np.pi)

pdbframexyz_old = np.array(pdbframe.ix[:, ['x', 'y', 'z']])
pdbframexyz_new1 = np.round(np.dot(pdbframexyz_old, np.transpose(rotation_matrix1)), decimals=3)
pdbframexyz_new2 = np.round(np.dot(pdbframexyz_old, np.transpose(rotation_matrix2)), decimals=3)

newpdbframe1 = pdbframe.copy()
newpdbframe1.ix[:, ['x', 'y', 'z']] = pdbframexyz_new1
newpdbframe2 = pdbframe.copy()
newpdbframe2.ix[:, ['x', 'y', 'z']] = pdbframexyz_new2

zcoor1 = newpdbframe1.ix[planeAtomIndexes[0], 'z']
zcoor2 = newpdbframe2.ix[planeAtomIndexes[0], 'z']

newpdbframe1.ix[:, 'z'] = np.array(newpdbframe1.ix[:, 'z']) - zcoor1
newpdbframe2.ix[:, 'z'] = np.array(newpdbframe2.ix[:, 'z']) - zcoor2

csv2pdb(newpdbframe1, newpdbfilename1)
csv2pdb(newpdbframe2, newpdbfilename2)



