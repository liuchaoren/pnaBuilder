__author__ = 'Chaoren'

import pandas as pd
from database.geoOperations import *;
import numpy as np

def unitVectorCalculator(pdbframe, layerInfo):
    # calculate the average unit vector which is perpendicular to the bp planes
    unitVectorPlanes = []
    for eachlayer in layerInfo:
        (layername, layerPointsIndex) = eachlayer
        bpPlane = pdbframe.ix[layerPointsIndex, ['x', 'y', 'z']]
        bpPlaneVector = planeUnitVector(bpPlane)
        unitVectorPlanes.append(bpPlaneVector)
    # return unitVectorPlanes
    return np.mean(unitVectorPlanes, axis=0)


def riseCalculator(pdbframe, layerInfo, unitVectorPlane):
    # calculate the average bp layer rise
    riseLayers = []
    for eachlayer in layerInfo:
        (layername, layerPointsIndex) = eachlayer
        pointAtPlane = pdbframe.ix[layerPointsIndex[0], ['x', 'y', 'z']]
        if layername == "layer1":
            pointAtPlanePre = pointAtPlane
        else:
            rise = np.abs(np.dot(pointAtPlane - pointAtPlanePre, unitVectorPlane))
            riseLayers.append(rise)
            pointAtPlanePre = pointAtPlane
    return np.mean(np.array(riseLayers))
    # return riseLayers

def slide(pdbframe, layerInfo, unitVectorPlane):
    # return the shifting of the starting points of the two vectors in the two nearest neighbor bp plane
    slideLen = []
    slideAngle = []
    for eachlayer in layerInfo:
        (layername, layerPointsIndex) = eachlayer
        pointAtPlane = pdbframe.ix[layerPointsIndex[0], ['x', 'y', 'z']]
        pointsAtPlane = np.array(pdbframe.ix[layerPointsIndex[:2], ['x', 'y', 'z']])
        vectorAtPlane = pointsAtPlane[1] - pointsAtPlane[0]
        if layername == "layer1":
            pointAtPlanePre = pointAtPlane
            vectorAtPlanePre = vectorAtPlane
        else:
            shiftofVectorStartPoint = np.array(pointAtPlane - pointAtPlanePre)
            para, perpen = decomposeVector(shiftofVectorStartPoint, unitVectorPlane)
            slideLen.append(np.linalg.norm(perpen))
            slideAngle.append(angleBetweenvectors(perpen,vectorAtPlanePre))
        pointAtPlanePre = pointAtPlane
        vectorAtPlanePre = vectorAtPlane
    # return (np.array(slideLen), np.array(slideAngle))
    return (np.mean(np.array(slideLen)), np.mean(np.array(slideAngle)))


def angle(pdbframe, layerInfo):
    # calculate the rotation angle between nearest neighbor layers
    angles = []
    for eachlayer in layerInfo:
        (layername, layerPointsIndex) = eachlayer
        pointsAtPlane = np.array(pdbframe.ix[layerPointsIndex[:2], ['x', 'y', 'z']])
        vectorAtPlane = pointsAtPlane[0] - pointsAtPlane[1]
        if layername == "layer1":
            vectorAtPlanePre = vectorAtPlane
        else:
            angle = angleBetweenvectors(vectorAtPlane, vectorAtPlanePre)
            if angle > np.pi / 4:
                raise Exception
            else:
                angles.append(angle)
            vectorAtPlanePre = vectorAtPlane
    return np.mean(np.array(angles))
    # return np.array(angles)
