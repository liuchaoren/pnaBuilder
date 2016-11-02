__author__ = 'Chaoren'
import numpy as np


def rotationAround(u, theta):
    # rotate matrix: rotate theta degree around u, u is a unit vector
    c = np.cos(theta)
    s = np.sin(theta)
    ux = u[0]
    uy = u[1]
    uz = u[2]
    return np.array([[c + ux**2 * (1 - c), ux * uy * (1 - c) - uz * s, ux * uz * (1 - c) + uy * s], [uy * ux * (1 - c) + uz * s, c + uy**2 * (1 - c), uy * uz * (1 - c) - ux * s],
            [uz * ux * (1 - c) - uy * s, uz * uy * (1 - c) + ux * s, c + uz ** 2 * (1 - c)]])


def planeUnitVector(threePointsFrame):
    # given a 3x3 frame with each line being a point, return a unit vector perpendicular to this plane
    threePointsArray = np.array(threePointsFrame)
    vector1 = threePointsArray[1] - threePointsArray[0]
    vector2 = threePointsArray[2] - threePointsArray[0]
    crossProduct = np.cross(vector1, vector2)
    unitvector = crossProduct / np.linalg.norm(crossProduct)
    return unitvector

def angleBetweenvectors(vector1, vector2):
    # angle between two vectors, vector1 and vector can be non-normalized
    costheta = np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
    if costheta > 1:
        costheta = 1
    elif costheta < -1:
        costheta = -1
    return np.arccos(costheta)


def decomposeVector(a, u):
    # give vector a and u (u is unit length), return two vectors, the first one is part of a in u direction and the second one is part perpendicular to u
    paralellPart = np.dot(a,u) * u
    return (paralellPart, a - paralellPart)     # (parallel part, perpendicular part)
