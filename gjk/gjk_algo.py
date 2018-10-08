import numpy as np
from collections import deque


def getMinkowskyDiff(polygon1, polygon2):
    minkowsky_diff = []
    for p1 in polygon1:
        for p2 in polygon2:
            minkowsky_diff.append(np.array(p1) - np.array(p2))
    return minkowsky_diff



def cross(a, b):
    if isinstance(a, list): a = np.array(a)
    if isinstance(b, list): b = np.array(b)
    if a.shape[0] == 2:
        a = np.array([a[0], a[1], 0])
    if b.shape[0] == 2:
        b = np.array([b[0], b[1], 0])
    return np.cross(a, b)


def tripleCross(a, b, c):
    crs = cross(cross(a, b), c)
    return np.array([crs[0], crs[1]])


class Sphere(object):
    def __init__(self, center, radius):
        self.center = np.array(center)
        self.radius = radius

class Circle(object):
    def __init__(self, center, radius):
        self.center = np.array(center)
        self.radius = radius


def getFarthestPointInDirection(shape, d):
    d_norm = d / np.linalg.norm(d)
    if isinstance(shape, Sphere):
        return shape.center + d_norm * shape.radius
    if isinstance(shape, Circle):
        return shape.center + d_norm * shape.radius

    maxP = None
    maxVal = None
    for p in shape:
        dp = np.dot(np.array(p), d)
        if maxVal is None or dp > maxVal:
            maxVal = dp
            maxP = p
    return np.array(maxP)


def support(shape1, shape2, d):
    d = np.array(d)
    p1 = getFarthestPointInDirection(shape1, d)
    p2 = getFarthestPointInDirection(shape2, -d)
    return p1 - p2


def containsOrigin(shape):
    n = len(shape)
    sign = None
    for i in xrange(n):
        a = np.array(shape[i])
        b = np.array(shape[i+1] if i+1 < n else shape[0])
        if sign is None:
            sign = cross(b - a, -a)[2]  # 0 - a
        elif sign * cross(b - a, -a)[2] < 0:
            return False
    return True

def getEdge(tri):
    if len(tri) == 2:
        return tri
    
    cases = [
        (tri[0], tri[1], tri[2]),
        (tri[1], tri[2], tri[0]),
        (tri[2], tri[0], tri[1]),
    ]
    
    edge = None
    for case in cases:
        dedge = cross(case[1] - case[0], case[2] - case[1])[2]
        dorig = cross(case[1] - case[0], -case[1])[2]
        if dedge * dorig < 0:
            edge = (case[0], case[1])
            break

    return edge


class GJK(object):
    def __init__(self, shape1, shape2, init_d=(-1, 0)):
        self.shape1 = shape1
        self.shape2 = shape2
        self.simplex = [support(shape1, shape2, init_d)]
        self.d = -1.0 * np.array(init_d)
    
    def step(self):
        self.simplex.append(support(self.shape1, self.shape2, self.d))
        if np.dot(self.simplex[-1], self.d) <= 0:
            return "No overlap"

        if len(self.simplex) >= 3:
            newTri = self.simplex[-3:]
            if containsOrigin(newTri):
                return "Overlap"
        
        a, b = getEdge(self.simplex)
        self.simplex = [a, b]
        self.d = tripleCross(b-a, -a, b-a)
        return "In progress"



## --------- 3D Version ----------- ##

def containsOrigin3(shape):
    # only support tetrahedron
    if isinstance(shape, list):
        assert(len(shape) == 4)
    
    cases = [
        ([0, 1, 2], 3),
        ([0, 1, 3], 2),
        ([0, 2, 3], 1),
        ([1, 2, 3], 0)
    ]

    for case in cases:
        # get innerNormal
        face = [shape[i] for i in case[0]]
        vert = shape[case[1]]
        faceNormal = getFaceNormal(face) 
        edge = vert - face[0]
        if np.dot(faceNormal, edge) < 0:
            faceNormal = -faceNormal

        # check origin inside
        if np.dot(faceNormal, -face[0]) < 0:
            return False

    return True


def getFaceNormal(tri):
    return np.cross(tri[1] - tri[0], tri[2] - tri[1])


def getFace(tet):
    if len(tet) == 3:
        faceNormal = getFaceNormal(tet)
        if np.dot(faceNormal, -tet[0]) < 0:
            faceNormal = -faceNormal
        return tet, faceNormal

    cases = [
        ([0, 1, 2], 3),
        ([0, 1, 3], 2),
        ([0, 2, 3], 1),
        ([1, 2, 3], 0)
    ]

    tri, normal = None, None
    for case in cases:
        # get outerNormal
        face = [tet[i] for i in case[0]]
        vert = tet[case[1]]
        faceNormal = getFaceNormal(face)
        edge = vert - face[0]
        if np.dot(faceNormal, edge) > 0:
            faceNormal = -faceNormal
        
        # check origin direction
        if np.dot(faceNormal, -face[0]) > 0:
            tri = face
            normal = faceNormal
            break

    return tri, normal


class GJK3(object):
    def __init__(self, shape1, shape2, init_d=np.array([-1, 0, 0]), init_d2=np.array([0, -1, 0])):
        self.shape1 = shape1
        self.shape2 = shape2
        self.simplex = [support(shape1, shape2, init_d), support(shape1, shape2, init_d2)]
        self.d = -np.array(init_d)
    
    def step(self):
        self.simplex.append(support(self.shape1, self.shape2, self.d))
        if np.dot(self.simplex[-1], self.d) <= 0:
            return "No overlap"

        if len(self.simplex) >= 4:
            newTet = self.simplex[-4:]
            if containsOrigin3(newTet):
                return "Overlap"

        (a, b, c), self.d = getFace(self.simplex)
        self.simplex = [a, b, c]
        return "In progress"