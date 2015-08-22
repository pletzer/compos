#!/usr/bin/env python

import Triangulate
import ellipt2d
import DirichletBound
import copy
import math

class Compos2d:

    def __init__(self):
    
        self.refEq = None
    
        self.eq = None

    def setReference(self, xyPoints):
    
        self.gridRef = self._triangulatePoints(xyPoints)
        numPoints = len(xyPoints)
        
        # Laplace equation solver
        self.refEq = ellipt2d.ellipt2d(self.gridRef, '1.0', '0.0', '0.0')
        amat, bvec = self.refEq.stiffnessMat()
        
        #
        # rho field
        #
        
        self.rhoRefMat = copy.deepcopy(amat)
        self.rhoRefBVec = copy.deepcopy(bvec)
        
        # zero Dirichlet on 0,
        # one Dirichlet on 1: numPoints + 2
        bcData = {i: 1.0 for i in range(1, numPoints + 2)}
        bcData[0] = 0.0
        bcData[numPoints + 1] = 0.0
        bc = DirichletBound.DirichletBound(bcData)
        self.refEq.dirichletB(bc, self.rhoRefMat, self.rhoRefBVec)
        
        #
        # theta field
        #
        
        self.theRefMat = copy.deepcopy(amat)
        self.theRefBVec = copy.deepcopy(bvec)
        
        # zero Dirichlet on 0,
        # one Dirichlet on 1: numPoints + 2
        bcData = {0: 0.0, 1: 0.0, numPoints: 1.0, numPoints+1: 1.0}
        bcData[numPoints + 1] = 1.0
        bc = DirichletBound.DirichletBound(bcData)
        self.refEq.dirichletB(bc, self.theRefMat, self.theRefBVec)

    def setSpecimen(self, xyPoints):
        pass

    def findSpecimenPoint(self, referencePoint):
        pass

    def _triangulatePoints(self, xyPoints):
        
        numPoints = len(xyPoints)

        xs = [p[0] for p in xyPoints]
        ys = [p[1] for p in xyPoints]
        xmin = min(xs)
        xmax = max(xs)
        ymin = min(ys)
        ymax = max(ys)
        
        # max cell area
        areaMax = (xmax - xmin)*(ymax - ymin)/100.0
                    
        # mid point
        xmid = reduce(lambda x, y: x + y, xs) / float(numPoints)
        ymid = reduce(lambda x, y: x + y, ys) / float(numPoints)
        
        # boundary points. contour starts in the middle, goes counterclockwise
        # around. To close the loop we need to add points that are very close
        # to (xmid, ymid) and xyPoints[-1]
        eps = 1.2456e-6*math.sqrt(areaMax)
        pts = [(xmid, ymid)] + xyPoints + \
              [(xyPoints[-1][0], xyPoints[-1][1]-eps)] + [(xmid, ymid-eps)]
        numPts = len(pts)
        segs = [(i, i+1) for i in range(numPts-1)] + [(numPts-1, 0)]
            
        # triangulate
        self.refTri = Triangulate.Triangulate()
        self.refTri.set_points(pts)
        self.refTri.set_segments(segs)
        self.refTri.triangulate(area=areaMax)

        return self.refTri.get_nodes()


#####################################################

def test1():
    cs = Compos2d()
    pts = [(1., 0.), (1., 1.), (0., 1.), (0., 0.)]
    cs.setReference(pts)

if __name__ == '__main__':
    test1()
