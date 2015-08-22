#!/usr/bin/env python

import Triangulate
import ellipt2d
import DirichletBound
import cell
import superlu
import copy
import math

class Compos2d:

    def __init__(self):
    
        self.refData = {}
        self.spcData = {}

    def setReference(self, xyPoints):
        
        self.refData = self._buildSparseSystem(xyPoints)
        self.refData['rho'] = superlu.solve(self.refData['rhoMat'],
                                            self.refData['rhoBVec'])
        self.refData['the'] = superlu.solve(self.refData['theMat'],
                                            self.refData['theBVec'])
    
    def setSpecimen(self, xyPoints):
    
        self.spcData = self._buildSparseSystem(xyPoints)
        self.spcData['rho'] = superlu.solve(self.spcData['rhoMat'],
                                            self.spcData['rhoBVec'])
        self.spcData['the'] = superlu.solve(self.spcData['theMat'],
                                            self.spcData['theBVec'])

    def findSpecimenPoint(self, referencePoint):
    
        rhoRef = self._interp(self.refData, referencePoint, 'rho')
        theRef = self._interp(self.refData, referencePoint, 'the')
   
    def _interp(self, data, xy, var):
        
        cells = data['cells']
        v = data[var]
        return cells.interp(v, xy[0], xy[1])

    def _triangulatePoints(self, xyPoints):
        
        numPoints = len(xyPoints)

        xs = [p[0] for p in xyPoints]
        ys = [p[1] for p in xyPoints]
        xmin = min(xs)
        xmax = max(xs)
        ymin = min(ys)
        ymax = max(ys)
        
        # max cell area
        areaMax = (xmax - xmin)*(ymax - ymin) / 1000.0
                    
        # mid point
        xmid = reduce(lambda x, y: x + y, xs) / float(numPoints)
        ymid = reduce(lambda x, y: x + y, ys) / float(numPoints)
        
        # boundary points. contour starts in the middle, goes counterclockwise
        # around. To close the loop we need to add points that are very close
        # to (xmid, ymid) and xyPoints[-1]
        
        eps = 1.2456e-1 * math.sqrt(areaMax)
        
        pts = [(xmid, ymid+eps), (xyPoints[0][0], xyPoints[0][1]+eps)] + \
               xyPoints[1:] + \
              [(xyPoints[0][0], xyPoints[0][1]-eps), (xmid, ymid-eps)]
        numPts = len(pts)
        
        # to indicate that the points are boundary points. we will need
        # to keep track of boundary points after triangulation
        mrkrs = [1 for i in range(numPts)]
        
        segs = [(i, i+1) for i in range(numPts-1)] + [(numPts-1, 0)]
        
        # pass the BCs as attributes so Steiner points can interpolate
        x0, y0 = xyPoints[0]
        the0 = 0.5 * math.atan2(y0 - ymid, x0 - xmid) / math.pi
        def getTheta(p):
            # the is a poloidal-like variable in the range 0 to 1
            the = 0.5 * math.atan2(p[1] - ymid, p[0] - xmid) / math.pi
            the -= the0
            if the < 0:
                the += 1.0
            return the
        
        attrs = [(0., 0.)] + \
                [(1., getTheta(p)) for p in xyPoints] + \
                [(1., 1.), (0., 1.)]
        #print attrs
            
        # triangulate
        tri = Triangulate.Triangulate()
        tri.set_points(pts, mrkrs)
        tri.set_segments(segs)
        tri.set_attributes(attrs)
        tri.triangulate(area=areaMax)
        
        attrs = tri.get_attributes()
        grid = tri.get_nodes()
        
        return grid, attrs

    def _buildSparseSystem(self, xyPoints):
    
        grid, bcs = self._triangulatePoints(xyPoints)
        #grid.plot()
        numPoints = len(xyPoints)
    
        # Laplace equation solver
        eq = ellipt2d.ellipt2d(grid, '1.0', '0.0', '0.0')
        amat, bvec = eq.stiffnessMat()
    
        #
        # rho field
        #
    
        rhoMat = copy.deepcopy(amat)
        rhoBVec = copy.deepcopy(bvec)
    
        bcData = {}
        for i, g in grid.data.items():
            if g[2]:
                # boundary
                #print 'x, y = {} rho Dirichlet BC is {}'.format(grid[i][0], bcs[i][0])
                bcData[i] = bcs[i][0]
        bc = DirichletBound.DirichletBound(bcData)
        # modify the matrix and source vector
        eq.dirichletB(bc, rhoMat, rhoBVec)
    
        #
        # theta field
        #
    
        theMat = copy.deepcopy(amat)
        theBVec = copy.deepcopy(bvec)
    
        bcData = {}
        for i, g in grid.data.items():
            if g[2]:
                # boundary
                #print 'x, y = {} theta Dirichlet BC is {}'.format(grid[i][0], bcs[i][1])
                bcData[i] = bcs[i][1]
        bc = DirichletBound.DirichletBound(bcData)
        eq.dirichletB(bc, theMat, theBVec)

        return {'eq': eq, 'grid': grid, 'cells': cell.cell(grid),
            'rhoMat': rhoMat, 'rhoBVec': rhoBVec,
            'theMat': theMat, 'theBVec': theBVec}

    def plot(self, refFlag=True, var='rho'):
        
        import tkplot
        from Tkinter import Tk, Frame, Canvas
        root = Tk()
        frame = Frame(root)
        frame.pack()
        width, height = 500, 450
        canvas = Canvas(bg="white", width=width, height=height)
        canvas.pack()
        title = 'reference: ' + var
        data = self.refData
        if not refFlag:
            'specimen: ' + var
            data = self.spcData
        tkplot.tkplot(canvas, data['grid'], data[var], 0,0,1,
                            title=title, WIDTH=width, HEIGHT=height)
        root.mainloop()

#####################################################

def test1():
    cs = Compos2d()
    pts = [(1., 0.5), (1., 1.), (0., 1.), (0., 0.), (1., 0.)]
    cs.setReference(pts)
    cs.plot(refFlag=True, var='rho')
    cs.plot(refFlag=True, var='the')

if __name__ == '__main__':
    test1()
