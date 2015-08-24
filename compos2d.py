#!/usr/bin/env python

import Triangulate
import ellipt2d
import DirichletBound
import cell
import superlu
import copy
import math

import vtk
import numpy

class Compos2d:

    def __init__(self):
        """
        Constructor
        """
    
        self.refData = {}
        self.spcData = {}
    
        self.refDic = {}
        self.spcDic = {}

    def setReference(self, xyPoints):
        """
        Set reference coordinates
        @param xyPoints list of unique (x, y) tuples, ordered counterclockwise
        """
        
        refDic = self._buildSparseSystem(xyPoints)
        refDic['rho'] = superlu.solve(refDic['rhoMat'],
                                      refDic['rhoBVec'])
        refDic['the'] = superlu.solve(refDic['theMat'],
                                      refDic['theBVec'])
        self._toPolyData(refDic, self.refData)
        self.refDic = refDic
    
    def setSpecimen(self, xyPoints):
        """
        Set specimen coorindates
        @param xyPoints list of unique (x, y) tuples, ordered counterclockwise
        """
    
        spcDic = self._buildSparseSystem(xyPoints)
        spcDic['rho'] = superlu.solve(spcDic['rhoMat'],
                                      spcDic['rhoBVec'])
        spcDic['the'] = superlu.solve(spcDic['theMat'],
                                      spcDic['theBVec'])
        self._toPolyData(spcDic, self.spcData)
        self.spcDic = spcDic

    def findSpecimenPoint(self, referencePoint, h=0.01, niter=10, tol=1.e-6):
        """
        Find the x, y position in the specimen corresponding to a given position
        in the reference
        @param referencePoint point in the reference object
        @param h small excursion used to compute the finite difference Jacobian
        @param niter max number of Newton iterations
        @param tol max error in the fields
        @return specimen position, error, and number of iterations
        """
        
        mat = numpy.zeros((2,2), numpy.float64)
    
        # compute field values on the reference object
        refStarData = self.starInterpolation(self.refData, referencePoint, h)
        refField = self._averageStarData(refStarData)

        # initial guess
        spcPos = numpy.array(referencePoint)
        error = float('inf')
        iter = 0
        while error > tol and iter <= niter:
    
            # compute the Jacobian at the specimen location
            spcStarData = self.starInterpolation(self.spcData, spcPos, h)
            spcField = self._averageStarData(spcStarData)
            print '.... spcStarData = ', spcStarData, ' spcPos = ', spcPos, ' h = ', h
            mat[0, 0] = (spcStarData['e'][0] - spcStarData['w'][0])/(2. * h) # d rho / dx
            mat[0, 1] = (spcStarData['n'][0] - spcStarData['s'][0])/(2. * h) # d rho / dy
            mat[1, 0] = (spcStarData['e'][1] - spcStarData['w'][1])/(2. * h) # d the / dx
            mat[1, 1] = (spcStarData['n'][1] - spcStarData['s'][1])/(2. * h) # d the / dy
        
            # make sure the determinant is non-zero
            det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
            if det == 0:
                print '*** refStarData = ', refStarData
                print '*** spcStarData = ', spcStarData
                print '*** mat = ', mat
                print '*** cannot invert mat! zero determinant'
                break
            matInv = numpy.linalg.inv(mat)
        
            # solve for the correction on the specimen
            deltas = numpy.dot(matInv, refField - spcField)
            newSpcPos = spcPos + deltas
            
            # update the specimen position
            newSpcStarData = self.starInterpolation(self.spcData, newSpcPos, h)
            spcField = self._averageStarData(newSpcStarData)
            
            # update the error
            newError = math.sqrt(numpy.dot(refField - spcField, refField - spcField))
    
            if newError < error:
                # accept new position
                spcPos = newSpcPos
                error = newError
            else:
                # backtrack and try again
                spcPos += 0.3 * (newSpcPos - spcPos)
            
            iter += 1
    
        return spcPos, error, iter
    
    def _averageStarData(self, starData):
        """
        Compute the average field from the star stencil
        @param starData dictionary containing est, west, north, and south field values
        @return field
        """
        field = numpy.zeros((2,), numpy.float64)
        for k, f in starData.items():
            field += numpy.array(f)
        field /= float(len(starData))
        return field
    
    def _toPolyData(self, xDic, xData):
        """
        Convert elltip2d data to VTK vtkPolyData object
        @param xDic dictionalry of ellipt2d objects
        @param xData dictionary of VTK objects (ouput)
        """
        grid = xDic['grid']
        cells = xDic['cells']
        numPoints = len(grid)
        numCells = len(cells.data)
        
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(numPoints)
        for i in range(numPoints):
            x, y = grid[i][0]
            points.SetPoint(i, x, y, 0.0)
        
        pdata = vtk.vtkPolyData()
        pdata.SetPoints(points)
        pdata.Allocate(numCells, 1)
        ptIds = vtk.vtkIdList()
        npts = 3 # triangles
        ptIds.SetNumberOfIds(3)
        for i in range(numCells):
            cell = cells.data[i]
            for j in range(npts):
                ptIds.SetId(j, cell[j])
            pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)
        
        pointData = pdata.GetPointData()
        numComps = 2 # 2 fields
        vArray = vtk.vtkDoubleArray()
        vArray.SetNumberOfComponents(2) # rho and the
        vArray.SetNumberOfTuples(numPoints)
        rhoArray, theArray = xDic['rho'], xDic['the']
        for i in range(numPoints):
            vals = rhoArray[i], theArray[i]
            #print '.... {} setting rho and the: {}'.format(i, vals)
            vArray.SetTuple(i, vals)
        pointData.SetScalars(vArray)

        # store
        xData['points'] = points
        xData['polydata'] = pdata
        xData['array'] = vArray

    def starInterpolation(self, data, xy, h):
        """
        Interpolate along a star stencil 
        @param data either self.refData or self.spcData
        @param xy x and y coordinates at center of stencil
        @param h excursion from the center
        @return {'w': value, 'e': value, 'n': value, 's': value}
        """
        
        lineX = vtk.vtkLineSource()
        lineX.SetPoint1(xy[0] - h, xy[1], 0.0)
        lineX.SetPoint2(xy[0] + h, xy[1], 0.0)
        lineX.SetResolution(1)

        lineY = vtk.vtkLineSource()
        lineY.SetPoint1(xy[0], xy[1] - h, 0.0)
        lineY.SetPoint2(xy[0], xy[1] + h, 0.0)
        lineY.SetResolution(1)

        probeX = vtk.vtkProbeFilter()
        if vtk.VTK_MAJOR_VERSION >= 6:
            probeX.SetSourceData(data['polydata'])
        else:
            probeX.SetSource(data['polydata'])
        probeX.SetInputConnection(lineX.GetOutputPort())
        probeX.Update()

        probeY = vtk.vtkProbeFilter()
        if vtk.VTK_MAJOR_VERSION >= 6:
            probeY.SetSourceData(data['polydata'])
        else:
            probeY.SetSource(data['polydata'])
        probeY.SetInputConnection(lineY.GetOutputPort())
        probeY.Update()
        
        res = {}
        
        # west and east
        res['w'] = probeX.GetOutput().GetPointData().GetArray(0).GetTuple(0)
        res['e'] = probeX.GetOutput().GetPointData().GetArray(0).GetTuple(1)
        
        # south and north
        res['s'] = probeY.GetOutput().GetPointData().GetArray(0).GetTuple(0)
        res['n'] = probeY.GetOutput().GetPointData().GetArray(0).GetTuple(1)
        
        return res

    def _triangulatePoints(self, xyPoints):
        """
        Triangulate a set of points
        @param xyPoints list of unique (x, y) tuples, ordered counterclockwise
        @return node object and attribute values
        """
        
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
        """
        Build sparse matrix system
        @param xyPoints list of unique (x, y) tuples, ordered counterclockwise
        @return ellipt2d objects in a dictionary
        """
    
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
        """
        Plot solution 
        @param refFlag set to True if reference solution, False for specimen
        @param var either 'rho' or 'the'
        """
        
        import tkplot
        from Tkinter import Tk, Frame, Canvas
        root = Tk()
        frame = Frame(root)
        frame.pack()
        width, height = 500, 450
        canvas = Canvas(bg="white", width=width, height=height)
        canvas.pack()
        title = 'reference: ' + var
        data = self.refDic
        if not refFlag:
            'specimen: ' + var
            data = self.spcDic
        tkplot.tkplot(canvas, data['grid'], data[var], 0,0,1,
                            title=title, WIDTH=width, HEIGHT=height)
        root.mainloop()

#####################################################

def test1():
    
    cs = Compos2d()
    
    pts = [(1., 0.5), (1., 1.), (0., 1.), (0., 0.), (1., 0.)]
    cs.setReference(pts)
    print cs.starInterpolation(cs.refData, xy=(0.5, 0.7), h=0.01)

    # add perturbations
    n = len(pts)
    spcPts = [ (pts[i][0] + 0.10*math.sin(2*i*2*math.pi/float(n)),
                pts[i][1] + 0.05*math.cos(2*i*2*math.pi/float(n)))
                for i in range(n)]
    cs.setSpecimen(spcPts)
    print cs.starInterpolation(cs.spcData, xy=(0.5, 0.7), h=0.01)
    
    refPos = (0.3, 0.4)
    spcPos, error, iter = cs.findSpecimenPoint(refPos, tol=1.e-6, niter=10, h=0.01)
    print 'refPos = {} spcPos = {} error = {} iter = {}'.format(refPos, spcPos, error, iter)
    
    #cs.plot(refFlag=True, var='rho')
    #cs.plot(refFlag=False, var='the')

if __name__ == '__main__':
    test1()
