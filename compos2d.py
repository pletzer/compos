#!/usr/bin/env python

import Triangulate
import ellipt2d
import DirichletBound
import cell
import superlu
import copy
import math

import vtk

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

    def findSpecimenPoint(self, referencePoint):
    
        pass
        #rhoRef = self._interp(self.refDic, referencePoint, 'rho')
        #theRef = self._interp(self.refDic, referencePoint, 'the')
    
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
        Build sparse matri system 
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
            data = self.spcData
        tkplot.tkplot(canvas, data['grid'], data[var], 0,0,1,
                            title=title, WIDTH=width, HEIGHT=height)
        root.mainloop()

#####################################################

def test1():
    cs = Compos2d()
    pts = [(1., 0.5), (1., 1.), (0., 1.), (0., 0.), (1., 0.)]
    cs.setReference(pts)
    print cs.starInterpolation(cs.refData, xy=(0.5, 0.7), h=0.01)
    #cs.plot(refFlag=True, var='rho')
    #cs.plot(refFlag=True, var='the')

if __name__ == '__main__':
    test1()
