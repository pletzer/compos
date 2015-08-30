#!/usr/bin/env python

import compos2d
import sys
import math


cs = compos2d.Compos2d()
sys.path.append('examples')

# input data
import bone1
pts = bone1.xy
cs.setReference(pts)

# modify reference
n = len(pts)
spcPts = [ (pts[i][0] + 0.10*math.sin(2*i*2*math.pi/float(n)),
            pts[i][1] + 0.05*math.cos(2*i*2*math.pi/float(n)))
                for i in range(n)]
cs.setSpecimen(spcPts)

refPos = (1.5, 1.8)
spcPos, error, iter, history = cs.findSpecimenPoint(refPos, tol=1.e-6, niter=10, h=0.01)
print 'refPos = {} spcPos = {} error = {} iter = {}'.format(refPos, spcPos, error, iter)

cs.plot(refFlag=True, title="reference", xyPath=[refPos])
cs.plot(refFlag=False, title="specimen", xyPath=history)
