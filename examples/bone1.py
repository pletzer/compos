#!/usr/bin/en python


xy = [
      (0., 0.),
      (-0.2, 1.),
      (-1., 2.5),
      (-0.7, 2.8),
      (0.3, 3.2),
      (1.7, 3.),
      (2., 2.5),
      (2., 1.7),
      (1.7, 1.3),
      (1.3, 0.8),
      (1., 0.),
     ]

if __name__ == '__main__':
    from matplotlib import pylab
    xs = [p[0] for p in xy] + [xy[0][0]]
    ys = [p[1] for p in xy] + [xy[0][1]]
    pylab.plot(xs, ys)
    pylab.show()


