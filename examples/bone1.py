#!/usr/bin/en python


xy = [
      (0.4, 0.),
      (0.5, 0.3),
      (0.8, 0.8),
      (1., 1.2),
      (0.9, 1.4),
      (0.8, 1.5),
      (0.5, 1.7),
      (0.2, 1.8),
      (-0.3, 1.8),
      (-0.5, 1.5),
      (-0.7, 1.3),
      (-0.85, 1.2),
      (-0.83, 1.0),
      (-0.8, 0.7),
      (-0.3, 0.),
     ]

if __name__ == '__main__':
    from matplotlib import pylab
    xs = [p[0] for p in xy] + [xy[0][0]]
    ys = [p[1] for p in xy] + [xy[0][1]]
    pylab.plot(xs, ys)
    pylab.show()


