#! /usr/bin/env python

"""
Tkinter color plot utility (for real functions)
"""

try:
    from tkinter import Tk, Frame, Button, Canvas, BOTTOM # python 3
except:
    from Tkinter import Tk, Frame, Button, Canvas, BOTTOM # python 2
import cell
import math

def colormap(rho, the):
    r = rho
    g = the
    b = 1. - rho
    return "#%02x%02x%02x" % (r*255, g*255, b*255)

class Plot:

    def __init__(self, root, grid, title='', width=400, height=400):
        self.width = width
        self.height = height
        self.title = title
        self.grid = grid
    
        self.frame = Frame(root)
        self.frame.pack()
        self.canvas = Canvas(bg="white", width=width, height=height)
        self.canvas.pack()

    def draw(self, rho, the):
        
        border_x, border_y = 0.2, 0.2
        xmin, ymin, xmax, ymax = self.grid.boxsize()
        scale = min((1.-border_x)*self.width/(xmax-xmin), \
                    (1.-border_y)*self.height/(ymax-ymin))
        a = max(border_x * self.width/2., (self.width-scale*(xmax-xmin))/2.)
        b = max(border_y * self.height/2.,(self.height-scale*(ymax-ymin))/2.)
    
        box_pix = (a, self.height-scale*(ymax-ymin) - b,
                   scale*(xmax-xmin) + a, self.height - b)
        # draw the box
        self.addBox((xmin, ymin, xmax, ymax), box_pix)
    
        # color plot
        cells = cell.cell(self.grid)
        for index in range(len(cells.data)):
                ia, ib, ic = cells.data[index]
                xa, ya = scale*(self.grid.x(ia)-xmin) + a, \
                         self.height -scale*(self.grid.y(ia)-ymin) - b
                xb, yb = scale*(self.grid.x(ib)-xmin) + a, \
                         self.height -scale*(self.grid.y(ib)-ymin) - b
                xc, yc = scale*(self.grid.x(ic)-xmin) + a, \
                         self.height -scale*(self.grid.y(ic)-ymin) - b

                rhoAbc = (rho[ia] + rho[ib] + rho[ic])/3.0
                theAbc = (the[ia] + the[ib] + the[ic])/3.0
                color = colormap(rhoAbc, theAbc)
                self.canvas.create_polygon(xa, ya, xb, yb, xc, yc, fill=color)
    
        # add title
        x, y = self.width/3., self.height/15.0
        self.canvas.create_text(x, y, text=str(self.title),
                                font =("Helvetica", 14),
                                fill='black')


    def addBox(self, xxx_todo_changeme, xxx_todo_changeme1):
        """
        Add Box. Min/Max coordinates are in pixels.
        """
        (xmin,ymin,xmax,ymax) = xxx_todo_changeme
        (xmin_pix, ymin_pix, xmax_pix, ymax_pix) = xxx_todo_changeme1
        self.canvas.create_line(xmin_pix, ymax_pix, xmax_pix, ymax_pix)
        self.canvas.create_text(xmin_pix   ,ymax_pix+ 8, text=('%4.1f' % xmin))
        self.canvas.create_text(xmin_pix-12,ymax_pix   , text=('%4.1f' % ymin))
        self.canvas.create_line(xmax_pix, ymax_pix, xmax_pix, ymin_pix)
        self.canvas.create_text(xmax_pix   ,ymax_pix+ 8, text=('%4.1f' % xmax))
        self.canvas.create_line(xmax_pix, ymin_pix, xmin_pix, ymin_pix)
        self.canvas.create_text(xmin_pix-12,ymin_pix   , text=('%4.1f' % ymax))
        self.canvas.create_line(xmin_pix, ymin_pix, xmin_pix, ymax_pix)


