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
    r = rho**3
    g = the
    b = 1. - rho**3
    return "#%02x%02x%02x" % (r*255, g*255, b*255)

class Plot:

    def __init__(self, root, grid, title='', boxsize=[], width=800, height=800):
        
        if boxsize:
            self.xmin, self.ymin, self.xmax, self.ymax = boxsize
        else:
            self.xmin, self.ymin, self.xmax, self.ymax = grid.boxsize()
        self.width = width
        self.height = height
        self.title = title
        self.grid = grid
    
        self.frame = Frame(root)
        self.frame.pack()
        self.canvas = Canvas(bg="white", width=width, height=height)
        self.canvas.pack()

    def draw(self, rho, the, xyPath=[]):
        
        border_x, border_y = 0.2, 0.2
        scale = min((1.-border_x)*self.width/(self.xmax-self.xmin), \
                    (1.-border_y)*self.height/(self.ymax-self.ymin))
        a = max(border_x * self.width/2., (self.width-scale*(self.xmax-self.xmin))/2.)
        b = max(border_y * self.height/2.,(self.height-scale*(self.ymax-self.ymin))/2.)
    
        box_pix = (a, self.height-scale*(self.ymax-self.ymin) - b,
                   scale*(self.xmax-self.xmin) + a, self.height - b)
        # draw the box
        self.addBox((self.xmin, self.ymin, self.xmax, self.ymax), box_pix)
    
        # color plot
        cells = cell.cell(self.grid)
        for index in range(len(cells.data)):
                ia, ib, ic = cells.data[index]
                xa, ya = scale*(self.grid.x(ia)-self.xmin) + a, \
                         self.height -scale*(self.grid.y(ia)-self.ymin) - b
                xb, yb = scale*(self.grid.x(ib)-self.xmin) + a, \
                         self.height -scale*(self.grid.y(ib)-self.ymin) - b
                xc, yc = scale*(self.grid.x(ic)-self.xmin) + a, \
                         self.height -scale*(self.grid.y(ic)-self.ymin) - b

                rhoAbc = (rho[ia] + rho[ib] + rho[ic])/3.0
                theAbc = (the[ia] + the[ib] + the[ic])/3.0
                color = colormap(rhoAbc, theAbc)
                self.canvas.create_polygon(xa, ya, xb, yb, xc, yc, fill=color)
        
        # add path if present
        index = 0
        for (x, y) in xyPath:
            xa, ya = scale*(x - self.xmin) + a, \
                self.height -scale*(y - self.ymin) - b
            self.canvas.create_oval(xa-5, ya-5,
                                    xa+5, ya+5)
            #self.canvas.create_text(xa+2, ya+1, text=str(index))
            index += 1
    
        # add title
        x, y = self.width/3., self.height/15.0
        self.canvas.create_text(x, y, text=str(self.title),
                                font =("Helvetica", 24),
                                fill='black')


    def addBox(self, xxx_todo_changeme, xxx_todo_changeme1):
        """
        Add Box. Min/Max coordinates are in pixels.
        """
        (self.xmin,self.ymin,self.xmax,self.ymax) = xxx_todo_changeme
        (self.xmin_pix, self.ymin_pix, self.xmax_pix, self.ymax_pix) = xxx_todo_changeme1
        self.canvas.create_line(self.xmin_pix, self.ymax_pix, self.xmax_pix, self.ymax_pix)
        self.canvas.create_text(self.xmin_pix   ,self.ymax_pix+ 8, text=('%4.1f' % self.xmin))
        self.canvas.create_text(self.xmin_pix-12,self.ymax_pix   , text=('%4.1f' % self.ymin))
        self.canvas.create_line(self.xmax_pix, self.ymax_pix, self.xmax_pix, self.ymin_pix)
        self.canvas.create_text(self.xmax_pix   ,self.ymax_pix+ 8, text=('%4.1f' % self.xmax))
        self.canvas.create_line(self.xmax_pix, self.ymin_pix, self.xmin_pix, self.ymin_pix)
        self.canvas.create_text(self.xmin_pix-12,self.ymin_pix   , text=('%4.1f' % self.ymax))
        self.canvas.create_line(self.xmin_pix, self.ymin_pix, self.xmin_pix, self.ymax_pix)


