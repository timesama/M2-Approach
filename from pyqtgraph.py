from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

pg.mkQApp()

# Axis
a2 = pg.AxisItem("left")
a3 = pg.AxisItem("left")
a4 = pg.AxisItem("left")  
a5 = pg.AxisItem("left")
a6 = pg.AxisItem("left")

# ViewBoxes
v2 = pg.ViewBox()
v3 = pg.ViewBox()
v4 = pg.ViewBox()
v5 = pg.ViewBox()
v6 = pg.ViewBox()

# main view
pw = pg.GraphicsView()
pw.setWindowTitle('pyqtgraph example: multiple y-axis')
pw.show()

# layout
l = pg.GraphicsLayout()
pw.setCentralWidget(l)

# add axis to layout
## watch the col parameter here for the position

l.addItem(a2, row=1, col=5, rowspan=1, colspan=1)
l.addItem(a3, row=1, col=4, rowspan=1, colspan=1)
l.addItem(a4, row=1, col=3, rowspan=1, colspan=1)
l.addItem(a5, row=1, col=2, rowspan=1, colspan=1)
l.addItem(a6, row=1, col=1, rowspan=1, colspan=1)

# Blank axis used for aligning things
ax = pg.AxisItem(orientation='bottom')
ax.setPen('#000000')
pos = (2,2)
l.addItem(ax, *pos)

# plotitem and viewbox
## at least one plotitem is used which holds its own viewbox and left axis

pI = pg.PlotItem()
v1 = pI.vb  # reference to viewbox of the plotitem

l.addItem(pI, row=1, col=8, rowspan=2, colspan=1)  # add plotitem to layout

# split off 1st axis and put to side
pI.axis_left = pI.getAxis('left')
pos = (1,7)
l.addItem(pI.axis_left, *pos)