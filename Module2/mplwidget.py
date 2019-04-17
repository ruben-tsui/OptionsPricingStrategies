# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 13:59:22 2019

@author:
Understanding optics with Python - [Multidisciplinary & applied optics] Ammar, Ahmed_ Lakshminarayanan, Vasudevan_ Varadharajan, L. Srinivasa - CRC 2017.pdf
"""
import os
os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "0"

'''MATPLOTLIB WIDGET'''
# Python Qt5 bindings for GUI objects
from PyQt5.QtWidgets import QSizePolicy, QWidget, QVBoxLayout
# import the Qt5Agg FigureCanvas object, that binds Figure to
# Qt5Agg backend. It also inherits from QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# Matplotlib Toolbar
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
# Matplotlib Figure object
from matplotlib.figure import Figure
from matplotlib import rcParams, gridspec
rcParams['font.size'] = 9
#rcParams['toolbar'] = 'None'
rcParams["figure.edgecolor"] = 'darkblue'

import matplotlib.pyplot as plt
#plt.tight_layout(pad = 1.25)


class MplCanvas(FigureCanvas):
    """
    Class to represent the FigureCanvas widget
    """

    def __init__(self):
        # setup Matplotlib Figure and Axis
        #self.fig = Figure()
        self.fig = plt.figure(figsize=(6, 6))
        plt.subplots_adjust(left=0.05, bottom=0.05, right=1, top=0.95, wspace=0, hspace=0)
        #plt.tight_layout(pad = 1.1)
        #self.fig.tight_layout(pad = 1.5)

        #self.ax = self.fig.add_subplot(111)
        #########
        '''
        fig1, axes = plt.subplots(nrows=2, ncols=2, sharey=True, gridspec_kw={'width_ratios':[2, 1], 'height_ratios': [2,1]})
        fig1.set_figheight(5)
        fig1.set_figwidth(8)
        '''
        #gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[2,1]) 
        gs = plt.GridSpec(3, 3, wspace=0.1, hspace=0.3)
        #plt.axis('tight')
        ax1 = self.fig.add_subplot(gs[:2, :2])
        ax2 = self.fig.add_subplot(gs[:2, 2:], sharey=ax1)
        ax3 = self.fig.add_subplot(gs[2, :2])
        ax4 = self.fig.add_subplot(gs[2, 2:])
        #self.fig.subplots(2,2, gridspec_kw={'width_ratios':[2, 1], 'height_ratios': [2,1]})
        #ax1 = self.fig.add_subplot(2, 2, 1, gs[0])
        #ax2 = self.fig.add_subplot(2, 2, 2, sharey = ax1, gs[1])
        #ax3 = self.fig.add_subplot(2, 2, 3, gs[2])
        #ax4 = self.fig.add_subplot(2, 2, 4, gs[3])
        #'''
        #ax1 = plt.subplot(gs[0,0])
        #ax2 = plt.subplot(gs[0,1])
        #ax3 = plt.subplot(gs[1,0])
        #ax4 = plt.subplot(gs[1,1])
        #'''
        #plt.tight_layout()
        #ax1 = axes[0][0]
        #ax2 = axes[0][1]
        #self.fig = fig1
        axes = [[ax1, ax2], [ax3, ax4]]
        self.axes = axes
        #########
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        # initialization of the canvas
        FigureCanvas.__init__(self, self.fig)
        # we define the widget as expandable
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        # notify the system of updated policy
        FigureCanvas.updateGeometry(self)

class MPL_WIDGET(QWidget):
    """Widget defined in Qt Designer"""
    def  __init__(self, parent=None):
        # initialization of Qt MainWindow widget
        QWidget.__init__(self, parent)
        # set the canvas to the Matplotlib widget
        self.canvas = MplCanvas()
        # create a navigation toolbar for our plot canvas
        #self.navi_toolbar = NavigationToolbar(self.canvas, self)
        # create a vertical box layout
        self.vbl = QVBoxLayout()
        # add mpl widget to vertical box
        self.vbl.addWidget(self.canvas)
        # add the navigation toolbar to vertical box
        #self.vbl.addWidget(self.navi_toolbar)
        # set the layout to vertical box
        self.setLayout(self.vbl)
        #
        self.canvas.fig.set_visible(False)


