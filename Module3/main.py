# -*- coding: utf-8 -*-

import time, datetime, os
#os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"

import sys
#import design
#from PyQt5.QtCore import Qt
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QHBoxLayout, QSlider, QLineEdit, QLabel, QTableWidgetItem, QSizePolicy
from matplotlib.widgets import Slider, Button, RadioButtons

# Libraries for drawing figures
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.ticker import FuncFormatter
# Ensure using PyQt5 backend
import matplotlib as mpl
mpl.use('QT5Agg')
import matplotlib.pyplot as plt
#plt.style.use('ggplot')
plt.style.use('seaborn')
#plt.tight_layout(pad = 1.25)

# Libraries for stochastic simulation
import numpy as np
from numpy import log, exp
from scipy.stats import skew, kurtosis
from math import ceil, floor, sqrt
from scipy.special import binom
from scipy.stats import norm
import bs_call_base as bsm

class Window(QtWidgets.QDialog):
    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        # load the .ui file created in Qt Creator directly
        print(f"screen width = {screen_width}")
        if screen_width >= 1920:
            Interface = uic.loadUi('mod3.ui', self)
        else:
            Interface = uic.loadUi('mod3.ui', self)

        # other useful into
        self.titles = {'ΔC': r'$\Delta_C$', 'ΔP': r'$\Delta_P$', 
                       'ΓC': r'$\Gamma_C$', 'ΓP': r'$\Gamma_P$',
                       'VC': r'$Vega_C$', 'VP': r'$Vega_P$', 
                       'ΘC': r'$\Theta_C$', 'ΘP': r'$\Theta_P$', 
                       'ρC': r'$\rho_C$', 'ρP': r'$\rho_P$'}
        self.greek = 'Δ' # current greek
        # default values for all parameters
        self.K = 110
        self.r = 0.02
        self.q = 0.01
        self.σ = 0.20
        self.T = 1
        self.slope = self.slider_slope
 
        self.radioDelta.clicked.connect(self.on_click_radioDelta)
        self.radioGamma.clicked.connect(self.on_click_radioGamma)
        self.radioVega.clicked.connect(self.on_click_radioVega)
        self.radioTheta.clicked.connect(self.on_click_radioTheta)
        self.radioRho.clicked.connect(self.on_click_radioRho)
        
        self.slider_K.valueChanged.connect(self.on_change_K)
        self.slider_r.valueChanged.connect(self.on_change_r)
        self.slider_q.valueChanged.connect(self.on_change_q)
        self.slider_sigma.valueChanged.connect(self.on_change_sigma)
        self.slider_T.valueChanged.connect(self.on_change_T)
        self.slider_slope.valueChanged.connect(self.on_change_slope)
        self.button_reset.clicked.connect(self.reset)

    def reset(self):
        self.K      = 100
        self.slider_K.setValue(self.K)
        self.label_K.setText(f"K = {str(self.K)}")
        ####
        self.r      = 0.02
        self.slider_r.setValue(self.r * 100)
        self.label_r.setText(f"r = {str(self.r)}")
        ####
        self.q      = 0.01
        self.slider_q.setValue(self.q * 100)
        self.label_q.setText(f"q = {str(self.q)}")
        ####
        self.σ      = 0.2
        self.slider_sigma.setValue(self.σ * 100)
        self.label_sigma.setText(f"σ = {str(self.σ)}")
        ####
        self.T      = 1
        self.slider_T.setValue(self.T * 100)
        self.label_T.setText(f"T = {str(self.T)}")
        ####
        self.slope  = 85
        self.slider_slope.setValue(self.slope)
        self.label_slope.setText(f"Slope = {str(self.slope)}")
        ####
        self.showGraph(graphtype=self.greek)

    def on_change_K(self, value):
        self.K = value
        self.label_K.setText(f"K = {str(self.K)}")
        self.showGraph(self.greek)
    def on_change_slope(self, value):
        self.slope = value
        self.label_slope.setText(f"Slope = {str(self.slope)}")
        self.showGraph(self.greek)
    def on_change_sigma(self, value):
        self.σ = value/100
        self.label_sigma.setText(f"σ = {str(self.σ)}")
        self.showGraph(self.greek)
    def on_change_r(self, value):
        self.r = value/100
        self.label_r.setText(f"r = {str(self.r)}")
        self.showGraph(self.greek)
    def on_change_q(self, value):
        self.q = value/100
        self.label_q.setText(f"q = {str(self.q)}")
        self.showGraph(self.greek)
    def on_change_T(self, value):
        self.T = value/10
        self.label_T.setText(f"T = {str(self.T)}")
        self.showGraph(self.greek)

    def on_click_radioDelta(self):
        self.greek = 'Δ'
        self.showGraph(graphtype=self.greek)
    def on_click_radioGamma(self):
        self.greek = 'Γ'
        self.showGraph(graphtype=self.greek)
    def on_click_radioTheta(self):
        self.greek = 'Θ'
        self.showGraph(graphtype=self.greek)
    def on_click_radioRho(self):
        self.greek = 'ρ'
        self.showGraph(graphtype=self.greek)
    def on_click_radioVega(self):
        self.greek = 'V'
        self.showGraph(graphtype=self.greek)



    def showGraph(self, graphtype):

        mpl = self.oMplCanvas.canvas
        [[s1, s2], [s3, s4]] = mpl.axes

        #S0 = 100; K = 100; r = 0.025; q = 0.005; σ = 0.2; T = 1.0
        K      = self.K
        r      = self.r
        q      = self.q
        σ      = self.σ
        T      = self.T
        slope  = self.slope 
        x      = np.linspace(40, 160, 51)
        y_call = bsm.bs_call_price(x, K, r, q, σ, T)
        y_put  = bsm.bs_put_price(x, K, r, q, σ, T)
        z_call = bsm.payoff_call(x, K)
        z_put  = bsm.payoff_put(x, K)

        # call parabola
        print(f"Parabola K={K}, r={r}, q={q}, σ={σ}, T={T}")
        a, b, d = bsm.taylor_parabola(slope, K, r, q, σ, T, category='call')
        print(f'a={a:.2f}, b={b:.2f}, d={d:.2f}')
        taylor2_call = a*x*x + b*x + d
        # put parabola
        if self.check_Taylor.isChecked():
            a, b, d = bsm.taylor_parabola(slope, K, r, q, σ, T, category='put')
            taylor2_put  = a*x*x + b*x + d

        greek_call = None; greek_put = None
        if graphtype == 'Δ':
            greek_call = delta_call =  np.exp(-q*T) * bsm.N(  bsm.d1(x, K, r, q, σ, T) )
            greek_put  = delta_put  = -np.exp(-q*T) * bsm.N( -bsm.d1(x, K, r, q, σ, T) )
        elif graphtype == 'Γ':
            greek_call = greek_put = gamma_call =  bsm.gamma(x, K, r, q, σ, T, category='call')
        elif graphtype == 'V':
            greek_call = greek_put = vega_put   =  bsm.vega(x, K, r, q, σ, T, category='put')
        elif graphtype == 'Θ':
            greek_call = theta_call =  bsm.theta(x, K, r, q, σ, T, category='call')
            greek_put  = theta_put  =  bsm.theta(x, K, r, q, σ, T, category='put')
        elif graphtype == 'ρ':
            greek_call = rho_call =  bsm.rho(x, K, r, q, σ, T, category='call')
            greek_put  = rho_put  =  bsm.rho(x, K, r, q, σ, T, category='put')

        s1.clear()
        s1.set_title('Call', fontsize=14, color='brown')
        s1.yaxis.set_tick_params(labelright='off', labelleft='on')
        s1.grid(True)
        s1.plot(x, y_call, 'r', lw=2.5, label='present value')
        s1.axvline(x=slope, color="g", alpha=0.75, linewidth=2)
        s1.plot(x, z_call, 'b-.', lw=1.5, label='payoff')
        if self.check_Taylor.isChecked():
            s1.plot(x, taylor2_call, 'purple', lw=1.5, label='taylor 2nd order')

        s2.clear()
        s2.set_title('Put', fontsize=14, color='brown')
        s2.grid(True)
        s2.plot(x, y_put, 'r', lw=2.5, label='present value')
        s2.axvline(x=slope, color="g", alpha=0.75, linewidth=2)
        s2.plot(x, z_put, 'b-.', lw=1.5, label='payoff')
        if self.check_Taylor.isChecked():
            s2.plot(x, taylor2_put, 'purple', lw=1.5, label='taylor 2nd order')

        s3.clear()
        s3.set_title(self.titles[graphtype+'C'], fontsize=14, color='blue')
        s3.grid(True)
        s3.plot(x, greek_call, 'r', lw=2.5, label='delta call')
        s3.axvline(x=slope, color="g", alpha=0.75, linewidth=2)

        s4.clear()
        s4.set_title(self.titles[graphtype+'P'], fontsize=14, color='blue')
        s4.grid(True)
        s4.plot(x, greek_put, 'r', lw=2.5, label='delta put')
        s4.axvline(x=slope, color="g", alpha=0.75, linewidth=2)

        mpl.fig.set_visible(True)
        mpl.draw()

 




# This creates the GUI window:
if __name__ == '__main__':

    import sys
    app = QtWidgets.QApplication(sys.argv)
    ### load appropriately sized UI based on screen resolution detected
    screen = app.primaryScreen()
    screen_width=screen.size().width()
    window = Window()
    #window.showGraph('Δ')
    window.reset()
    window.show()
    sys.exit(app.exec_())

