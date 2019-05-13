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
#plt.style.use('seaborn')
#plt.tight_layout(pad = 1.25)

# Libraries for stochastic simulation
import numpy as np
from numpy import log, exp
from scipy.stats import skew, kurtosis
from math import ceil, floor, sqrt
from scipy.special import binom
from scipy.stats import norm
import bs_call_base as bsm

def shrink(v, fac):
    '''
    parameters:
        v   = a numpy array
        fac = shrink factor
    purpose:
        shrink the array by fac by removing the head and tail elements from the array
    example:
        v = np.array([1,2,3,4,5,6,7,8,9,10])  # a vector of size 10
        shrink(v, 0.2) should return
        array([2,3,4,5,6,7,8,9]) # a vector of size 10*(1-0.2) = 8 
    '''
    s = int(v.size * (fac/2)) # no. of elements removed from each end 
    return v[s:-s]
    

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
                       'VC': r'$V_C$', 'VP': r'$V_P$', 
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

        self.check_Taylor.clicked.connect(self.on_click_Taylor)
        self.check_Taylor1.clicked.connect(self.on_click_Taylor1)
        self.check_Taylor2.clicked.connect(self.on_click_Taylor2)
        self.check_Theta_T2M.clicked.connect(self.on_click_Theta_T2M)

    def on_click_Theta_T2M(self):
        print("Theta time-to-maturity clicked!")
        if self.greek == 'Θ':
            if self.check_Theta_T2M.isChecked():
                self.oMplCanvas_Theta_T2M.show()
            else:
                self.oMplCanvas_Theta_T2M.hide()
        self.showGraph(graphtype=self.greek)

    def on_click_Taylor(self):
        self.showGraph(graphtype=self.greek)
    def on_click_Taylor1(self):
        self.showGraph(graphtype=self.greek)
    def on_click_Taylor2(self):
        self.showGraph(graphtype=self.greek)

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
        self.slope  = 100
        self.slider_slope.setValue(self.slope)
        self.label_slope.setText(f"S0 = {str(self.slope)}")
        ####
        self.showGraph(graphtype=self.greek)

    def on_change_K(self, value):
        self.K = value
        self.label_K.setText(f"K = {str(self.K)}")
        self.showGraph(self.greek)
    def on_change_slope(self, value):
        self.slope = value
        self.label_slope.setText(f"S0 = {str(self.slope)}")
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

        mpl2 = self.oMplCanvas_Theta_T2M.canvas
        s5 = mpl2.axes[0]

        if self.greek != 'Θ' or not self.check_Theta_T2M.isChecked():
            window.oMplCanvas_Theta_T2M.hide()
        else:
            window.oMplCanvas_Theta_T2M.show()

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
        a, b = bsm.taylor_linear(slope, K, r, q, σ, T, category='call')
        taylor1_call  = a*x + b
        # put parabola
        a, b, d = bsm.taylor_parabola(slope, K, r, q, σ, T, category='put')
        taylor2_put  = a*x*x + b*x + d
        a, b = bsm.taylor_linear(slope, K, r, q, σ, T, category='put')
        taylor1_put  = a*x + b

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

        fac = 0.6 # shrink factor for Taylor polynomials

        s1.clear()
        s1.set_title('Call', fontsize=14, color='brown')
        s1.yaxis.set_tick_params(labelright=False, labelleft=True)
        s1.set_xlabel('$S_0$')
        s1.grid(True)
        s1.plot(x, y_call, 'r', lw=1.5, label='value')
        s1.axvline(x=slope, color="g", alpha=0.5, linewidth=1)
        s1.plot(x, z_call, 'b-.', lw=1.5, label='payoff')
        s1.plot(slope, bsm.bs_call_price(slope, K, r, q, σ, T), 'black', marker="o")
        s1.annotate(f"Δ = {bsm.delta(slope, K, r, q, σ, T, category='call'):.3f}", xy=[slope-20, bsm.bs_call_price(slope, K, r, q, σ, T)+5])
        s1.annotate(f"Γ = {bsm.gamma(slope, K, r, q, σ, T, category='call'):.3f}", xy=[slope-20, bsm.bs_call_price(slope, K, r, q, σ, T)+0])
        if self.check_Taylor.isChecked():
            if self.check_Taylor1.isChecked():
                s1.plot(shrink(x, fac), shrink(taylor1_call, fac), 'orange', lw=1.0, label='Taylor 1st order')
            if self.check_Taylor2.isChecked():
                s1.plot(shrink(x, fac), shrink(taylor2_call, fac), 'grey', lw=1.0, label='Taylor 2nd order')
        legend = s1.legend(loc='upper left', shadow=True, fontsize='medium')


        s2.clear()
        s2.set_title('Put', fontsize=14, color='brown')
        s2.set_xlabel('$S_0$')
        s2.grid(True)
        s2.plot(x, y_put, 'r', lw=1.5, label='value')
        s2.axvline(x=slope, color="g", alpha=0.5, linewidth=1)
        s2.plot(x, z_put, 'b-.', lw=1.5, label='payoff')
        s2.plot(slope, bsm.bs_put_price(slope, K, r, q, σ, T), 'black', marker="o")
        s2.annotate(f"Δ = {bsm.delta(slope, K, r, q, σ, T, category='put'):.3f}", xy=[slope+5, bsm.bs_put_price(slope, K, r, q, σ, T)+5])
        s2.annotate(f"Γ = {bsm.gamma(slope, K, r, q, σ, T, category='put'):.3f}", xy=[slope+5, bsm.bs_put_price(slope, K, r, q, σ, T)+0])
        if self.check_Taylor.isChecked():
            if self.check_Taylor1.isChecked():
                s2.plot(shrink(x, fac), shrink(taylor1_put, fac), 'orange', lw=1.0, label='Taylor 1st order')
            if self.check_Taylor2.isChecked():
                s2.plot(shrink(x, fac), shrink(taylor2_put, fac), 'grey', lw=1.0, label='Taylor 2nd order')
        legend = s2.legend(loc='upper right', shadow=True, fontsize='medium')

        s3.clear()
        s3.get_shared_y_axes().join(s3, s4)
        s3.set_title(self.titles[graphtype+'C'], fontsize=14, color='brown')
        s3.grid(True)
        s3.plot(x, greek_call, 'brown', lw=1.5, label='')
        s3.axvline(x=slope, color="g", alpha=1.0, linewidth=1)
        s3.axhline(y=0, color="black", alpha=1, linewidth=1)

        s4.clear()
        s4.set_title(self.titles[graphtype+'P'], fontsize=14, color='brown')
        s4.grid(True)
        s4.plot(x, greek_put, 'brown', lw=1.5, label='')
        s4.axvline(x=slope, color="g", alpha=1.0, linewidth=1)
        s4.axhline(y=0, color="black", alpha=1, linewidth=1)

        mpl.fig.set_visible(True)
        mpl.draw()

        ## time-to-maturity Theta graph
        s5.clear()
        if self.check_Theta_T2M.isChecked() and self.greek == 'Θ':
            self.oMplCanvas_Theta_T2M.show()
            #print("Cool!")
            S0 = slope
            s5.set_title("Variation of Θ with time to maturity", fontsize=10, color='brown')
            s5.grid(True)
            s5.get_shared_y_axes().remove(s3)
            t = np.linspace(0, 0.99, 100) # time to maturity
            K_vec = [S0-10, S0, S0+10] 
            Θ_vec = [None, None, None]
            Θ_vec[0] = bsm.theta(S0, K_vec[0], r, q, σ, T-t, category='call')
            Θ_vec[1] = bsm.theta(S0, K_vec[1], r, q, σ, T-t, category='call')
            Θ_vec[2] = bsm.theta(S0, K_vec[2], r, q, σ, T-t, category='call')
            s5.plot(t, Θ_vec[0], 'brown', lw=0.7, label=f'K={K_vec[0]}')
            s5.plot(t, Θ_vec[1], 'green', lw=0.7, label=f'K={K_vec[1]}')
            s5.plot(t, Θ_vec[2], 'grey', lw=0.7, label=f'K={K_vec[2]}')
            legend = s5.legend(loc='lower left', shadow=True, fontsize='small')

        mpl2.fig.set_visible(True)
        mpl2.draw()

# This creates the GUI window:
if __name__ == '__main__':

    import sys
    app = QtWidgets.QApplication(sys.argv)
    ### load appropriately sized UI based on screen resolution detected
    screen = app.primaryScreen()
    screen_width=screen.size().width()
    window = Window()
    #window.showGraph('Δ')
    window.oMplCanvas_Theta_T2M.hide()
    window.reset()
    window.show()
    sys.exit(app.exec_())

