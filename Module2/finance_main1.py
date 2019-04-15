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

from mpl_toolkits.axes_grid1 import make_axes_locatable

# Libraries for stochastic simulation
import numpy as np
from numpy import log, exp
from scipy.stats import skew, kurtosis
from math import ceil, floor, sqrt
from scipy.special import binom
from scipy.stats import norm

###########################
## Simulation modules
#from Simulation import binomial_randomwalk_multiplicative_simulate as brwm
from Simulation import geometric_brownian_motion_simulate as gbm

fontspec1 = {'family': 'sans serif',
             'name':   'Arial',
             'color':  'darkred',
             'weight': 'normal',
             'size': 8,
            }

def A(k, u, n, p):
    ''' helper function to evaluate analytical statistics
    '''
    return 1 + (u**k - 1)*p

def rw_mean1a(y0, u, n, p, K):
    Δx = log(u)
    x0 = log(y0)
    K  = log(K)
    def summand(m):
        return binom(n, m) * p**m * (1-p)**(n-m) * ( (2*m - n)*Δx + x0 - K )
    # lower index
    index0 = max(0, ceil( (K - (x0 - n*Δx)) / (2*Δx) ) )
    sequence = [summand(i) for i in range(index0, n+1)]
    return sum(sequence)

def rw_mean2a(y0, u, n, p, K):
    Δx = log(u)
    x0 = log(y0)
    K  = log(K)
    def summand(m):
        return binom(n, m) * p**m * (1-p)**(n-m) * ( K - x0 - (2*m - n)*Δx )
    # upper index
    indexN = min(n, ceil( (K - (x0 - n*Δx)) / (2*Δx) ) )
    sequence = [summand(i) for i in range(0, indexN+1)]
    return sum(sequence)

def grw_mean1a(y0, u, n, p, K):
    def summand(m):
        return binom(n, m) * p**m * (1-p)**(n-m) * ( y0*u**(2*m-n) - K )
    # lower index
    index0 = max( 0, ceil( (log(K) - (log(y0) - n*log(u))) / (2*log(u)) ) ) 
    sequence = [summand(i) for i in range(index0, n+1)]
    return sum(sequence)   

def grw_mean2a(y0, u, n, p, K):
    def summand(m):
        return binom(n, m) * p**m * (1-p)**(n-m) * ( K - y0*u**(2*m-n) )
    # upper index
    indexN = min( n, ceil( (log(K) - (log(y0) - n*log(u))) / (2*log(u)) ) ) 
    sequence = [summand(i) for i in range(0, indexN+1)]
    return sum(sequence)    

def theoretical_statistics_BM(x0, α, σ, t, K):
    μ = α + σ*σ/2
    theo_mean = x0 + μ * t
    theo_var  = σ*σ * t
    theo_skew = 0
    theo_kurt = 3
    theo_mean1a = σ*σ*t * norm.pdf((x0+μ*t-K)/(σ*sqrt(t))) + (x0+μ*t-K) * norm.cdf((x0+μ*t-K)/(σ*sqrt(t)))
    theo_mean2a = σ*σ*t * norm.pdf((-x0-μ*t+K)/(σ*sqrt(t))) + (-x0-μ*t+K) * norm.cdf((-x0-μ*t+K)/(σ*sqrt(t)))
    return (theo_mean, theo_var, theo_skew, theo_kurt, theo_mean1a, theo_mean2a)

def theoretical_statistics_GBM(y0, α, σ, t, K):
    μ   = α + σ*σ/2
    ν_t = (μ - σ*σ/2) * t
    β   = (log(K/y0) - ν_t) / (σ*sqrt(t))
    theo_mean = y0 * exp(μ * t)
    theo_var  = y0*y0 * ( exp(σ*σ * t) - 1 ) * exp(2*μ*t)
    theo_skew = ( exp(σ*σ * t) + 2 ) * sqrt( exp(σ*σ * t) - 1 )
    theo_kurt = exp(4*σ*σ*t) + 2*exp(3*σ*σ*t) + 3*exp(2*σ*σ*t) - 3
    theo_mean1a = y0*exp(μ*t)*norm.cdf( -β + σ*sqrt(t) ) - K*norm.cdf( -β )
    theo_mean2a = K*norm.cdf( β ) - y0*exp(μ*t)*norm.cdf( β - σ*sqrt(t) )
    return (theo_mean, theo_var, theo_skew, theo_kurt, theo_mean1a, theo_mean2a)

def theoretical_statistics_GRW(y0, u, n, p, K):
    theo_mean = (y0/u**n) * A(2, u, n, p)**n
    theo_var  = (y0/u**n)**2 * ( A(4, u, n, p)**n - A(2, u, n, p)**(2*n) )
    theo_skew = ( (y0/u**n)**3 * ( A(6, u, n, p)**n - 3 * A(4, u, n, p)**n * A(2, u, n, p)**n + 2*A(2, u, n, p)**(3*n) ) ) / theo_var**(3/2)
    theo_kurt = (y0/u**n)**4 * ( A(8, u, n, p)**n - 4 * A(6, u, n, p)**n * A(2, u, n, p)**n + 6 * A(4, u, n, p)**n * A(2, u, n, p)**(2*n) - 3 * A(2, u, n, p)**(4*n) ) / theo_var**2
    theo_mean1a = grw_mean1a(y0, u, n, p, K)
    theo_mean2a = grw_mean2a(y0, u, n, p, K)
    return (theo_mean, theo_var, theo_skew, theo_kurt, theo_mean1a, theo_mean2a)

def theoretical_statistics_RW(y0, u, n, p, K):
    theo_mean = log(y0) + n * log(u) * (2*p - 1)  # Δx = log(y0)
    theo_var  = 4 * log(u)**2 * n*p*(1 - p)
    theo_skew = 8 * log(u)**3 * n*p*(1 - p)*(1 - 2*p) / (theo_var**(3/2))
    theo_kurt = 16 * log(u)**4 * n*p*(1 - p) * ( 1 + 3*(n - 2)*p - 3*(n - 2)*p**2 ) / (theo_var**2)
    theo_mean1a = rw_mean1a(y0, u, n, p, K)
    theo_mean2a = rw_mean2a(y0, u, n, p, K)
    return (theo_mean, theo_var, theo_skew, theo_kurt, theo_mean1a, theo_mean2a)


class Window(QtWidgets.QDialog):
#class Window(QtWidgets.QScrollArea):
    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        # load the .ui file created in Qt Creator directly
        uic.loadUi('mod2.ui', self)
        #initialize(self)
        self.DEBUG = True
        self.dirty = True # this flag indicates SIMULATION is REQUIRED
        self.deltaTdict = {0:0.1, 1:0.05, 2:0.01, 3:0.005}
        # initialization parameters
        self.randomwalk_params = {'S0': 100, 'μ': 0.05, 'σ': 0.2, 'Δt': 0.01, 'T': 1, 'N':100000, 'P': 30, 'K': 100}
        self.S = None # this will store the current set of simulated values
        self.X = None
        self.Y = None
        self.t1t2_max = None
        self.lastDeltaT = self.randomwalk_params['Δt']
        self.plot_title = None
        self.Model = 'GRW' # eith 'RW' or 'GRW'
        self.pushButton_ShowRandomWalk.clicked.connect(self.on_click_ShowGeometricBrownianMotion)
        #self.figure = self.MplCanvas
        self.oSlider_mu.valueChanged.connect(self.changeValue_mu)
        self.lcdNum_prob.display(self.oSlider_mu.value())
        #self.oSlider_mu.valueChanged['int'].connect(self.lcdNum_prob.display)
        self.label_sigma.setText(str(self.oSlider_sigma.value()))
        #self.delta_x = log(self.oSlider_sigma.value())
        self.T_Max = 200 # maximum number of time step to be simulated
        #self.self.oSlider_sigma.valueChanged['int'].connect(self.label_sigma.text)
        ### random walk parameters
        additive_model = {}
        multiplicative_model = self.randomwalk_params
        additive_model['S0'] = log(multiplicative_model['S0']) 
        additive_model['σ']  = multiplicative_model['σ'] 
        additive_model['Δt'] = multiplicative_model['Δt'] 
        additive_model['K']  = log(multiplicative_model['K']) 
        additive_model['μ']  = multiplicative_model['μ'] 
        additive_model['T']  = multiplicative_model['T'] 
        additive_model['N']  = multiplicative_model['N'] 
        additive_model['P']  = multiplicative_model['P'] 
        self.randomwalk_params_additive = additive_model
        #self.slider_.setText(f"{self.randomwalk_params['']}")
        self.label_S0.setText(f"x0 = {log(self.randomwalk_params['S0']):.2f}\ny0 = {self.randomwalk_params['S0']}")
        self.label_P.setText(f"Paths = {self.randomwalk_params['P']}")
        self.label_T.setText(f"T = {self.randomwalk_params['T']}")
        self.label_K.setText(f"ln(K) ={log(self.randomwalk_params['K'])}\nK = {self.randomwalk_params['K']}")
        self.label_sigma.setText(f"σ = {self.randomwalk_params['σ']:.2f}")
        self.label_mu.setText(f"α = {self.randomwalk_params['μ']:.2f}")
        #self.slider_.setText(f"{self.randomwalk_params['']}")
        self.label_deltaT.setText(f"Δt = {self.randomwalk_params['Δt']}")
        
        self.comboBox_SimSize.currentTextChanged.connect(self.changeValue_comboBox_SimSize)
        self.oSlider_sigma.valueChanged['int'].connect(self.changeValue_oSlider_sigma)

        self.oSlider_TimeStep.setValue(self.randomwalk_params['T'])
        self.oSlider_TimeStep.valueChanged.connect(self.changeValue_oSlider_TimeStep)

        self.oSlider_PathsShown.valueChanged.connect(self.changeValue_oSlider_PathsShown)
        self.oSlider_deltaT.valueChanged.connect(self.changeValue_oSlider_deltaT)

        # clear canvas
        #self.oMplCanvas.canvas.clear()
        # t1 and t2 sliders
        self.reset_t1t2()
        magScale = self.lastDeltaT / self.randomwalk_params['Δt']
        self.slider_t1.valueChanged.connect(self.changeSliderT1)
        self.slider_t2.valueChanged.connect(self.changeSliderT2)
        ###################
        self.slider_S0.valueChanged.connect(self.changeSliderS0)
        
        self.slider_K.valueChanged.connect(self.changeSliderK)

        #self.btn_fixY.clicked.connect(self.on_click_btn_fixY)

        ###################

        ###### delete ax4   
        cv = self.oMplCanvas.canvas 
        #ax3 = cv.axes[1][0]
        ax4 = cv.axes[1][1]
        #cv.fig.delaxes(ax3)
        cv.fig.delaxes(ax4)
        
        # hide descriptive stats result plot
        #self.oMplCanvas2.Show(False)
        mpltext = self.oMplCanvas2.canvas
        mpltext.fig.delaxes(mpltext.axes[0][1])
        mpltext.fig.delaxes(mpltext.axes[1][0])
        mpltext.fig.delaxes(mpltext.axes[1][1])
        
        ##
        self.OutputStats = {'T': [], 't1': [], 't2': []}

        ## 
        self.btn_Copy2Clipboard.clicked.connect(self.on_click_Copy2Clipboard)
        
        ## RW/GRW radio buttons
        self.radioGRW.clicked.connect(self.on_click_radioGRW)
        self.radioRW.clicked.connect(self.on_click_radioRW)


    def reset_t1t2(self):
        magScale = self.lastDeltaT / self.randomwalk_params['Δt']
        t1t2_max = int(1/self.randomwalk_params['Δt'])
        self.t1t2_max = t1t2_max
        self.slider_t1.setMaximum(t1t2_max)
        self.slider_t1.setSingleStep(t1t2_max//50)
        self.slider_t1.setMinimum(0)
        #self.slider_t1.setValue(self.slider_t1.value() * magScale)
        #self.label_t1.setText(str(self.slider_t1.value()))
        #self.label_t1.setText('1')
        ###################
        self.slider_t2.setMaximum(t1t2_max)
        self.slider_t2.setSingleStep(t1t2_max//50)
        self.slider_t2.setMinimum(0)
        #self.slider_t2.setValue(self.slider_t2.value() * magScale)
        #self.label_t2.setText(str(self.slider_t2.value()))
        #self.label_t2.setText('1')


    def on_click_radioGRW(self):
        #if self.S == None:
        #    return
        #else:
        #if not self.radioGRW.isChecked():
        if isinstance(self.S, np.ndarray):
            self.S = self.Y
            if self.DEBUG: print(f"on_click_radioGRW")
            self.showGraph()
        
    def on_click_radioRW(self):
        #if self.S == None:
        #    return
        #else:
        #if not self.radioRW.isChecked():
        if isinstance(self.S, np.ndarray):
            self.S = self.X
            if self.DEBUG: print(f"on_click_radioRW")
            self.showGraph()
        
    def on_click_Copy2Clipboard(self):
        cp = QApplication.clipboard()
        descr_stats_labels = ['E(X)', 'Var(X)', 'Skew(X)', 'Kurt(X)', 'E[(X - K)+]', 'E[(K - X)+]']
        txt = ''
        for i in range(6):
            txt += descr_stats_labels[i] + "\t" + str(self.OutputStats['T'][i]) + "\n"
        cp.setText(txt)
        
        #FigureCanvas.__init__(mpltext.fig)
        #FigureCanvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def changeSliderS0(self, value):
        self.label_S0.setText(f"x0 = {log(value):.2f}\ny0 = {value}")
        self.randomwalk_params['S0'] = value

    def changeSliderK(self, value):
        self.label_K.setText(f"log(K) = {log(value):.2f}\nK = {value}")
        self.randomwalk_params['K'] = value

    def changeSliderT1(self, value):
        if self.dirty: # do nothing if dirty
            return None
        self.label_t1.setText(str(value/self.t1t2_max))
        #self.t1t2_max = int(1/self.randomwalk_params['Δt'])
        if self.DEBUG: print(f"{__name__}")
        self.showGraph()
        
    def changeSliderT2(self, value):
        if self.dirty: # do nothing if dirty
            return None
        self.label_t2.setText(str(value/self.t1t2_max))
        #self.t1t2_max = int(1/self.randomwalk_params['Δt'])
        if self.DEBUG: print(f"changeSliderT2")
        self.showGraph()

    def changeValue_oSlider_TimeStep(self, value):
        self.randomwalk_params['T'] = value
        self.label_T.setText(f"T = {self.randomwalk_params['T']}")
        self.slider_t1.setMaximum(self.oSlider_TimeStep.value())
        self.slider_t2.setMaximum(self.oSlider_TimeStep.value())
        #print(f"'T' clicked: {value}")
        
        if isinstance(self.S, np.ndarray):
            if self.DEBUG: print(f"changeValue_oSlider_TimeStep")
            self.showGraph()

    def changeValue_oSlider_PathsShown(self, value):
        self.randomwalk_params['P'] = value
        self.label_P.setText(f"Paths = {self.randomwalk_params['P']}")
        if isinstance(self.S, np.ndarray):
            if self.DEBUG: print(f"changeValue_oSlider_PathsShown")
            self.showGraph()
        
    def changeValue_oSlider_deltaT(self, value):
        self.dirty = True
        self.lastDeltaT = self.randomwalk_params['Δt']
        #self.randomwalk_params['Δt'] = 10**(-value) # 1 => 0.1; 2 => 0.01, etc.
        #self.label_deltaT.setText(f"Δt = {10**(-value)}")
        self.randomwalk_params['Δt'] = self.deltaTdict[value] # 0 => 0.1; 1 => 0.05, etc.
        self.label_deltaT.setText(f"Δt = {self.deltaTdict[value]}")

    def changeValue_oSlider_sigma(self, value):
        self.dirty = True
        self.randomwalk_params['σ'] = value/10
        self.label_sigma.setText(f"σ = {value/10:.2f}")
        
    def changeValue_comboBox_SimSize(self, value):
        self.dirty = True
        self.randomwalk_params['N'] = int(value)

    def changeValue_mu(self, value):
        self.dirty = True
        self.randomwalk_params['μ'] = value/100
        self.label_mu.setText(f"α = {value/100:.2f}")


    def on_click_ShowGeometricBrownianMotion(self):
        # print("I'm clicked!")
        S0 = self.randomwalk_params['S0']
        μ  = self.randomwalk_params['μ']
        σ  = self.randomwalk_params['σ']
        Δt = self.randomwalk_params['Δt']
        T  = self.randomwalk_params['T']
        N  = self.randomwalk_params['N']
        P  = self.randomwalk_params['P']
        s  = 12345
        params = {'S0': S0, 'μ': μ, 'σ': σ, 'Δt': Δt, 'T': T, 'N': N, 'P': P, 'seed': s}
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        self.Y = gbm(params).round(5)
        self.X = log(self.Y).round(5)
        self.dirty = False
        if self.DEBUG: print(f"on_click_ShowGeometricBrownianMotion")
        self.showGraph()
        QApplication.restoreOverrideCursor()

    #def on_click_btn_fixY(self):
    #    mpl = self.oMplCanvas.canvas
    #    mainplot = mpl.axes[0][0]
    #    mainplot.set_ylim(top=400)
    #    mpl.draw()
    #    #mpl.flush_events()
    #    return None
        
    def showGraph(self):
        #self.oMplCanvas.close()
        mpl = self.oMplCanvas.canvas
        mainplot = mpl.axes[0][0]
        mainplot.clear()
        #mainplot.autoscale(enable=True)
        #mainplot.autoscale(enable=False)
        #mainplot.set_adjustable(adjustable='box')
        #mainplot.set_xlim(auto=True)
        #mainplot.set_ylim(auto=True)
        P = self.randomwalk_params['P']
        T = self.randomwalk_params['T']
        Δt = self.randomwalk_params['Δt']
        xvals = np.linspace(0.0, 1.0, int(1/Δt)+1)
        #TimeStepShown = min(T, ?) 
        if self.radioGRW.isChecked():
            self.S = self.Y
            self.plot_title = "Geometric Brownian Motion"
        elif self.radioRW.isChecked():
            self.S = self.X
            self.plot_title = "Brownian Motion"
        #mainplot.plot(self.S[:T+1, :P])
        mainplot.plot(xvals, self.S[:, :P], linewidth=0.5)
        mainplot.grid(True)
        mainplot.set_xlabel('time step', fontsize=8, color='brown')
        mainplot.set_ylabel('price', fontsize=8, color='brown')
        mainplot.tick_params(axis = 'both', which = 'major', labelsize = 6)
        #mainplot.draw()
        #divider = make_axes_locatable(mainplot)
        #####
        #axHist = divider.append_axes("right", 1.25, pad=0.1, sharey=mainplot)
        axHist = mpl.axes[0][1]
        axHist.clear()
        axHist.set_xlim(auto=True)
        axHist.set_ylim(auto=True)
        # Turn off tick labels
        #axHist.set_yticklabels([])
        #axHist.set_xticklabels([])        
        #axHist.hist(S[-1, :N], bins=51, density=True, orientation='horizontal', rwidth=2)
        N = self.randomwalk_params['N']
        axHist.hist(self.S[-1, :N], density=True, bins=101, orientation='horizontal', alpha=0.80, color='brown')
        axHist.yaxis.set_ticks_position("right")
        axHist.xaxis.set_major_formatter(FuncFormatter('{0:.1%}'.format))
        axHist.grid(True)
        axHist.tick_params(axis = 'both', which = 'major', labelsize = 6)
        #μ = self.oSlider_mu.value()/100
        SimSize = int(self.comboBox_SimSize.currentText())
        #mainplot.set_title(self.plot_title + f" $p$ = {p:.2f}", fontsize=12, color='brown')
        mainplot.set_title(self.plot_title, fontsize=12, color='brown')
        #mainplot.set_title(r'$z = \sqrt{x^2+y^2}$', fontsize=14, color='r')
        axHist.set_xlabel('probability', fontsize=8, color='brown')
        #####
        # Now take care of sliders t1 and t2
        self.reset_t1t2()
        #self.slider_t1.setEnabled(False)
        #self.slider_t1.setValue(0)
        #self.slider_t1.setEnabled(True)
        t1 = self.slider_t1.value()
        t2 = self.slider_t2.value()
        mainplot.axvline(x=t1/self.t1t2_max, color="red")
        mainplot.axvline(x=t2/self.t1t2_max, color="blue")
        ax3 = mpl.axes[1][0]
        ax3.clear()
        ax3.set_xlabel('price', fontsize=8, color='brown')
        ax3.tick_params(axis = 'both', which = 'major', labelsize = 6)
        #ax3.yaxis.set_major_formatter(FuncFormatter('{0:.1}'.format))
        MaxStep = self.S.shape[0]-1
        ax3.hist(self.S[min(t1+0, MaxStep), :], density=True, bins=101, orientation='vertical', rwidth=2, color='red', alpha=0.5)
        ax3.hist(self.S[min(t2+0, MaxStep), :], density=True, bins=101, orientation='vertical', rwidth=2, color='blue', alpha=0.5)
        #####
        ## populate tableWidget with descriptive statistics
        table = self.tableWidget
        table.setColumnCount(1)     #Set three columns
        table.setRowCount(6)
        table.setHorizontalHeaderLabels(["Value"])
        descr_stats_labels = ['E(X)', 'Var(X)', 'Skew(X)', 'Kurt(X)', 'E[(X - K)+]', 'E[(K - X)+]']
        table.setVerticalHeaderLabels(descr_stats_labels)
        #for i in range(0, table.rowCount()):
        #    table.setItem(i, 0, QTableWidgetItem(descr_stats_labels[i]))
        samp_mean = self.S[-1, :].mean()
        table.setItem(0, 0, QTableWidgetItem(f"{samp_mean:.2f}"))
        samp_var = self.S[-1, :].var()
        table.setItem(1, 0, QTableWidgetItem(f"{samp_var:.2f}"))
        samp_skew = skew(self.S[-1, :])
        table.setItem(2, 0, QTableWidgetItem(f"{samp_skew:.2f}"))
        samp_kurt = kurtosis(self.S[-1, :], fisher=False)   # 3.0 is subtracted if fisher=True
        table.setItem(3, 0, QTableWidgetItem(f"{samp_kurt:.2f}"))
        #table.setVerticalHeaderLabels(None)
        #table.set_visible(False)
        ## Compute E[(X_n - K)^+] 
        T = self.randomwalk_params['T']
        K = self.randomwalk_params['K']
        if self.radioRW.isChecked():
            K = log(K)
        N = self.randomwalk_params['N']
        X = self.S[T, :] - K
        mean1 =  X[X > 0].sum()/N
        X = self.S[T, :]
        mean1a = (np.maximum(X, K) - K).mean()   ## (a - b)+ = max(a, b) - b
        ## Compute E[(K - X_n)^+]
        X = K - self.S[T, :]
        mean2 =  X[X > 0].sum()/N
        X = self.S[T, :]
        mean2a = (np.maximum(K, X) - X).mean()
        del(X)
        table.setItem(4, 0, QTableWidgetItem(f"{mean1:.2f}"))
        table.setItem(5, 0, QTableWidgetItem(f"{mean2:.2f}"))
        table.resizeColumnsToContents()
        #####
        self.OutputStats['T'] = [samp_mean, samp_var, samp_skew, samp_kurt, mean1a, mean2a]
        #####
        ## deal with descriptive statistics of main plot
        mpltext = self.oMplCanvas2.canvas
        FigureCanvas.updateGeometry(mpltext)
        simulated_statistics = ['STATISTICS', '  ', '$Simulated$:',
                r'$\mathrm{E}(X_T) = ' + f'{samp_mean:.2f}$',
                r'$\mathrm{Var}(X_T) = ' + f'{samp_var:.2f}$',
                r'$\mathrm{Skew}(X_T) = ' + f'{samp_skew:.2f}$',
                r'$\mathrm{Kurt}(X_T) = ' + f'{samp_kurt:.2f}$',
                r'$\mathrm{E}[(X_T - K)^+] = ' + f'{mean1a:.2f}$',
                r'$\mathrm{E}[(K - X_T)^+] = ' + f'{mean2a:.2f}$',
                ]
        #simulated_statistics.extend([r'line 8', r'line 9', r'line 10', r'line 11', r'line 12', r'line 13'])
        ### Get theoretical values
        y0 = self.randomwalk_params['S0']
        n  = self.randomwalk_params['T']
        K  = self.randomwalk_params['K']
        μ  = self.randomwalk_params['μ']
        σ  = self.randomwalk_params['σ']
        Δt = self.randomwalk_params['Δt']
        K  = self.randomwalk_params['K']
        print(self.randomwalk_params)
        if self.radioGRW.isChecked():
            #(theo_mean, theo_var, theo_skew, theo_kurt, theo_mean1a, theo_mean2a) = theoretical_statistics_GRW(y0, u, n, p, K)
            (theo_mean, theo_var, theo_skew, theo_kurt, theo_mean1a, theo_mean2a) = theoretical_statistics_GBM(y0, μ, σ, T, K)
        elif self.radioRW.isChecked():
            #(theo_mean, theo_var, theo_skew, theo_kurt, theo_mean1a, theo_mean2a) = theoretical_statistics_RW(y0, u, n, p, K)
            (theo_mean, theo_var, theo_skew, theo_kurt, theo_mean1a, theo_mean2a) = theoretical_statistics_BM(y0, μ, σ, T, K)
        else:
            theo_mean, theo_var, theo_skew, theo_kurt, theo_mean1a, theo_mean2a = (0,0,0,0,0,0)
        analytical_statistics = ['  ', '$Analytical$:',
                r'$\mathrm{E}(X_T) = ' + f'{theo_mean:.2f}$',
                r'$\mathrm{Var}(X_T) = ' + f'{theo_var:.2f}$',
                r'$\mathrm{Skew}(X_T) = ' + f'{theo_skew:.2f}$',
                r'$\mathrm{Kurt}(X_T) = ' + f'{theo_kurt:.2f}$',
                r'$\mathrm{E}[(X_T - K)^+] = ' + f'{theo_mean1a:.2f}$',
                r'$\mathrm{E}[(K - X_T)^+] = ' + f'{theo_mean2a:.2f}$',
                ]
        #statistics = list(simulated_statistics.append([r'-----']))
        statistics = list(simulated_statistics)
        statistics.extend(analytical_statistics)
        #print('stats: ', statistics)
        #mpltext.fig = plt.figure()
        plt.style.use('seaborn-white')
        ax = mpltext.axes[0][0]
        ax.clear()
        #ax.plot([0, 0], 'r')
        #ax.set_title('Test')
        ax.axis([0, 3, -len(statistics), 0])
        ax.set_yticks(-np.arange(len(statistics)))
        for i, s in enumerate(statistics):
            ax.text(0.1, -(i+1)*1.5, s, fontsize=14)
        plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        mpltext.fig.set_visible(True)
        mpltext.draw()
        #mpl2.close()
        #mpl2.fig.show(False)
        #####
        mpl.fig.set_visible(True)
        if self.checkBox_yLim.isChecked():
            new_yLim = int(self.txt_yLim.toPlainText())
            mainplot.set_ylim(top=new_yLim)
        mpl.draw()

        # Matplotlib canvas class to create figure

## The code below show bounds of x-coordinate values and y-coordinate values
# window.oMplCanvas.canvas.ax.dataLim.get_points()


# This creates the GUI window:
if __name__ == '__main__':

    import sys
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    #window.showGraph()
    window.show()
    #window.initialize()
    #sys.exit(app.exec_())
