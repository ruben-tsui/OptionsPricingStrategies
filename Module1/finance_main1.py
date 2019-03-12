# -*- coding: utf-8 -*-

import time, datetime, os
os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"

import sys
#import design
#from PyQt5.QtCore import Qt
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QHBoxLayout, QSlider, QLineEdit, QLabel, QTableWidgetItem, QSizePolicy
from matplotlib.widgets import Slider, Button, RadioButtons

# Imports
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.ticker import FuncFormatter

# Ensure using PyQt5 backend
import matplotlib as mpl
mpl.use('QT5Agg')

import matplotlib.pyplot as plt
plt.style.use('ggplot')
#plt.style.use('seaborn')
#plt.tight_layout(pad = 1.25)

from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
from numpy import log, exp
from scipy.stats import skew, kurtosis

###########################
## Simulation modules
from Simulation import binomial_randomwalk_multiplicative_simulate as brwm

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

class Window(QtWidgets.QDialog):
    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        # load the .ui file created in Qt Creator directly
        uic.loadUi('finance1b.ui', self)
        #initialize(self)
        self.S = None # this will store the current set of simulated values
        self.X = None
        self.Y = None
        self.plot_title = None
        self.Model = 'GRW' # eith 'RW' or 'GRW'
        self.pushButton_ShowRandomWalk.clicked.connect(self.on_click_ShowRandomWalk)
        #self.figure = self.MplCanvas
        self.oSlider_prob.valueChanged.connect(self.changeValue_prob)
        self.lcdNum_prob.display(self.oSlider_prob.value())
        self.oSlider_prob.valueChanged['int'].connect(self.lcdNum_prob.display)
        self.label_u.setText(str(self.oSlider_u.value()))
        self.T_Max = 200 # maximum number of time step to be simulated
        #self.self.oSlider_u.valueChanged['int'].connect(self.label_u.text)
        ### random walk parameters
        # initialization parameters
        self.randomwalk_params = {'S0': 100, 'p': 50, 'u': 105, 'd': 10000/105, 'T': 50, 'N':100000, 'P': 25, 'K': 100}
        additive_model = {}
        multiplicative_model = self.randomwalk_params
        additive_model['S0'] = log(multiplicative_model['S0']) 
        additive_model['u'] = log(multiplicative_model['u']) 
        additive_model['d'] = log(multiplicative_model['d']) 
        additive_model['K'] = log(multiplicative_model['K']) 
        additive_model['p'] = multiplicative_model['p'] 
        additive_model['T'] = multiplicative_model['T'] 
        additive_model['N'] = multiplicative_model['N'] 
        additive_model['P'] = multiplicative_model['P'] 
        self.randomwalk_params_additive = additive_model
        #self.slider_.setText(f"{self.randomwalk_params['']}")
        self.label_S0.setText(f"X0 = {log(self.randomwalk_params['S0'])}\nY0 = {self.randomwalk_params['S0']}")
        self.label_P.setText(f"Paths = {self.randomwalk_params['P']}")
        self.label_T.setText(f"n = {self.randomwalk_params['T']}")
        self.label_K.setText(f"ln(K) ={log(self.randomwalk_params['K'])}\nK = {self.randomwalk_params['K']}")
        self.label_u.setText(f"Δx=\nu = {self.randomwalk_params['u']/100:.2f}")
        self.label_p.setText(f"p = {self.randomwalk_params['p']/100:.2f}")
        #self.slider_.setText(f"{self.randomwalk_params['']}")
        
        self.comboBox_SimSize.currentTextChanged.connect(self.changeValue_comboBox_SimSize)
        self.oSlider_u.valueChanged['int'].connect(self.changeValue_oSlider_u)

        self.oSlider_TimeStep.setValue(self.randomwalk_params['T'])
        self.oSlider_TimeStep.valueChanged.connect(self.changeValue_oSlider_TimeStep)

        self.oSlider_PathsShown.valueChanged.connect(self.changeValue_oSlider_PathsShown)
        
        # clear canvas
        #self.oMplCanvas.canvas.clear()
        # t1 and t2 sliders        
        self.slider_t1.setSingleStep(1)
        self.slider_t1.setMaximum(self.oSlider_TimeStep.value())
        self.slider_t1.setMinimum(1)
        self.slider_t1.setValue(5)
        self.slider_t1.valueChanged.connect(self.changeSliderT1)
        self.label_t1.setText(str(self.slider_t1.value()))
        ###################
        self.slider_t2.setSingleStep(1)
        self.slider_t2.setMaximum(self.oSlider_TimeStep.value())
        self.slider_t2.setMinimum(1)
        self.slider_t2.setValue(15)
        self.slider_t2.valueChanged.connect(self.changeSliderT2)
        self.label_t2.setText(str(self.slider_t2.value()))
        ###################
        self.slider_S0.valueChanged.connect(self.changeSliderS0)
        
        self.slider_K.valueChanged.connect(self.changeSliderK)
        
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

    def on_click_radioGRW(self):
        #if self.S == None:
        #    return
        #else:
        if isinstance(self.S, np.ndarray):
            self.S = self.Y
            self.showGraph()
        
    def on_click_radioRW(self):
        #if self.S == None:
        #    return
        #else:
        if isinstance(self.S, np.ndarray):
            self.S = self.X
            self.showGraph()
        
    def on_click_Copy2Clipboard(self):
        cp = QApplication.clipboard()
        #cp.setText('小謙小若')
        descr_stats_labels = ['E(X)', 'Var(X)', 'Skew(X)', 'Kurt(X)', 'E[(X - K)+]', 'E[(K - X)+]']
        txt = ''
        for i in range(6):
            txt += descr_stats_labels[i] + "\t" + str(self.OutputStats['T'][i]) + "\n"
        cp.setText(txt)
        
        #FigureCanvas.__init__(mpltext.fig)
        #FigureCanvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def changeSliderS0(self, value):
        self.label_S0.setText(f"X0 = {log(value)}\nY0 = {value}")
        self.randomwalk_params['S0'] = value

    def changeSliderK(self, value):
        self.label_K.setText(f"log(K) = {log(value)}\nK = {value}")
        self.randomwalk_params['K'] = value

    def changeSliderT1(self, value):
        self.label_t1.setText(str(value))
        self.showGraph()
        
    def changeSliderT2(self, value):
        self.label_t2.setText(str(value))
        self.showGraph()

    def changeValue_oSlider_TimeStep(self, value):
        self.randomwalk_params['T'] = value
        self.label_T.setText(f"n = {self.randomwalk_params['T']}")
        self.slider_t1.setMaximum(self.oSlider_TimeStep.value())
        self.slider_t2.setMaximum(self.oSlider_TimeStep.value())
        #print(f"'T' clicked: {value}")
        
        if isinstance(self.S, np.ndarray):
            self.showGraph()

    def changeValue_oSlider_PathsShown(self, value):
        self.randomwalk_params['P'] = value
        self.label_P.setText(f"Paths = {self.randomwalk_params['P']}")
        if isinstance(self.S, np.ndarray):
            self.showGraph()
        
    def changeValue_oSlider_u(self, value):
        self.randomwalk_params['u'] = value
        self.label_u.setText(f"Δx=\nu = {value/100:.2f}")
        
    def changeValue_comboBox_SimSize(self, value):
        self.randomwalk_params['N'] = int(value)

    def changeValue_prob(self, value):
        self.randomwalk_params['p'] = value
        self.label_p.setText(f"p = {value/100:.2f}")
        
    def on_click_ShowRandomWalk(self):
        # print("I'm clicked!")
        S0 = self.randomwalk_params['S0']
        p  = self.randomwalk_params['p']/100
        u  = self.randomwalk_params['u']/100
        d  = 1/u
        T  = self.randomwalk_params['T']
        N  = self.randomwalk_params['N']
        P  = self.randomwalk_params['P']
        s  = 12345
        params = {'S0': S0, 'p': p, 'u': u, 'd': d, 'T': T, 'N': N, 'P': P, 'seed': s}
        self.slider_t1.setMaximum(self.randomwalk_params['T'])
        self.slider_t2.setMaximum(self.randomwalk_params['T'])
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        self.Y = brwm(params)
        self.X = log(self.Y)
        self.showGraph()
        QApplication.restoreOverrideCursor()
    
        
    def showGraph(self):
        #self.oMplCanvas.close()
        mpl = self.oMplCanvas.canvas
        mainplot = mpl.axes[0][0]
        mainplot.clear()
        mainplot.autoscale(enable=True)
        #mainplot.set_adjustable(adjustable='box')
        #mainplot.set_xlim(auto=True)
        mainplot.set_ylim(auto=True)
        P = self.randomwalk_params['P']
        T = self.randomwalk_params['T']
        #TimeStepShown = min(T, ?) 
        if self.radioGRW.isChecked():
            self.S = self.Y
            self.plot_title = f"Geometric Random Walk"
        elif self.radioRW.isChecked():
            self.S = self.X
            self.plot_title = f"Additive Random Walk"
        mainplot.plot(self.S[:T+1, :P])
        mainplot.grid(True)
        mainplot.set_xlabel('time step', fontsize=8, color='brown')
        mainplot.set_ylabel('price', fontsize=8, color='brown')
        mainplot.tick_params(axis = 'both', which = 'major', labelsize = 6)
        divider = make_axes_locatable(mainplot)
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
        axHist.hist(self.S[-1, :N], density=True, bins=101, orientation='horizontal', rwidth=2)
        axHist.yaxis.set_ticks_position("right")
        axHist.xaxis.set_major_formatter(FuncFormatter('{0:.1%}'.format))
        axHist.grid(True)
        axHist.tick_params(axis = 'both', which = 'major', labelsize = 6)
        p = self.oSlider_prob.value()/100
        SimSize = int(self.comboBox_SimSize.currentText())
        mainplot.set_title(self.plot_title + f" $p$ = {p:.2f}", fontsize=12, color='brown')
        #mainplot.set_title(r'$z = \sqrt{x^2+y^2}$', fontsize=14, color='r')
        axHist.set_xlabel('prob', fontsize=8, color='brown')
        #####
        t1 = self.slider_t1.value()
        t2 = self.slider_t2.value()
        mainplot.axvline(x=t1)
        mainplot.axvline(x=t2, color="blue")
        ax3 = mpl.axes[1][0]
        ax3.clear()
        ax3.set_xlabel('price', fontsize=8, color='brown')
        ax3.tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax3.yaxis.set_major_formatter(FuncFormatter('{0:.1%}'.format))
        MaxStep = self.S.shape[0]-1
        ax3.hist(self.S[min(t1+1, MaxStep), :], density=True, bins=101, orientation='vertical', rwidth=2, color='red', alpha=0.5)
        ax3.hist(self.S[min(t2+1, MaxStep), :], density=True, bins=101, orientation='vertical', rwidth=2, color='blue', alpha=0.5)
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
        samp_kurt = kurtosis(self.S[-1, :])
        table.setItem(3, 0, QTableWidgetItem(f"{samp_kurt:.2f}"))
        #table.setVerticalHeaderLabels(None)
        #table.set_visible(False)
        ## Compute E[(X_n - K)^+] 
        T = self.randomwalk_params['T']
        K = self.randomwalk_params['K']
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
                r'$\mathrm{E}(X_n) = ' + f'{samp_mean:.2f}$',
                r'$\mathrm{Var}(X_n) = ' + f'{samp_var:.2f}$',
                r'$\mathrm{Skew}(X_n) = ' + f'{samp_skew:.2f}$',
                r'$\mathrm{Kurt}(X_n) = ' + f'{samp_kurt:.2f}$',
                r'$\mathrm{E}[(X_n - K)^+] = ' + f'{mean1a:.2f}$',
                r'$\mathrm{E}[(K - X_n)^+] = ' + f'{mean2a:.2f}$',
                ]
        #simulated_statistics.extend([r'line 8', r'line 9', r'line 10', r'line 11', r'line 12', r'line 13'])
        ### Get theoretical values
        y0 = self.randomwalk_params['S0']
        n  = self.randomwalk_params['T']
        u  = self.randomwalk_params['u']/100.0
        p  = self.randomwalk_params['p']/100.0
        print(f"y0={y0}, n={n}, u={u}, p={p}")
        theo_mean = (y0/u**n) * A(2, u, n, p)**n
        theo_var  = (y0/u**n)**2 * ( A(4, u, n, p)**n - A(2, u, n, p)**(2*n) )
        theo_skew = ( (y0/u**n)**3 * ( A(6, u, n, p)**n - 3 * A(4, u, n, p)**n * A(2, u, n, p)**n + 2*A(2, u, n, p)**(3*n) ) ) / theo_var**(3/2)
        theo_kurt = (y0/u**n)**4 * ( A(8, u, n, p)**n - 4 * A(6, u, n, p)**n * A(2, u, n, p)**n + 6 * A(4, u, n, p)**n * A(2, u, n, p)**(2*n) - 3 * A(2, u, n, p)**(4*n) ) / theo_var**2
        theo_mean1a = 0
        theo_mean2a = 0
        analytical_statistics = ['  ', '$Analytical$:',
                r'$\mathrm{E}(X_n) = ' + f'{theo_mean:.2f}$',
                r'$\mathrm{Var}(X_n) = ' + f'{theo_var:.2f}$',
                r'$\mathrm{Skew}(X_n) = ' + f'{theo_skew:.2f}$',
                r'$\mathrm{Kurt}(X_n) = ' + f'{theo_kurt:.2f}$',
                r'$\mathrm{E}[(X_n - K)^+] = ' + f'{theo_mean1a:.2f}$',
                r'$\mathrm{E}[(K - X_n)^+] = ' + f'{theo_mean2a:.2f}$',
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
    sys.exit(app.exec_())
