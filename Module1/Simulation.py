# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import math

def binomial_randomwalk_multiplicative_simulate(params):
    '''
    Multiplicative binomial random walk simulation
    S0 = 100   # initial price
    u  = 1.1   # "up" factor
    d  = 0.9   # "down" factor, currently set to the multiplicative inverse of u
    p  = 0.5   # probability of "up"
    T  = 30    # time-step size
    N  = 50000 # sample size (no. of simulations)
    '''
    S0 = params['S0']
    p  = params['p']
    u  = params['u']
    d  = 1/u
    T  = params['T']
    N  = params['N']
    P  = params['P']
    s  = params['seed']
    np.random.seed(s)
    # Simulating N paths with T time steps
    S = np.zeros((T + 1, N))
    S[0] = S0
    for t in range(1, T + 1):
        z = np.random.rand(N)  # pseudorandom numbers
        S[t] = S[t - 1] * ( (z<=p)*u + (z>p)*d )
          # vectorized operation per time step over all paths
    return S

def returnSampleDescriptiveStatistics(S, K, decimalPlaces=2):
    '''
    Returns mean, variance, skewness, kurtosis, etc.
    '''
    samp_mean = f"{S[-1, :].mean():.2f}"
    samp_var = f"{S[-1, :].var():.2f}"
    samp_skew = f"{skew(S[-1, :]):.2f}"
    samp_kurt = f"{kurtosis(S[-1, :]):.2f}"
    
    #samp_X_minus_K = f"{S[-1, :].mean():.2f}"
    #samp_K_minus_X = 


#def random_walk

def generate_path(S0, r, sigma, T, M):
	''' 
	Source: Hilpisch 2017, p. 236
	
	Function to simulate a geometric Brownian motion.
	Parameters
	==========
		S0: float
		initial index level
		r: float
		constant risk-less short rate
		sigma: float
		instantaneous volatility
		T: float
		date of maturity (in year fractions)
		M: int
		number of time intervals
	Returns
	=======
		path: pandas DataFrame object simulated path
	'''
	# length of time interval

	dt = float(T) / M
	# random numbers
	np.random.seed(100000)
	rn = np.random.standard_normal(M + 1)
	rn[0] = 0 # to keep the initial value
	# simulation of path
	path = S0 * np.exp(np.cumsum((r - 0.5 * sigma ** 2) * dt + sigma * math.sqrt(dt) * rn))
	# setting initial value
	path = pd.DataFrame(path, columns=['index'])
	return path


