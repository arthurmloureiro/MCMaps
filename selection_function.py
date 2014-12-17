########################################
# This is the selection function module
########################################
import numpy as np
def n_bar_func(gridr,a_,b_):
	"""
	number of parameters = 2
	"""
	return a_*np.exp(-b_*gridr)
#initial guesses for the selection function parameters
n_bar0 = [7.0,0.8]
bb = [0.01,0.01]
