########################################
# This is the selection function module
########################################
import numpy as np
def selection_func(gridr,a_,b_):
	"""
	number of parameters = 2
	"""
	return a_*np.exp(-b_*gridr)
#initial guesses for the selection function parameters
n_bar0 = [False,8.0,1.0]
bb = [False,0.04,0.01]
