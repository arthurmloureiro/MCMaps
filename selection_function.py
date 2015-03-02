########################################
# This is the selection function module
########################################
import numpy as np
#def selection_func(gridr,a_,b_):
#	"""
#	number of parameters = 2
#	"""
#	return a_*np.exp(-b_*gridr)
def selection_func(x,y,z,n0,c1,c2,c3):
	"""
	number of parameters = 3
	n0 = mean number of galaxies/cell
	c1=b
	c2=c2
	c3=k0
	"""
	c0 = 1.
	part1 = (c0 + c2*np.outer(np.sin(c3*x),np.sin(c3*y)))
	return n0*np.einsum('ij,l->ijl',part1,np.exp(-c1*z))
