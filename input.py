#! /usr/bin/env python
# -*- coding: utf-8 -*-
#################################################
# This is the input file for the pk_fkp.py pack #
#################################################

# path to camb
camb_path = "~/Documents/camb/"

# Map File name - The map to be analyzed 
map_file = "map_fid.dat"

#Chain file name. Saves the chain in this file if everything goes wrong and the final figure
chain_name = "DE_chain"

# Cell physical size in Mpc h^-1
cell_size = 20.00

#Number of cells in x,y and z directions. It must be integer.
n_x = 32
n_y = 32
n_z = 32

#Number of realizations. It must be integer.
num_realiz = 10

# Bias
bias = 1.000

# Number of bins to estimate P(k). It must be integer.
num_bins = 18

################################################
# Which kind of realizations to vary? 
#(1) Gaussian + Poissonian (2) Poissonian Only
################################################
realiz_type = 1

#Number of Parameters to estimate. Cosmological + Selection Function
ndim = 5

#Number of Walkers. It must be integer, bigger than 2*ndim and even
nwalkers = 12

#Number of steps for each walker
nsteps = 10

#Burn in - Chain points to be thrown away
nburn = 1

#Number of cores to paralelize the chains:
ncores = 4
