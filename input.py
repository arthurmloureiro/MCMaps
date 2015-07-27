#! /usr/bin/env python
# -*- coding: utf-8 -*-
#################################################
# This is the input file for the MCMaps code    #
#################################################

# path to camb
camb_path = "~/Documents/camb/"
#camb_path = "~/Documents/Dropbox/Mestrado/camb"
# Map File name - The map to be analyzed 
map_file = "128Gauss_All_New1_cs_10.npy"

#Chain file name.
chain_name = "H0_CDM_w_wa_selec_Gauss128_ALL_NEW"

# Cell physical size in Mpc h^-1
cell_size = 10.00

#Number of cells in x,y and z directions. It must be integer.
n_x = 256
n_y = 256
n_z = 256

# Bias
bias = 1.000

# Number of bins to estimate P(k). It must be integer.
num_bins = 50

################################################
# Which kind of realizations to vary? 
#(1) Gaussian + Poissonian (2) Poissonian Only
################################################
realiz_type = 1

#Number of realizations. It must be integer.
num_realiz = 1

#Number of Parameters to estimate. Cosmological + Selection Function
ndim = 7

#Number of Walkers. It must be integer, bigger than 2*ndim and even
nwalkers = 64

#Number of steps for each walker
nsteps = 1000

#Burn in - Chain points to be thrown away
nburn = 20000

#Number of cores:
ncores = 4

##########################
#Cosmological parameters 
##########################
#[Boolean(if False, assumes the fixed fiducial value), Initial Guess, lower prior, higher prior]:
hubble = [True, 65., 40., 95.,]
omega_lambda = [False, 0.7,0.0,1.0]
omega_cdm = [True, 0.24, 0.05,0.6]
omega_baryon = [False, 0.0462,0.006,0.1]
omega_neutrino = [False, 0.001, 0., 0.02]
w = [True, -1.0, -4.0,-0.453]
w_a = [True, 0.1, -4.0, 0.75]
n_s = [False, 0.96, 0.80,1.20]
tau = [False, 0.09, 0.04,0.2]

#initial guesses for the selection function parameters
n_bar0 = [True,30.,0.01,12]
bb = [True,62.,30.,80.]
c2 = [True,38,20.,45.]
k0 = [False,1.3,0.5,2.]
