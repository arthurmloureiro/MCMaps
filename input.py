#! /usr/bin/env python
# -*- coding: utf-8 -*-
#################################################
# This is the input file for the pk_fkp.py pack #
#################################################

# path to camb
camb_path = "~/Documents/camb/"
#camb_path = "~/Documents/Dropbox/Mestrado/camb"
# Map File name - The map to be analyzed 
map_file = "map_TESTE.dat"

#Chain file name.
chain_name = "2_params_a15_big"

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
ndim = 2

#Number of Walkers. It must be integer, bigger than 2*ndim and even
nwalkers = 50

#Number of steps for each walker
nsteps = 200

#Burn in - Chain points to be thrown away
nburn = 0

#Number of cores to paralelize the chains:
ncores = 4

##########################
#Cosmological parameters 
##########################
#[Boolean(if False, assumes the fixed fiducial value), Mean, Std for the gaussian prior]:
hubble = [True, 60., 8.0]
omega_lambda = [False, 0.7,0.2]
omega_cdm = [True, 0.20,0.1]
omega_baryon = [False, 0.0462,0.0046]
omega_neutrino = [False, 0.001, 0.01]
w = [False, -1.0, 1.0]
n_s = [False, 0.96, 0.096]
tau = [False, 0.09,0.009]


