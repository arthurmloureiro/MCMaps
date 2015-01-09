#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
------------
Arthur E. da Mota Loureiro
------------
"""
from __future__ import print_function
import numpy as np
import grid3D as gr
import sys
from time import time
from scipy import interpolate
import pylab as pl
from matplotlib import cm
import fkp_class as fkpc
import gauss_pk_class as pkgauss
import os
import emcee
import triangle
import uuid
from input import *
from selection_function import *
#####################################################
# the input.py file contains the initial information
#####################################################
path = os.path.expanduser(camb_path)
Map_flat = np.loadtxt(map_file)
Map = Map_flat.reshape(n_x,n_y,n_z)

L_x = n_x*cell_size ; L_y = n_y*cell_size ; L_z = n_z*cell_size 		     # size of the box
box_vol = L_x*L_y*L_z								     # Box's volume
L_max = np.sqrt(L_x*L_x + L_y*L_y + L_z*L_z)	

print("Generating the k-space Grid...\n")
grid = gr.grid3d(n_x,n_y,n_z,L_x,L_y,L_z)					     # generates the grid
grid_bins = gr.grid3d(num_bins, num_bins, num_bins, L_x,L_y,L_z)		     # generates the bins grid
# multiplying the grid for the cell_size will give us a grid in physical units 

##########################################
#	Generating the Bins Matrix M^a_{ijl}
##########################################

nn = int(np.sqrt(n_x**2 + n_y**2 + n_z**2))
kk_bar = np.fft.fftfreq(nn)
kmaxbar = np.sqrt(1.7)*np.max(abs(kk_bar))
dk_bar = kmaxbar/(num_bins+0.0001)*np.ones(num_bins)
k_bar = (dk_bar/2)+np.linspace(0.,kmaxbar-dk_bar[2],num_bins)
M = np.asarray([0.5*(np.sign((k_bar[a]+dk_bar[a]/2)-grid.grid_k[:,:,:n_z/2+1])+1.)*0.5*(np.sign(grid.grid_k[:,:,:n_z/2+1]-(k_bar[a]-dk_bar[a]/2))+1.)for a in range(len(k_bar))])
#sys.exit(-1)
################################################################
#	Assuming a Fiducial selection function n(r) = n_0 exp(-r/b)
################################################################
n_bar_matrix_fid = n_bar_func(grid.grid_r, 8.0,0.04)
#########################################
#	FKP of the data to get the P_data(k)
#########################################
fkp_stuff = fkpc.fkp_init(num_bins,n_bar_matrix_fid,bias,cell_size,n_x,n_y,n_z,M)
FKP1 = fkp_stuff.fkp(Map)
P_data = fkp_stuff.P_ret.real
Sigma_data = fkp_stuff.sigma.real/10 #this 1./10 factor is a problem!

#kk,pp = np.loadtxt("fid_new_matterpower.dat", unpack=1)
#pl.figure()
#pl.loglog()
#pl.xlim(1E-3,2.0)
#pl.ylim(5E2,1E5)
#pl.plot(kk,pp)
#pl.errorbar(k_bar*(2.*np.pi*n_x/L_x),P_data,yerr=Sigma_data,label='Mean Estimated P(k)')
#pl.axvline(x=np.max(np.max(grid.k_x)*(2.*np.pi*n_x/L_x)),color='k',label='Nyquist f.')
#pl.show()
#sys.exit(-1)

##############################################################
#	Let us difine the **p.d.f.** s and some useful functions:
##############################################################
def A_k(P_):									 
	#################################################################################
	# The Gaussian Amplitude
	# Zero Medium and STD=SQRT(2*P(k)*Volume)
	# It must have the 2 factor to take the complex part into account after the iFFT
	#################################################################################
	return np.random.normal(0.0,np.sqrt(2.*P_*box_vol))		
def phi_k(P_): 									     	
	######################
	# Random regular phase
	######################
	return (np.random.random(len(P_)))*2.*np.pi	
def delta_k_g(P_):								     
	########################################	
	# The density contrast in Fourier Space
	########################################
	return A_k(P_)*np.exp(1j*phi_k(P_))	
def delta_x_ln(d_,sigma_,bias_):
	###############################
	# The log-normal density field
	###############################
	return np.exp(bias_*d_ - ((bias_**2)*(sigma_))/2.0) -1.

##############################################################################
#	This is the P_theory function, it will call camb and run GaussianPoisson 
#	realizations of the density field. Then it fkp the final galaxy maps
#	resulting in a avarage fkp power spectrum
##############################################################################
def P_theory(q):
	"""
	q = parameters
	o_cdm = q[0]
	hubble1 = q[1]
	w1 = q[2]
	a = q[3]
	b = abs(q[4])
	c=q[5]
	"""
    #a = n_bar
    #b = 0.01
	##############################
	# modifing camb and running it
	##############################
	#path = os.path.expanduser("~/Documents/camb/params.ini")
	#numb2 = str(np.random.randint(0,250000))
	numb = uuid.uuid4()
	#numb2 = uuid.uuid4()
	#numb=str(1)
	powername = "realiz"+str(numb)
	new_file="realiz"+str(numb)+"_params.ini"
	os.system('touch '+new_file)
	params_file = open("params_realiz.ini")
	linhas = params_file.readlines()
	params_file.close()
	temp=linhas
	ct_p = 0
	temp[0] = "output_root = %s\n" %(powername)
	init = time()
	if hubble[0]==True:
		temp[11] = "hubble         = %.4f\n" %(q[ct_p])
		ct_p = ct_p + 1
	if omega_lambda[0] == True:
		temp[16] = "omega_lambda         = %.4f\n" %(q[ct_p])
		ct_p = ct_p + 1
	if omega_cdm[0] == True:
		temp[15] = "omega_cdm      = %.4f\n" %(q[ct_p])
		ct_p = ct_p + 1
	if omega_baryon[0] == True:
		temp[14] = "omega_baryon   = %.4f\n" %(q[ct_p])
		ct_p = ct_p + 1	
	if omega_neutrino[0] == True:
		temp[17] = "omega_neutrino = %.4f\n" %(q[ct_p])
		ct_p = ct_p + 1
	if w[0] == True:
		temp[12] = "w              = %4f\n" %(q[ct_p])
		ct_p = ct_p + 1
	if n_s[0] == True:
		temp[30] = "scalar_spectral_index(1)           = %4f\n" %(q[ct_p])
		ct_p = ct_p + 1
	if tau[0] == True:
		temp[39] = "re_optical_depth           = %4f\n" %(q[ct_p])
		ct_p = ct_p + 1
	final = time()
	tempo = final - init
	#print("time final P Theory=  ", tempo)
	a = q[ct_p]
	b = q[ct_p+1]
    #temp[29] = " scalar_amp(1)             = %.2e\n" %(A_s)
	out=open(new_file, 'w')
	wtimein = time()
	for i in temp:
		out.write(i)
	out.close()
	wtimef = time()
	wtimet = wtimef - wtimein
	#print("Tempo para escrever arquivo do camb: ", wtimet)
	tcambi = time()
	os.system(path+"/camb " + new_file + " 1>o.txt 2>e.txt")
	tcambf = time()
	tcambtotal = tcambf - tcambi
	#print("Tempo para rodar o CAMB = ", tcambtotal)
	#os.system("~/Documents/Dropbox/Mestrado/camb/camb " + new_file + " 1>o.txt 2>e.txt")
	#############################
	# Reading new camb file
	#############################
	filepower = powername+"_matterpower.dat"
	k_camb, Pk_camb = np.loadtxt(filepower, unpack=1)
	k_camb = k_camb[65:]
	Pk_camb = Pk_camb[65:]
	#os.system("rm " + new_file + " " + filepower)
	final = time()
	#print("tempo camb="+str(final-init))
	#################################
	# Calculating the P_Gauss(k) grid
	#################################
	pkg = pkgauss.gauss_pk(k_camb,Pk_camb,grid.grid_k,cell_size,L_max)
	k_flat = grid.grid_k.flatten()*(2.*np.pi*n_x/L_x)        #this norm factor is so we can have physical unities
	Pk_flat = pkg.Pk_gauss_interp(k_flat)
	p_matrix = Pk_flat.reshape((n_x,n_y,n_z))
	p_matrix[0][0][0] = 1. 						     # Needs to be 1.
	###########################################
	# Generating the selection function matrix
	###########################################
	n_bar_matrix = n_bar_func(grid.grid_r,a,b)
	fkp_mcmc = fkpc.fkp_init(num_bins,n_bar_matrix,bias,cell_size,n_x,n_y,n_z, M)
	######################################
	# here goes the realization's loop
	######################################
	P_all = np.zeros((num_bins, num_realiz))
	sigma_all = np.zeros((num_bins, num_realiz))
	timeloopi = time()
	if realiz_type == 1:
		#print "Doing both Gaussian + Poissonian realizations... \n"
		for m in range(num_realiz):
				#########################
				# gaussian density field
				#########################
				delta_x_gaus = ((delta_k_g(p_matrix).size)/box_vol)*np.fft.ifftn(delta_k_g(p_matrix))	#the iFFT
				var_gr = np.var(delta_x_gaus.real)
				#var_gi = np.var(delta_x_gaus.imag)
				delta_xr_g = delta_x_gaus.real
				#delta_xi_g = delta_x_gaus.imag
				###########################
				# Log-Normal Density Field
				###########################
				delta_xr = delta_x_ln(delta_xr_g, var_gr,bias)
				#delta_xi = delta_x_ln(delta_xi_g, var_gi,bias)
				#######################
				#poissonian realization
				#######################
				N_r = np.random.poisson(n_bar_matrix*(1.+delta_xr))#*(cell_size**3.))			     # This is the final galaxy Map
				#N_i = np.random.poisson(n_bar_matrix*(1.+delta_xi))#*(cell_size**3.))
				##########################################
				#                  FKP                   #
				##########################################
				#print m
				FKP2 = fkp_mcmc.fkp(N_r)
				P_all[:,m] = fkp_mcmc.P_ret.real
				sigma_all[:,m] = fkp_mcmc.sigma.real
				#print m
		
		#print "\nDone.\n"
	elif realiz_type == 2:
		#print "Doing Poissonian realizations only \n"
		#########################
		# gaussian density field
		#########################
		delta_x_gaus = ((delta_k_g(p_matrix).size)/box_vol)*np.fft.ifftn(delta_k_g(p_matrix))	#the iFFT
		var_gr = np.var(delta_x_gaus.real)
		#var_gi = np.var(delta_x_gaus.imag)
		delta_xr_g = delta_x_gaus.real
		#delta_xi_g = delta_x_gaus.imag
		###########################
		# Log-Normal Density Field
		###########################
		delta_xr = delta_x_ln(delta_xr_g, var_gr,bias)
		#delta_xi = delta_x_ln(delta_xi_g, var_gi,bias)
		for m in range(num_realiz):
                
				N_r = np.random.poisson(n_bar_matrix*(1.+delta_xr))#*(cell_size**3.))			     # This is the final galaxy Map
				#N_i = np.random.poisson(n_bar_matrix*(1.+delta_xi))#*(cell_size**3.))
                
				FKP = fkp_mcmc.fkp(N_r)
                
				P = fkp_mcmc.P_ret.real
				sigma = fkp_mcmc.sigma.real
				P_all[:,m] = P
				sigma_all[:,m] = sigma
        		#print "\nDone.\n" 
	else:
		print("Error, invalid option for realization's type \n")
		sys.exit(-1)
	timeloopf = time()
	timeloopt = timeloopf - timeloopi
	#print("Tempo dos Loops = ", timeloopt)
	P_av = (1./num_realiz)*np.sum(P_all, axis=1)
	P_sig = np.sqrt((1./num_realiz)*(np.sum(P_all**2, axis=1))-P_av**2)
	#P_sig =1.
	P_avsig = (1./num_realiz)*np.sum(sigma_all, axis=1)
	os.system('rm ' + powername +"*")
	return P_av.real, P_sig.real, P_avsig.real

########################################################
#	Baysian functions: prior, likelihood and posterior
########################################################
def ln_gaussian(mean, des,x):
	return -0.5*((x-mean)/des)**2

def ln_prior(q):
	cc = 0
	pH,pLamb,pCDM,pBaryon,pNu,pW,pN_s,pTau,pN_bar,pBB=0,0,0,0,0,0,0,0,0,0
	if hubble[0]==True:
		if 40. < q[cc] < 95.:
			pH = ln_gaussian(hubble[1],hubble[2],q[cc])
		else: 
			return -np.inf
		cc = cc + 1
	if omega_lambda[0]==True:
		if  0.0 < q[cc] < 1.:
			pLamb = ln_gaussian(omega_lambda[1],omega_lambda[2],q[cc])
		else: 
			return -np.inf
		cc = cc + 1
	if omega_cdm[0]==True:
		if  0.05 < q[cc] < 0.6:
			pCDM = ln_gaussian(omega_cdm[1],omega_cdm[2],q[cc])
		else: 
			return -np.inf
		cc = cc + 1
	if omega_baryon[0]==True:
		if  0.006 < q[cc] < 0.1:
			pBaryon = ln_gaussian(omega_baryon[1],omega_baryon[2],q[cc])
		else: 
			return -np.inf
		cc = cc + 1
	if omega_neutrino[0]==True:
		if  0. < q[cc] < 0.02:
			pNu = ln_gaussian(omega_neutrino[1],omega_neutrino[2],q[cc])
		else: 
			return -np.inf
		cc = cc + 1
	if w[0]==True:
		if  -2.0 < q[cc] < -0.3333:
			pW = ln_gaussian(w[1],w[2],q[cc])
		else: 
			return -np.inf
		cc = cc + 1
	if n_s[0]==True:
		if  0.80 < q[cc] < 1.20:
			pN_s = ln_gaussian(n_s[1],n_s[2],q[cc])
		else: 
			return -np.inf
		cc = cc + 1
	if tau[0]==True:
		if  0.04 < q[cc] < 0.2:
			pTau = ln_gaussian(tau[1],tau[2],q[cc])
		else: 
			return -np.inf
		cc = cc + 1
	if 4. < q[cc] < 12.:
		pN_bar = ln_gaussian(n_bar0[0],n_bar0[1],q[cc])
		cc = cc +1
	else:
		return -np.inf
	if 0.005 < q[cc] < 0.08:
		pBB = ln_gaussian(bb[0],bb[1],q[cc])
	else:
		return -np.inf
	return pH + pLamb + pCDM + pBaryon + pNu + pW + pN_s + pTau + pN_bar + pBB
def ln_likelihood(q,P_d,sig_d):
	"""
	Defining gaussian likelihood
	"""
	#o_cdm, hubble, a,b = q
	#o_cdm, hubble1, w1, a, b = q
	inicial = time()
	lp = ln_prior(q)
	final2 = time()
	tempo2 = final2 - inicial
	#print("Tempo do Prior = ", tempo2)
	if not np.isfinite(lp):
		return -np.inf
	else:
		ptimei = time()
		theory = P_theory(q)
		ptimef = time()
		ptimet = ptimef - ptimei
		#print("Tempo total do P_theory = ", ptimet)
		varr = sig_d**2 + theory[2]**2
		lk = -0.5*np.sum(((P_d-theory[0])**2)/(varr+1e-20))*(1./num_bins) #ALTEREI AQUI!
		return lk
def ln_post(q,P_d,sig_d):
	"""
	Posterior to be sampled
	"""
	return ln_prior(q) + ln_likelihood(q,P_d,sig_d).real

################################
#	MCMC using the emcee code
################################
#med = np.array([0.25,70., -1.0, 8.0,0.01])	#medium of each inicial step
#desv=np.array([0.03,9., 0.4,1.,0.001])	#deviation of each initial step
med=([])
desv=([])
if hubble[0]==True:
	med.append(hubble[1])
	desv.append(hubble[2])
if omega_lambda[0]==True:
	med.append(omega_lambda[1])
	desv.append(omega_lambda[2])
if omega_cdm[0]==True:
	med.append(omega_cdm[1])
	desv.append(omega_cdm[2])
if omega_baryon[0]==True:
	med.append(omega_baryon[1])
	desv.append(omega_baryon[2])
if omega_neutrino[0]==True:
	med.append(omega_neutrino[1])
	desv.append(omega_neutrino[2])
if w[0]==True:
	med.append(w[1])
	desv.append(w[2])
if n_s[0]==True:
	med.append(n_s[1])
	desv.append(n_s[2])
if tau[0]==True:
	med.append(tau[1])
	desv.append(tau[2])
med.append(n_bar0[0])
med.append(bb[0])
desv.append(n_bar0[1])
desv.append(bb[1])
med = np.array(med) ;  desv = np.array(desv)
starting_guesses = np.zeros((nwalkers,ndim))	#this generates random initial steps 

for i in range(nwalkers):
    starting_guesses[i] = np.random.normal(med,desv)
print(starting_guesses)

init = time()
sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_post, args=[P_data, Sigma_data], threads=ncores)
chain_name_file= chain_name + ".dat"
f = open(chain_name_file, "w")
f.close()

#for result in sampler.sample(starting_guesses, iterations=nsteps, storechain=True):
#for result in sampler.run_mcmc(starting_guesses, N=200):
#	position = np.array(result[0])
#	lnlike   = np.array(result[1])
#	f = open(chain_name_file, "a")
#	for k in range(nwalkers):
#		for d in range(ndim):
#			print(position[k][d], sep=" ", end=" ", file=f)
#		print(lnlike[k], file=f)
        #f.write("{0:4d} {1:s}\n".format(k, " ".join(str(position[k]))))
#	f.close()

sampler.run_mcmc(starting_guesses, N=nsteps)
print("done")
final = time()
print("tempo = "+str(final-init))

print("Mean acceptance fraction: {0:.3f}" .format(np.mean(sampler.acceptance_fraction)))
#samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
samples = sampler.flatchain
np.savetxt(chain_name_file, samples)
#fig = triangle.corner(samples,labels=["$H_0$","$\Omega_{cdm}$",'$n_0$', 'b'], truths=[72.,0.2538,-1., 8., 0.04])
#fig.savefig("fig_"+chain_name+".png")
os.system('rm realiz*')
