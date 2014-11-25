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
############################################################
#	Reading the input file, the map and converting unities
############################################################
map_file, cell_size, n_x, n_y, n_z, num_realiz, bias, num_bins, n_bar, realiz_type = np.loadtxt('input.dat', dtype=str)
cell_size = float(cell_size); n_x=int(n_x); n_y=int(n_y); n_z=int(n_z); num_realiz=int(num_realiz); bias=float(bias) ; num_bins=int(num_bins); realiz_type = int(realiz_type); n_bar = float(n_bar);

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
kk_bar = np.fft.fftfreq(n_x)
kminbar = np.min(abs(kk_bar))
kmaxbar = np.sqrt(1.7)*np.max(abs(kk_bar))
dk_bar = kmaxbar/(num_bins+0.0001)*np.ones(num_bins)
k_bar = (dk_bar/2)+np.linspace(0.,kmaxbar-dk_bar[2],num_bins)
M = np.asarray([0.5*(np.sign((k_bar[a]+dk_bar[a]/2)-grid.grid_k[:,:,:n_x/2+1])+1.)*0.5*(np.sign(grid.grid_k[:,:,:n_x/2+1]-(k_bar[a]-dk_bar[a]/2))+1.)for a in range(len(k_bar))])

################################################################
#	Assuming a Fiducial selection function n(r) = n_0 exp(-r/b)
################################################################
def n_bar_func(gridr,a_,b_):
    return a_*np.exp(-b_*gridr)
n_bar_matrix_fid = n_bar_func(grid.grid_r, n_bar,0.01)
#########################################
#	FKP of the data to get the P_data(k)
#########################################
fkp_stuff = fkpc.fkp_init(num_bins,n_bar_matrix_fid,bias,cell_size,n_x,M)
FKP = fkp_stuff.fkp(Map)
P_data = fkp_stuff.P_ret.real
Sigma_data = fkp_stuff.sigma.real/10 #this 1./10 factor is a problem!
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
    """
    o_cdm = q[0]
    hubble = q[1]
    if hubble > 100.00:
    	hubble = 100
    elif hubble < 40.00:
    	hubble = 40
    w = q[2]
    if w > 0.05:
    	w=0.05
    elif w < -1.20:
    	w=-1.20
    #A_s = q[2]
    #if A_s > 5.e-9:
    # 	A_s = 5e-9
    #elif A_s < 0.2e-9:
    #	A_s = 0.2e-9
    a = q[3]
    b = abs(q[4])
    #a = n_bar
    #b = 0.01
    init = time()
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
    temp[0] = "output_root = %s\n" %(powername)
    temp[11] = "hubble         = %.4f\n" %(hubble)
    temp[12] = "w              = %4f\n" %(w)
    temp[15] = "omega_cdm      = %.4f\n" %(o_cdm)
    #temp[29] = " scalar_amp(1)             = %.2e\n" %(A_s)
    out=open(new_file, 'w')
    for i in temp:
        out.write(i)
    out.close()
    #os.system("~/Documents/camb/camb " + new_file + " 1>o.txt 2>e.txt")
    os.system("~/Documents/Dropbox/Mestrado/camb/camb " + new_file + " 1>o.txt 2>e.txt")
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
    ######################################
    # here goes the realization's loop
    ######################################
    P_all = np.zeros((num_bins, num_realiz))
    sigma_all = np.zeros((num_bins, num_realiz))
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
                FKP = fkp_stuff.fkp(N_r)
                #print m
                P_all[:,m] = fkp_stuff.P_ret
                sigma_all[:,m] = fkp_stuff.sigma
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
                
                FKP = fkp_stuff.fkp(N_r)
                
                P = fkp_stuff.P_ret
                sigma = fkp_stuff.sigma
                P_all[:,m] = P
                sigma_all[:,m] = sigma
                

        #print "\nDone.\n" 
    else:
        #print "Error, invalid option for realization's type \n"
        sys.exit(-1)
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
	return -0.5*((mean-x)/des)**2

def ln_prior(q):
    """
	gaussian prior on the parameters
    """
    o_cdm, hubble, w,a,b = q
    ocdm_mean = 0.25
    ocdm_des = 0.025
    hubble_mean	= 72.
    hubble_des = 7.
    #A_mean = 2.1e-9
    #A_des = 2.1e-10
    a_mean = 8.
    a_des = 1.0
    b_mean	= 0.01
    b_des = 0.008
    w_mean = -1.0
    w_des = 0.1
    #o_cdm, hubble = q
    #a=n_bar
    #b=0.01
    #if 0.05 < o_cdm < 0.5 and 40. < hubble < 90. and 3. < a < 12. and 0.005 < b < 0.015:
    #    return 0.0
    #return -np.inf
    pO_cdm = ln_gaussian(ocdm_mean,ocdm_des,o_cdm)
    #if pO_cdm < 0.01:
    #	pO_cdm = 0.01
    pH = ln_gaussian(hubble_mean, hubble_des, hubble)
    pw = ln_gaussian(w_mean,w_des,w)
    #if pH > 90.00:
    #    pH = 90
    pa = ln_gaussian(a_mean, a_des, a)
    pb = ln_gaussian(b_mean, b_des, b)
    if 0.05 < o_cdm < 0.5  and 40.< hubble < 90. and -1.2 < w < 0.00 and 4.< a < 12 and 0.005 < b < 0.02:	
        return pO_cdm + pH + pw + pa + pb
        #return 0.0 
    return -np.inf
    
def ln_likelihood(q,P,sig):
    """
    Defining gaussian likelihood
    """
    #o_cdm, hubble, a,b = q
    o_cdm, hubble, w, a, b = q
    #a=n_bar
    #b=0.01
    #lp = ln_prior(q)
    #if not np.isfinite(lp):
    #   return -np.inf
    theory = P_theory(q)
    varr = sig**2 + theory[2]**2 #+ theory[2]**2
    lk = -0.5*np.sum(((P-theory[0])**2)/(varr+1e-20))*(1./num_bins) #ALTEREI AQUI!
    #print("====> ", P[1:], "\n -----> ", theory[0][1:])
    return lk
def ln_post(q,P,sig):
	"""
	Posterior to be sampled
	"""
	return ln_prior(q) + ln_likelihood(q,P,sig).real
    
################################
#	MCMC using the emcee code
################################
ndim = 5   #number of params
nwalkers = 50   #number of zombies?
nburn = 1     #burn in
nsteps = 300    #number of MCMC steps for each walker 
med = np.array([0.25,70., -1.0, 8.0,0.01])	#medium of each inicial step
desv=np.array([0.03,9., 0.4,1.,0.001])	#deviation of each initial step
starting_guesses = np.zeros((nwalkers,ndim))	#this generates random initial steps 
												#for each walker
for i in range(nwalkers):
    starting_guesses[i] = np.random.normal(med,desv)
print(starting_guesses)

init = time()
sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_post, args=[P_data, Sigma_data], threads=4)
chain_name = "DE_chain300_p5_w50.dat"
f = open(chain_name, "w")
f.close()

for result in sampler.sample(starting_guesses, iterations=nsteps, storechain=False):
	position = np.array(result[0])
	lnlike   = np.array(result[1])
	f = open(chain_name, "a")
	for k in range(nwalkers):
		for d in range(ndim):
			print(position[k][d], sep=" ", end=" ", file=f)
		print(lnlike[k], file=f)
        #f.write("{0:4d} {1:s}\n".format(k, " ".join(str(position[k]))))
	f.close()

sampler.run_mcmc(starting_guesses, nsteps)
print("done")
final = time()
print("tempo = "+str(final-init))

#emcee_trace = sampler.chain[:, nburn:, :].reshape(-1, ndim).T
samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
fig = triangle.corner(samples,labels=["$\Omega_{cdm}$", "$H_0$",'$w$', '$n_0$', 'b'],
                      truths=[0.2538, 72.,-1., 8., 0.01])
fig.savefig("fig_DEc300_p5_w50.png")
#os.system('rm realiz*')