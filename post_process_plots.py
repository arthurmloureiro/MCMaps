#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Post process the mcmc chains and plot the resulting data
"""
import numpy as np 
import triangle
import pylab as pl
from input import *

p = np.loadtxt(chain_name+'.dat', unpack=1)

tudo = ([])
label=([])
true_vals = ([])
for i in range(len(p)-1):
	tudo = np.append(tudo,p[i][nburn:])
N = len(p[0][nburn:])
tudo = tudo.reshape([ndim,N]).T
count = 0
if hubble[0] == True:
	label.append("$H_0$")
	true_vals.append(72)
	print 'mean H0 = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if omega_lambda[0] == True:
	label.append("$\Omega_{\Lambda}$")
	true_vals.append(0.7)
	print 'mean Omg_Lamb = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if omega_cdm[0] == True:
	label.append("$\Omega_{cdm}$")
	true_vals.append(0.2538)
	print 'mean Omg_cdm = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if omega_baryon[0] == True:
	label.append("$\Omega_{b}$")
	true_vals.append(0.0462)
	print 'mean Omg_b = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if omega_neutrino[0] == True:
	label.append("$\Omega_{b}$")
	true_vals.append(0.0)
	print 'mean Omg_neu = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if w[0] == True:
	label.append("$w$")
	true_vals.append(-1.0)
	print 'mean w = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if w_a[0] == True:
	label.append("$w_a$")
	true_vals.append(0.0)
	print 'mean w_a = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if n_s[0] == True:
	label.append("$n_s$")
	true_vals.append(0.96)
	print 'mean n_s = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if tau[0] == True:
	label.append(r"$\tau$")
	true_vals.append(0.09)
	print 'mean Tau = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1	
if n_bar0[0] == True:
	label.append(r"$n_0$")
	true_vals.append(30.0)
	print 'mean n0 = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if bb[0] == True:
	#label.append(r"$b$")
	#true_vals.append(1./16)
	label.append(r"$R$")
	true_vals.append(62.0)
	print 'mean c_1 = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if c2[0] == True:
	label.append(r"$\sigma_r$")
	#true_vals.append(0.2)
	true_vals.append(38.0)
	print 'mean sig_r = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1
if k0[0] == True:
	label.append(r"$k_0$")
	true_vals.append(1.2)
	print 'mean k0 = '+ str(np.mean(p[count][nburn:])) + ' std = ' + str(np.std(p[count][nburn:]))
	count = count + 1

fig = triangle.corner(tudo, labels=label, truths =true_vals,scale_hist=1,bins=20,color='k',quantiles=[0.16, 0.5, 0.84],verbose=1,plot_datapoints=0)
#fig.title("using j-pas selection function")
fig.suptitle("     Gaussian Selection function",fontsize=23)
fig_name= "fig_"+chain_name+".png"
fig.savefig(fig_name,transparent=0)
pl.show()
