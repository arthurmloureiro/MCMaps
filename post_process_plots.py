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
#p=np.loadtxt("H0_CDM_Selec_64_cz20_r1_giant2.dat", unpack=1)
tudo = ([])
for i in range(len(p)-1):
	tudo = np.append(tudo,p[i][nburn:])
N = len(p[0][nburn:])
tudo = tudo.reshape([ndim,N]).T
print 'mean H0 = ' + str(np.mean(p[0][nburn:])) + ' std = ' + str(np.std(p[0][nburn:]))
print 'mean Omega_cdm = ' + str(np.mean(p[1][nburn:])) + ' std = ' + str(np.std(p[1][nburn:]))


#fig = triangle.corner(tudo, labels=['$H_0$','$\Omega_{\Lambda}$','$\Omega_{cdm}$', '$\Omega_{b}$','$\Omega_{nu}$','$w$','n_s','$Tau$','$n_0$', 'b'], truths = [72,0.7,0.25,0.05,0.0,-1.,0.96,0.09,8,0.04],bins=15,color='k')
#fig = triangle.corner(tudo, labels=['$H_0$','$\Omega_{\Lambda}$','$\Omega_{cdm}$', '$\Omega_{b}$','$w$','$n_0$', 'b'], truths = [72,0.7,0.25,0.05,-1.,8.,0.04],bins=15,color='k')
fig = triangle.corner(tudo, labels=['$H_0$','$\Omega_{cdm}$','$n_0$','b'], truths = [72,0.2538,8.,0.05],bins=25,color='k')
fig_name= "fig_"+chain_name+".png"
fig.savefig(fig_name)
pl.show()
