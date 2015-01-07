#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Post process the mcmc chains and plot the resulting data
"""
import numpy as np 
import triangle
import pylab as pl

p = np.loadtxt('4_parameters_chain.dat', unpack=1)
#o1,h1,a1,b1,l1 = np.loadtxt('error_chain200_4p_50w.dat', unpack=1)
tudo = ([])
#tudo1 = ([])
burn=0
tudo = np.append(tudo,p[0][burn:])
tudo = np.append(tudo,p[1][burn:])
tudo = np.append(tudo,p[2][burn:])
tudo = np.append(tudo,p[3][burn:])
#tudo = np.append(tudo,p[4][burn:])
#tudo = np.append(tudo,p[5][burn:])
#tudo = np.append(tudo,p[6][burn:])
#tudo = np.append(tudo,p[7][burn:])
#tudo = np.append(tudo,p[8][burn:])
#tudo = np.append(tudo,p[9][burn:])
#tudo1 = np.append(tudo1,o1)
#tudo1 = np.append(tudo1,h1)
#tudo1 = np.append(tudo1,a1)
#tudo1 = np.append(tudo1,b1)

N = len(p[0][burn:])
#N1= len(o1)

tudo = tudo.reshape([4,N]).T
#tudo1 = tudo1.reshape([4,N1]).T
#print 'mean Omega_cdm = ' + str(np.mean(o[burn:])) + ' std = ' + str(np.std(o[burn:]))
#print 'mean H0 = ' + str(np.mean(h[burn:])) + ' std = ' + str(np.std(h[burn:]))
#print 'mean w = ' + str(np.mean(w[burn:])) + ' std= ' + str(np.std(w[burn:]))
#print 'mean a = ' + str(np.mean(a)) + ' std= ' + str(np.std(a))
#print 'mean b = ' + str(np.mean(b)) + ' std= ' + str(np.std(b))


#fig = triangle.corner(tudo, labels=['$H_0$','$\Omega_{\Lambda}$','$\Omega_{cdm}$', '$\Omega_{b}$','$\Omega_{nu}$','$w$','n_s','$Tau$','$n_0$', 'b'], truths = [72,0.7,0.25,0.05,0.0,-1.,0.96,0.09,8,0.04],bins=15,color='k')
#fig = triangle.corner(tudo, labels=['$H_0$','$\Omega_{\Lambda}$','$\Omega_{cdm}$', '$\Omega_{b}$','$w$','$n_0$', 'b'], truths = [72,0.7,0.25,0.05,-1.,8.,0.04],bins=15,color='k')
fig = triangle.corner(tudo, labels=['$H_0$','$\Omega_{cdm}$','$n_0$', 'b'], truths = [72,0.25,8.,0.04],bins=15,color='k')
fig.savefig('fig_4_params.png')
pl.show()
#fig1 = triangle.corner(tudo1)
#fig1.savefig('compara2.png')
#pl.figure()
#pl.title("$\Omega_{cdm}$")
#pl.hist(o,bins=50,normed=1,label="Fixed window")
#pl.hist(o1,bins=50,alpha=0.8,normed=1,label="marginalized window")
#pl.legend()
#pl.figure()
#pl.title("$H_0$")
#pl.hist(h,bins=50,normed=1,label="Fixed window")
#pl.hist(h1,bins=50,alpha=0.8,normed=1,label="marginalized window")
#pl.legend()
#pl.show()
