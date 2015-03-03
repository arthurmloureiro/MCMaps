#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pylab as pl
import matplotlib.cm as cm 

p1 = np.loadtxt("H0_CDM_w_wa_selec_big.dat", unpack=1)
nburn=120000
p1=p1[:,nburn:]
#p2 = np.loadtxt("H0_CDM_w_wa.dat", unpack=1)
j=np.where(p1[6]<=1.20) and np.where(p1[6]>1.19) and np.where(p1[5]<=0.20) and np.where(p1[5]>0.19) and np.where(p1[4]<=0.06250) and np.where(p1[4]>0.06249)
p2=p1[:,j[0]]

def compute_sigma_level(chain1, chain2, nbins=20):
	"""From a set of chains, bin by number of standard deviations"""
	L, xbins, ybins = np.histogram2d(chain1, chain2, nbins)
	L[L == 0] = 1E-16
	logL = np.log(L)

	shape = L.shape
	L = L.ravel()

	# obtain the indices to sort and unsort the flattened array
	i_sort = np.argsort(L)[::-1]
	i_unsort = np.argsort(i_sort)

	L_cumsum = L[i_sort].cumsum()
	L_cumsum /= L_cumsum[-1]
    
	xbins = 0.5 * (xbins[1:] + xbins[:-1])
	ybins = 0.5 * (ybins[1:] + ybins[:-1])

	return xbins, ybins, L_cumsum[i_unsort].reshape(shape)


medias1 = []
desvios1 = []
for i in range(len(p1)):
	medias1.append(np.percentile(p1[i,:],50))
	desvios1.append([np.percentile(p1[i,:],84.13) - np.percentile(p1[i,:],50), np.percentile(p1[i,:],50) - np.percentile(p1[i,:],15.869999999999997)])

medias2 = []
desvios2 = []
for i in range(len(p2)):
	medias2.append(np.percentile(p2[i,:],50))
	desvios2.append([np.percentile(p2[i,:],84.13) - np.percentile(p2[i,:],50), np.percentile(p2[i,:],50) - np.percentile(p1[i,:],15.869999999999997)])


pl.figure(figsize=(22,12))
pl.subplot(2,2,1)
#pl.subplot(1,2,1)
pl.grid(1)
pl.title("$H_0$",fontsize=25)
labelH01 = "$H_0 = " + "{0:.3f}".format(medias1[0]) + "^{+" + "{0:.3f}".format(desvios1[0][0])+"}_{-" + "{0:.3f}".format(desvios1[0][1])+" }$"
labelH02 = "$H_0 = " + "{0:.3f}".format(medias2[0]) + "^{+" + "{0:.3f}".format(desvios2[0][0])+"}_{-" + "{0:.3f}".format(desvios2[0][1])+" }$"
pl.hist(p1[0],bins=25, normed=1,label="Marg. $n(r)$ "+labelH01)
pl.hist(p2[0],bins=25, normed=1,alpha=0.7, label="Fix. $n(r)$ "+labelH02)
pl.axvline(72., label="True Value", linewidth=2.0, c="k")
pl.xticks(size=20)
pl.yticks(size=20)
pl.legend(loc=0,bbox_to_anchor=(1.05, 1.05),prop={'size':15},shadow=1)

pl.subplot(2,2,2)
#pl.subplot(122)
pl.grid(1)
pl.title("$\Omega_{cdm}$",fontsize=25)
labelCDM1 = "$\Omega_{cdm} = " + "{0:.3f}".format(medias1[1]) + "^{+" + "{0:.3f}".format(desvios1[1][0])+"}_{-" + "{0:.3f}".format(desvios1[1][1])+" }$"
labelCDM2 = "$\Omega_{cdm} = " + "{0:.3f}".format(medias2[1]) + "^{+" + "{0:.3f}".format(desvios2[1][0])+"}_{-" + "{0:.3f}".format(desvios2[1][1])+" }$"
pl.hist(p1[1],bins=25, normed=1,label="Marg. $n(r)$ "+labelCDM1)
pl.hist(p2[1],bins=25, normed=1,alpha=0.7, label="Fix. $n(r)$ " + labelCDM2)
pl.axvline(0.2538, label="True Value", linewidth=2.0, c="k")
pl.xticks(size=20)
pl.yticks(size=20)
pl.legend(loc=0,bbox_to_anchor=(1.05, 1.05),prop={'size':15},shadow=1)

pl.subplot(2,2,3)
#pl.subplot(1,2,1)
pl.grid(1)
pl.title("$w_0$",fontsize=25)
labelw01 = "$w_0 = " + "{0:.3f}".format(medias1[2]) + "^{+" + "{0:.3f}".format(desvios1[2][0])+"}_{-" + "{0:.3f}".format(desvios1[2][1])+" }$"
labelw02 = "$w_0 = " + "{0:.3f}".format(medias2[2]) + "^{+" + "{0:.3f}".format(desvios2[2][0])+"}_{-" + "{0:.3f}".format(desvios2[2][1])+" }$"
pl.hist(p1[2],bins=25, normed=1,label="Marg. $n(r)$ " +labelw01 )
pl.hist(p2[2],bins=25, normed=1,alpha=0.7, label="Fix. $n(r)$ "+labelw02)
pl.axvline(-1.0, label="True Value", linewidth=2.0, c="k")
pl.xticks(size=20)
pl.yticks(size=20)
pl.legend(loc=0,bbox_to_anchor=(1.05, 1.05),prop={'size':15},shadow=1)

pl.subplot(2,2,4)
#pl.subplot(1,2,2)
pl.grid(1)
pl.title("$w_a$",fontsize=25)
labelwa1 = "$w_a = " + "{0:.3f}".format(medias1[3]) + "^{+" + "{0:.3f}".format(desvios1[3][0])+"}_{-" + "{0:.3f}".format(desvios1[3][1])+" }$"
labelwa2 = "$w_a = " + "{0:.3f}".format(medias2[3]) + "^{+" + "{0:.3f}".format(desvios2[3][0])+"}_{-" + "{0:.3f}".format(desvios2[3][1])+" }$"
pl.hist(p1[3],bins=25, normed=1,label="Marginal. $n(r)$ " + labelwa1)
pl.hist(p2[3],bins=25, normed=1,alpha=0.7, label="Fixing $n(r)$ "+ labelwa2 )
pl.axvline(0.0, label="True Value", linewidth=2.0, c="k")
pl.xticks(size=20)
pl.yticks(size=20)
pl.legend(loc=0,bbox_to_anchor=(0.80, 1.05),prop={'size':15},shadow=1)

pl.savefig("fig_hist.png",transparent=1)
pl.show()
"""
class nf(float):
###############################
# forces to show contours in %
###############################
	def __repr__(self):
		str = '%.1f' % (self.__float__(),)
		if str[-1]=='0':
			return '%.1f' % self.__float__()
		else:
			return '%.1f' % self.__float__()

cmap=cm.jet
level=[0.683,0.955]
x1,y1,sig1 = compute_sigma_level(p1[0][25000:],p1[1][25000:])
pl.figure()
pl.subplot(111)
CS = pl.contour(x1,y1,sig1.T,colors=['b','g'],linewidths=(2.0,1.5,1.0),levels=level)
# Recast levels to new class
CS.levels = [nf(val*100) for val in CS.levels ]
# Label levels with specially formatted floats
if pl.rcParams["text.usetex"]:
     fmt = r'%r \%%'
else:
     fmt = '%r %%'
pl.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=10)

#pl.plot(p1[0],p1[1],'.y',alpha=0.2)
#CS=pl.contourf(x1,y1,sig1.T,levels=level,antialiased=0, cmap=cmap)
#pl.clabel(CS,inline=1)
#pl.legend(CS, ["1$\sigma$", "2$\sigma$"])
pl.show()
"""
