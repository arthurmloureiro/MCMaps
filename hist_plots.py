#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pylab as pl
import matplotlib.cm as cm 

p1 = np.loadtxt("H0_CDM_w_wa_selec_Gauss128_ALL_NEW.dat", unpack=1)
#nburn=50000
nburn=70000
p1=p1[:,nburn:]
fid2=45.00
fid1=28.00
fid3=8.0
#j=np.where(p1[5]<=0.0125) and np.where(p1[5]>0.01245) and np.where(p1[4]<=3.473) and np.where(p1[4]>3.468)
j1=np.where(p1[4]<=30.) 
j2=np.where(p1[4,j1[0]]>29.)
j3=np.where(p1[5,j2[0]]<=62.0)  
j4=np.where(p1[5,j3[0]]>60.0)
j5=np.where(p1[6,j4[0]]<=38.)
j6=np.where(p1[6,j5[0]]>37.)
#j=np.where(p1[5]<=0.1095) and np.where(p1[5]>0.1090) and np.where(p1[4]<=0.1459) and np.where(p1[4]>0.1455)
#j=np.where(p1[5]<=0.2068) and np.where(p1[5]>0.2063) and np.where(p1[4]<=0.001434) and np.where(p1[4]>0.001429)
#j=np.where(p1[5]<=8.0) and np.where(p1[5]>7.95) and np.where(p1[4]<=28.00) and np.where(p1[4]>27.95) and np.where(p1[6]<=8.0) and np.where(p1[6]>7.95)

p2=p1[:,j6[0]]


medias1 = []
desvios1 = []
maxlike_index1 = np.where(p1[7]==np.amax(p1[7]))
for i in range(len(p1)):
	medias1.append(np.percentile(p1[i,:],50))
	#medias1.append(p1[i,maxlike_index1][0,0])
	desvios1.append([np.percentile(p1[i,:],84) - np.percentile(p1[i,:],50), np.percentile(p1[i,:],50) - np.percentile(p1[i,:],16.)])
	#desvios1.append([np.abs(np.percentile(p1[i,:],84.13) - p1[i,maxlike_index1][0,0]), np.abs(p1[i,maxlike_index1][0,0] - np.percentile(p1[i,:],15.869999999999997))])

medias2 = []
desvios2 = []
maxlike_index2 = np.where(p2[7]==np.amax(p2[7]))
for i in range(len(p2)):
	medias2.append(np.percentile(p2[i,:],50))
	#medias2.append(p2[i,maxlike_index2][0,0])
	desvios2.append([np.percentile(p2[i,:],84.) - np.percentile(p2[i,:],50), np.percentile(p2[i,:],50) - np.percentile(p2[i,:],16)])
	#desvios2.append([np.abs(np.percentile(p2[i,:],84.13) - p2[i,maxlike_index2][0,0]), np.abs(p2[i,maxlike_index2][0,0] - np.percentile(p2[i,:],15.869999999999997))])


pl.figure(figsize=(20,20))
pl.subplot(2,2,1)
#pl.subplot(1,2,1)
pl.grid(1)
pl.title("$H_0$",fontsize=25)
labelH01 = "$H_0 = " + "{0:.3f}".format(medias1[0]) + "^{+" + "{0:.3f}".format(desvios1[0][0])+"}_{-" + "{0:.3f}".format(desvios1[0][1])+" }$"
labelH02 = "$H_0 = " + "{0:.3f}".format(medias2[0]) + "^{+" + "{0:.3f}".format(desvios2[0][0])+"}_{-" + "{0:.3f}".format(desvios2[0][1])+" }$"
pl.hist(p1[0],bins=25, normed=1,color='k',label="Marg. $n(r)$ "+labelH01)
pl.hist(p2[0],bins=25, normed=1,color='#4CDEC3',alpha=0.7, label="Fix. $n(r)$ "+labelH02)
pl.axvline(72., label="True Value", linewidth=2.0, c="r")
pl.xticks(size=20)
pl.yticks(size=20)
pl.legend(loc=9,bbox_to_anchor=(0., 1.05, 1., .102),prop={'size':20},shadow=1)

pl.subplot(2,2,2)
#pl.subplot(122)
pl.grid(1)
pl.title("$\Omega_{cdm}$",fontsize=25)
labelCDM1 = "$\Omega_{cdm} = " + "{0:.4f}".format(medias1[1]) + "^{+" + "{0:.5f}".format(desvios1[1][0])+"}_{-" + "{0:.5f}".format(desvios1[1][1])+" }$"
labelCDM2 = "$\Omega_{cdm} = " + "{0:.4f}".format(medias2[1]) + "^{+" + "{0:.5f}".format(desvios2[1][0])+"}_{-" + "{0:.5f}".format(desvios2[1][1])+" }$"
pl.hist(p1[1],bins=25, normed=1,color='k',label="Marg. $n(r)$ "+labelCDM1)
pl.hist(p2[1],bins=25, normed=1,color='#4CDEC3',alpha=0.7, label="Fix. $n(r)$ " + labelCDM2)
pl.axvline(0.2538, label="True Value", linewidth=2.0, c="r")
pl.xticks(size=20)
pl.yticks(size=20)
pl.legend(loc=9,bbox_to_anchor=(0., 1.05, 1., .102),prop={'size':20},shadow=1)

pl.subplot(2,2,3)
#pl.subplot(1,2,1)
pl.grid(1)
pl.title("$w_0$",fontsize=25)
labelw01 = "$w_0 = " + "{0:.3f}".format(medias1[2]) + "^{+" + "{0:.3f}".format(desvios1[2][0])+"}_{-" + "{0:.3f}".format(desvios1[2][1])+" }$"
labelw02 = "$w_0 = " + "{0:.3f}".format(medias2[2]) + "^{+" + "{0:.3f}".format(desvios2[2][0])+"}_{-" + "{0:.3f}".format(desvios2[2][1])+" }$"
pl.hist(p1[2],bins=25, normed=1,color='k',label="Marg. $n(r)$ " +labelw01 )
pl.hist(p2[2],bins=25, normed=1,color='#4CDEC3',alpha=0.7, label="Fix. $n(r)$ "+labelw02)
pl.axvline(-1.0, label="True Value", linewidth=2.0, c="r")
pl.xticks(size=20)
pl.yticks(size=20)
#pl.legend(loc=0,bbox_to_anchor=(0.35, 1.00),prop={'size':20},shadow=1)
pl.legend(loc=9,bbox_to_anchor=(0., 1.05, 1., .102),prop={'size':20},shadow=1)
pl.subplot(2,2,4)
#pl.subplot(1,2,2)
pl.grid(1)
pl.title("$w_a$",fontsize=25)
labelwa1 = "$w_a = " + "{0:.3f}".format(medias1[3]) + "^{+" + "{0:.3f}".format(desvios1[3][0])+"}_{-" + "{0:.3f}".format(desvios1[3][1])+" }$"
labelwa2 = "$w_a = " + "{0:.3f}".format(medias2[3]) + "^{+" + "{0:.3f}".format(desvios2[3][0])+"}_{-" + "{0:.3f}".format(desvios2[3][1])+" }$"
pl.hist(p1[3],bins=25, normed=1,color='k',label="Marginal. $n(r)$ " + labelwa1)
pl.hist(p2[3],bins=25, normed=1,color='#4CDEC3',alpha=0.7, label="Fixing $n(r)$ "+ labelwa2 )
pl.axvline(0.0, label="True Value", linewidth=2.0, c="r")
pl.xticks(size=20)
pl.yticks(size=20)
pl.legend(loc=9,bbox_to_anchor=(0., 1.05, 1., .102),prop={'size':20},shadow=1)

pl.savefig("fig_hist_Gauss128_ALL_NEW.png",transparent=0)
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
