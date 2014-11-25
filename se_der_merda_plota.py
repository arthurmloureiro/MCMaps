"""
Use se der merda nas cadeias e plote mesmo assim hue lololol brbrbr
"""
import numpy as np 
import triangle
import pylab as pl

o,h,w,a,b,l = np.loadtxt('DE_chain300_p5_w50.dat', unpack=1)
#o1,h1,a1,b1,l1 = np.loadtxt('error_chain200_4p_50w.dat', unpack=1)
tudo = ([])
#tudo1 = ([])
burn=7000
tudo = np.append(tudo,o[burn:])
tudo = np.append(tudo,h[burn:])
tudo = np.append(tudo,w[burn:])
tudo = np.append(tudo,a[burn:])
tudo = np.append(tudo,b[burn:])
#tudo1 = np.append(tudo1,o1)
#tudo1 = np.append(tudo1,h1)
#tudo1 = np.append(tudo1,a1)
#tudo1 = np.append(tudo1,b1)

N = len(o[burn:])
#N1= len(o1)

tudo = tudo.reshape([5,N]).T
#tudo1 = tudo1.reshape([4,N1]).T
print 'mean Omega_cdm = ' + str(np.mean(o[burn:])) + ' std = ' + str(np.std(o[burn:]))
print 'mean H0 = ' + str(np.mean(h[burn:])) + ' std = ' + str(np.std(h[burn:]))
print 'mean w = ' + str(np.mean(w[burn:])) + ' std= ' + str(np.std(w[burn:]))
#print 'mean a = ' + str(np.mean(a)) + ' std= ' + str(np.std(a))
#print 'mean b = ' + str(np.mean(b)) + ' std= ' + str(np.std(b))


fig = triangle.corner(tudo, labels=['$\Omega_{cdm}$', '$H_0$','$\omega$', '$n_0$', 'b'], truths = [0.25,72,-1.,8,0.01],bins=15,color='k')
fig.savefig('fig_DE_c300_p5_w50BURN.png')
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
