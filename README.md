#MCMaps
======
##Author: Arthur Eduardo da Mota Loureiro
##Nov/2014
Set up instructions:
---
- **1)** Compile camb without MPI so it can run serial. Compile it without the *-fopenmp* flag in the compiler line of *Makefile*
- **2)** Install *emcee* through *$sudo pip install emcee*
- **3)** Install the *Triangle Plot* pack through *$sudo pip install triangle_plot*
- **4)** Modify in the code where camb is

Table of content
---
- **mcmaps.py** - The main code, it uses *emcee* to estimate cosmological parameters and selection function parameters
- **MCMaps.ipynb** - Ipython Notebook of an old version of the code. *It needs update*
- **grid3D.py** - An updated version of this class, appearing also in other repos but now it does not have any physical unities. Gives us grids in real and Fourier space
- **gauss_pk_class.py** - a class that calculates the gaussian power spectrum. *It needs better comments*
- **fkp_class.py** - A class that calculates the fkp power spectrum
- **post_process_plot.py** - If something goes wrong, this small code plots what is left of the chains. It can also be used to check whot the chains are
- **input.py** - Main input file. Contains information about the map to analyse, the parameters to vary, etc
- **HighLExtrapTemplate_lenspotentialCls.dat** - CAMB needs this file
- **params_realiz.ini** - mcmaps reads this file and uses it to run camb
- **selection_function.py** - File containing information about the selection function and its parameters
Recent Modifications:
---
- 25/11/2014 - Repository created, it still need modifications on the *fkp_class* in order to avoid unphysical power spectra
- 26/11/2014 - Modifications on *fkp_class* made and fixed the *still running after done* problem
- 26/11/2014 Modifications on the input file. Now the MCMC nwalkers, nsteps, etc are a part of the input file
- 27/11/2014 The new *input.py* replaces *input.dat* and goes to the program through *from input import .*
- 17/11/2014 Generalized the code, it can accept any combination of cosmological parameters to prope from {hubble, omega-lambda, omega-cdm, omega-baryon, omega-neutrino, w, n-s, tau}. _5 parameters estimation for 200 points took 226,38 s in 4 cores_
- 07/01/2015 Realized there's a problem with the emcee code, need to
change the _.sample_ to _.runmcmc_ **need to realize who to save the
chains**