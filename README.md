MCMaps
======
=
Author: Arthur Eduardo da Mota Loureiro
Nov/2014
---

Table of content
---
- **mcmaps.py** - The main code, it uses *emcee* to estimate cosmological parameters and selection function parameters
- **MCMaps.ipynb** - Ipython Notebook of an old version of the code. *It needs update*
- **grid3D.py** - An updated version of this class, appearing also in other repos but now it does not have any physical unities. Gives us grids in real and Fourier space
- **gauss_pk_class.py** - a class that calculates the gaussian power spectrum. *It needs better comments*
- **fkp_class.py** - A class that calculates the fkp power spectrum
- **se_der_merda_plota.py** - If something goes wrong, this small code plots what is left of the chains. It can also be used to check whot the chains are.
- **input.dat** - Main input file that contains the map to be analyzed 
- **HighLExtrapTemplate_lenspotentialCls.dat** - CAMB needs this file

Recent Modifications:
---
- 25/11/2014 - Repository created, it still need modifications on the *fkp_class* in order to avoid unphysical power spectra
