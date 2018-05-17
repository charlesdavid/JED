## JED

* JED is a powerful tool for examining the dynamics of proteins from trajectories derived from MD or Geometric simulations (or NMA).
* Currently, there are two types of PCA: distance-pair and Cartesian, and three models: COV, CORR, and PCORR.
* The idea behind the development of JED was to not only analyze a trajectory, but also compare it to other dynamical analyses.
* To do this, JED uses subspace analysis to assess the similarity of essential vector spaces defined by each dynamical method.
* Additionally, JED computes the Free Energy of a trajectory from two PCs (weighted DVPs) so that a FE surface can be constructed.
* Finally, JED produces PDB frames and script files for viewing movies in PyMol(TM) of both individual modes and a superposition of Essential Modes.


### This repo is for the developers of JED to collaborate...
### and for users of JED to submit feedback

* Please report bugs by submitting a new issue tagged with BUG:  
	* Be sure to provide enough data to reproduce the error  
	* Include JED input files, JED errors, and Java stack traces  

* Please suggest new features by submitting a new issue tagged with NEW FEATURE  

* Please contact me if you wish to become a collaborator on this software:  

Charles David: cdavid@carolina.rr.com  

##### Last Update: May 2018
* New JavaDocs
* New Jars
* Minor code fixes
* Added FES to JED Driver
* Created new input files
* Fixed Pymol script output (now works with educational license)
