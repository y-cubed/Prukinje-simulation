Hituri K, Linne M-L (2013) Comparison of Models for IP3 Receptor
Kinetics Using Stochastic Simulations. PLoS ONE 8(4):
e59618. doi:10.1371/journal.pone.0059618

We provide scripts for the four models compared in our study:

Othmer HG, Tang Y. Oscillations and waves in a model of
InsP3-controlled calcium dynamics, London: Plenum Press, volume 259 of
Experimental and Theoretical Advances in Biological Pattern
Formation. pp. 277-300, 1993

- Dawson A, Lea E, Irvine R. Kinetic model of the inositol
  trisphosphate receptor that shows both steady-state and quantal
  patterns of Ca2+ release from intracellular stores. Biochem J
  370:621, 2003

- Fraiman D, Dawson SP (2004) A model of IP3 receptor with a luminal
calcium binding site: stochastic simulations and analysis. Cell
Calcium 35: 403-413, 2004

- Doi T, Kuroda S, Michikawa T, Kawato M (2005) Inositol
1,4,5-trisphosphate-dependent Ca2+ threshold dynamics detect spike
timing in cerebellar Purkinje cells. J Neurosci 25: 950-961, 2005

The work described in the study was done with STEPS version 1.1.2. on
a Linux computer and has not been tested on newer versions.
For more information about STEPS, please visit
http://steps.sourceforge.net/STEPS

The scripts require that you also have NumPy installed.  For more
information about NumPy, please visit http://www.numpy.org/

The Figures in our publication were plotted with MATLAB. The provided
MATLAB script (run_IP3R_P0.m) will run the simulations and plot the
Figure 2A in the our publication. The simulations will take several
hours. Simulation results are stored in ip3r_MODELNAME_op_res.dat
files and used Ca2+ concetrations in ip3r_MODELNAME_op_ca_concs.dat
files. Please read run_IP3R_P0.m for more information about running
the simulations and plotting the results.

If you do not have MATLAB available you can run the simulations from
command line:
- go to the folder IP3Rmodels

- run each of the simulations separately (each will take time from 10
  min to 2 h):
	python models/sim_doi.py
	python models/sim_fd.py
	python models/sim_dli.py
	python models/sim_ot.py
- the results will be stored in .dat files and you can use a software
  of your choise to plot them
