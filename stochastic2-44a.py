#stochastic2-44a.py

		#based on the paper "stochastic chemical reactions"
		#model:
		#X1bar ->c1-> Y1


import steps.model as smodel  
		#import package steps.model that contains all the definitions of 
		#the objects and functions you need to describe the physics and 
		#chemistory.

mdl = smodel.Model()  
		#mdl variable for discribing model.
		#create a top-level container object for our model this top model
		#container is required for all simulations in STEPS.

molX1bar = smodel.Spec('molX1bar', mdl)
molY1 = smodel.Spec('molY1',mdl)
		#create 2 steps.model.Spec objects corresponding to 2 chamical
		#spicies

vsys = smodel.Volsys('vsys', mdl)
		#create a volume system
		#volume systems art container objects that group a number of 
		#stoichimetric reaction rules.


c1reac_f = smodel.Reac('c1reac_f', vsys, lhs = [molX1bar], rhs = [molY1], kcst = 0.2)
		#create the reaction rules themselves
		#what is cicst = 0.3e6?

import steps.geom as swm
		#import the geometry module that contains the objects used to 
		#define the geometry, namely steps.geom.

wmgeom = swm.Geom()
		#generate parent container object


comp = swm.Comp('comp', wmgeom)
comp.addVolsys('vsys')
comp.setVol(1.6667e-27)
		#To this symple model, we only create one compartment and we 
		#store it in the variable called comp.
		
import steps.rng as srng
		#import package with alias srng.

r = srng.create('mt19937', 256)
		#we actually generate a random number generator we want.
		#using the function steps.rng.create()
		#STEPS currently only implements one pseudo RNG algorithm, 'mt19937'

r.initialize(23412)
		#initialize the random number generator with a seed value

import steps.solver as ssolver
		#to create the actual solver object, first of all we import the
		#package in which all solvers have been implemented

sim = ssolver.Wmdirect(mdl, wmgeom, r)

sim.reset()
		#reset the function on the solver project


		#we can start manipulating the "state" of simulation

sim.setCompConc('comp', 'molX1bar', 100)
sim.setCompConc('comp', 'molY1', 20)
		#this means we're setting the concentration of molX1bar to 31.4
		#micro m and molY1 to 22.3 micro m in our compartment comp.
		#we are setting these concentration value at simulation time t=0.

import numpy
		#we use Numpy to generate some auxilialy numerical arrays that will 
		#be used during simulation 

tpnt = numpy.arange(0.0, 20.001, 0.001)
		#time points
		#the range of numbers starts at 0.0 and runs to 2.0 seconds with 1ms 
		#intervals.
		#that gives us a total of 2001 "time points"

res = numpy.zeros([20001, 3])
		#res will be used to store the concentrarions of 'molX1bar' and 
		#'molY1' over time
		#that is why this array has 2001 rows and 3 columns
		#we use NumPy's zeros function, which not only allocates the array but
		#also initializes all elements to zero.

for t in range(0, 20001):
		sim.run(tpnt[t])
		res[t, 0] = sim.getCompCount('comp', 'molX1bar')
		res[t, 1] = sim.getCompCount('comp', 'molY1')
		#we loop over all time points using a range to generate indices.
		#then we use the basic solver function run to forward the simulation
		#until the time specified by the function's argument.

import pylab
		#plot values using Matplotlib 
pylab.plot(tpnt, res[:,0], label = 'X1bar')
pylab.plot(tpnt, res[:,1], label = 'Y1')

pylab.xlabel('Time(sec)')
pylab.ylabel('#molecules')
pylab.legend()
pylab.show()

		




























































