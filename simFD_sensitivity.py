#fraiman and dawson model finding sensitivity 

deterministic = False

if deterministic: NITER = 1
else: NITER = 100

#model container
import steps.model as smodel 
mdl = smodel.Model()

#spicies
Ca = smodel.Spec('Ca', mdl)
IP3 = smodel.Spec('IP3', mdl)
ER = smodel.Spec('ER', mdl)
cyt = smodel.Spec('cyt', mdl)

#receptor state objects
A00 = smodel.Spec('A00', mdl)
A10 = smodel.Spec('A10', mdl)
A11 = smodel.Spec('A11', mdl)
Oa = smodel.Spec('Oa', mdl)
Ob = smodel.Spec('Ob', mdl)
Oc = smodel.Spec('Oc', mdl)
Ia = smodel.Spec('Ia', mdl)
Ib = smodel.Spec('Ib', mdl) 
A01 = smodel.Spec('A01', mdl)
Pa = smodel.Spec('Pa', mdl)
Pb = smodel.Spec('Pb', mdl)
Pc = smodel.Spec('Pc', mdl)
Sa = smodel.Spec('Sa', mdl)
Sb = smodel.Spec('Sb', mdl)

#create the volume system
vsys = smodel.Volsys('vsys', mdl)

#surface system is similar to the valume system.
surfsys = smodel.Surfsys('ssys', mdl)


#the forward reaction without the forward binding reaction
A11_to_Oa_f = smodel.SReac('A11_to_Oa_f', surfsys, slhs = [A11], srhs = [Oa], kcst = 1800)
Oa_to_Ob_f = smodel.SReac('Oa_to_Ob_f', surfsys, slhs = [Oa], srhs = [Ob], kcst = 133)
Oc_to_Ia_f = smodel.SReac('Oc_to_Ia_f', surfsys, slhs = [Oc], srhs = [Ia], kcst = 630)
A01_to_Pa_f = smodel.SReac('A01_to_Pa_f', surfsys, slhs = [A01], srhs = [Pa], kcst = 0.3)
Pa_to_Pb_f = smodel.SReac('Pa_to_Pb_f', surfsys, slhs = [Pa], srhs = [Pb], kcst = 500)
Pc_to_Sa_f = smodel.SReac('Pc_to_Sa_f', surfsys, slhs = [Pc], srhs = [Sa], kcst = 3000)

#the backward reaction without the backrward unbinding reaction
A11_to_Oa_b = smodel.SReac('A11_to_Oa_b', surfsys, slhs = [Oa], srhs = [A11], kcst = 330)
Oa_to_Ob_b = smodel.SReac('Oa_to_Ob_b', surfsys, slhs = [Ob], srhs = [Oa], kcst = 1500) #!
Oc_to_Ia_b = smodel.SReac('Oc_to_Ia_b', surfsys, slhs = [Ia], srhs = [Oc], kcst = 400)
A01_to_Pa_b = smodel.SReac('A01_to_Pa_b', surfsys, slhs = [Pa], srhs = [A01], kcst = 700) #!
Pa_to_Pb_b = smodel.SReac('Pa_to_Pb_b', surfsys, slhs = [Pb], srhs = [Pa], kcst = 100)
Pc_to_Sa_b = smodel.SReac('Pc_to_Sa_b', surfsys, slhs = [Sa], srhs = [Pc], kcst = 250)



#the forward bindng reactions
A00_bind_IP3_f = smodel.SReac('A00_bind_IP3_f', surfsys, olhs = [IP3], slhs = [A00], srhs = [A10])

A10_bind_Ca_f = smodel.SReac('A10_bind_Ca_f', surfsys, olhs = [Ca], slhs = [A10], srhs = [A11])
Ob_bind_Ca_f = smodel.SReac('Ob_bind_Ca_f', surfsys, olhs = [Ca], slhs = [Ob], srhs = [Oc])
Ia_bind_Ca_f = smodel.SReac('Ia_bind_Ca_f', surfsys, olhs = [Ca], slhs = [Ia], srhs = [Ib])

A00_bind_Ca_f = smodel.SReac('A00_bind_Ca_f', surfsys, olhs = [Ca], slhs = [A00], srhs = [A01])
Pb_bind_Ca_f = smodel.SReac('Pb_bind_Ca_f', surfsys, olhs = [Ca], slhs = [Pb], srhs = [Pc])
Sa_bind_Ca_f = smodel.SReac('Sa_bind_Ca_f', surfsys, olhs = [Ca], slhs = [Sa], srhs = [Sb])

A01_bind_IP3_f = smodel.SReac('A01_bind_IP3_f', surfsys, olhs = [IP3], slhs = [A01], srhs = [A11])

#the backward unbinding reactions
A00_bind_IP3_b = smodel.SReac('A00_bind_IP3_b', surfsys, slhs = [A10], orhs = [IP3], srhs = [A00])

A10_bind_Ca_b = smodel.SReac('A10_bind_Ca_b', surfsys, slhs = [A11], orhs = [Ca], srhs = [A10])
Ob_bind_Ca_b = smodel.SReac('Ob_bind_Ca_b', surfsys, slhs = [Oc], orhs = [Ca], srhs = [Ob])
Ia_bind_Ca_b = smodel.SReac('Ia_bind_Ca_b', surfsys, slhs = [Ib], orhs = [Ca], srhs = [Ia])

A00_bind_Ca_b = smodel.SReac('A00_bind_Ca_b', surfsys, slhs = [A01], orhs = [Ca], srhs = [A00])
Pb_bind_Ca_b = smodel.SReac('Pb_bind_Ca_b', surfsys, slhs = [Pc], orhs = [Ca], srhs = [Pb])
Sa_bind_Ca_b = smodel.SReac('Sa_bind_Ca_b', surfsys, slhs = [Sb], orhs = [Ca], srhs = [Sa])

A01_bind_IP3_b = smodel.SReac('A01_bind_IP3_b', surfsys, slhs = [A11], orhs = [IP3], srhs = [A01])

#Ca ion passing through open IP3R channel
A00_Oa_Ca_channel_f = smodel.SReac('A00_0a_Ca_channel_f', surfsys, ilhs = [Ca], slhs = [Oa], orhs = [Ca], srhs = [Oa], kcst = 2e8)
A00_Ob_Ca_channel_f = smodel.SReac('A00_Ob_channel_f', surfsys, ilhs = [Ca], slhs = [Ob], orhs = [Ca], srhs = [Ob], kcst = 2e8)
A00_Oc_Ca_channel_f = smodel.SReac('A00_Oc_Ca_channel_f', surfsys, ilhs = [Ca], slhs = [Oc], orhs = [Ca], srhs = [Oc], kcst = 2e8)



#set constant
A00_bind_IP3_f.setKcst(6670e6)

A10_bind_Ca_f.setKcst(500e6)		
Ob_bind_Ca_f.setKcst(70e6)		
Ia_bind_Ca_f.setKcst(60e6)		

A00_bind_Ca_f.setKcst(5000e6)			
Pb_bind_Ca_f.setKcst(5000e6)			
Sa_bind_Ca_f.setKcst(5000e6)			

A01_bind_IP3_f.setKcst(1540e6)



A00_bind_IP3_b.setKcst(200)	

A10_bind_Ca_b.setKcst(667)				 
Ob_bind_Ca_b.setKcst(2000)					
Ia_bind_Ca_b.setKcst(16)				

A00_bind_Ca_b.setKcst(1)				
Pb_bind_Ca_b.setKcst(150)				
Sa_bind_Ca_b.setKcst(20)				

A01_bind_IP3_b.setKcst(18)

#geometry setup
import steps.geom as swm
wmgeom = swm.Geom()

#create the Cytosole compartment
cyt = swm.Comp('cyt',wmgeom)
cyt.addVolsys('vsys')
cyt.setVol(0.1e-18)

#Create the Emdoplasmic Reticulum compartment
ER = swm.Comp('ER', wmgeom)
ER.addVolsys('vsys')
ER.setVol(0.1e-18)

#ER it the inner compertment, cyr is the outer compartment
memb = swm.Patch('memb', wmgeom, ER, cyt)
memb.addSurfsys('ssys')
memb.setArea(0.21544368e-12)

print 'Inner compartment to memb is', memb.getIComp().getID()
print 'Outer compartment to memb is', memb.getOComp().getID()

#RNG setup
import steps.rng as srng
import pylab
import numpy

#finished creating the components
r = srng.create('mt19937', 1000)
r.initialize(7233)
#solver setup
import steps.solver as ssolver
if deterministic: sim = ssolover.Wmrk4(mdl, wmgeom, r)
else: sim = ssolver.Wmdirect(mdl, wmgeom, r)

# Number of iterations (defines how many times the model is simulated)
NITER = 20
		#you can change here

# Ca2+ concentrations in cytosol
ca_concs = numpy.arange(0.00, 1.01, 0.01)  
#time
tpnt = numpy.arange(0.0, 1.01, 0.01)

# array for simulation results
res = numpy.zeros([ca_concs.size, tpnt.size])

for i in xrange(ca_concs.size):

	print 'Round', i+1, '/', ca_concs.size
	temp_res = numpy.zeros([NITER, tpnt.size]) # temporary storage for results

	for j in xrange(NITER): 
		sim.reset()
        	sim.setPatchCount('memb', 'A00', 1) # number of naive receptor
        	sim.setCompConc('cyt', 'IP3', 10e-6) # [IP3] = 10 uM
        	sim.setCompClamped('cyt', 'IP3', 1)
        	sim.setCompConc('cyt', 'Ca', ca_concs[i])
        	sim.setCompClamped('cyt', 'Ca', 1)
        	sim.setCompConc('ER', 'Ca', 150e-6) 
        	sim.setCompClamped('ER', 'Ca', 1)

    
        for t in xrange(tpnt.size):
            		sim.run(tpnt[t])
            		o1 = sim.getPatchCount('memb', 'Oa')
            		o2 = sim.getPatchCount('memb', 'Ob')
            		o3 = sim.getPatchCount('memb', 'Oc')  
            		temp_res[j,t] = o1 + o2 + o3 

            		#print 'o1 is', sim.getPatchCount('memb', 'Oa')
            		#print 'o2 is', sim.getPatchCount('memb', 'Ob')
            		#print 'o3 is', sim.getPatchCount('memb', 'Oc')
                    


	# calculate the mean and standard deviation of the simulation results   
	#temp = numpy.mean(temp_res[:,2001:]) # take only into account results after 20 s
	#res[i,0] = numpy.mean(temp, 0)
    #res[i,1] = numpy.std(temp, 0)



# save the results (means and stds)
numpy.savetxt('ip3r_fd_time_results.dat', res)

print #res[:,0]

