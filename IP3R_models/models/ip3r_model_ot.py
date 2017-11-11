# ip3r_model.py
#
# Katri Hituri
# IP3R model by Othmer and Tang (1993), parameter values from Diambra and Guisoni 2005
###########

# Import packages

import steps.model as smodel
import steps.geom as sgeom



# *MODEL*
#

mdl = smodel.Model()
volsys = smodel.Volsys('vsys', mdl) # Volume system
surfsys = smodel.Surfsys('ssys', mdl) # Surface system


# CHEMICaL SPECIES 
Ca = smodel.Spec('Ca', mdl)				# Calcium
IP3 = smodel.Spec('IP3', mdl)				# IP3

# IP3 receptor states
R000 = smodel.Spec('R000', mdl)     #naive state
R100 = smodel.Spec('R100', mdl)
Ropen = smodel.Spec('Ropen', mdl) #open state
R111 = smodel.Spec('R111', mdl)


# REACTIONS


#1 R000 + IP3 <=> R100
reac1_f = smodel.SReac('R000_R100', surfsys, olhs=[IP3], slhs=[R000], srhs=[R100])
reac1_b = smodel.SReac('R110_R000', surfsys, slhs=[R100], orhs=[IP3], srhs=[R000])

#2 R100 + Ca <=> Ropen
reac2_f = smodel.SReac('R100_Ropen', surfsys, olhs=[Ca], slhs=[R100], srhs=[Ropen])
reac2_b = smodel.SReac('Ropen_R100', surfsys, slhs=[Ropen], orhs=[Ca], srhs=[R100])

#3 Ropen + Ca <=> R111
reac3_f = smodel.SReac('Ropen_R111', surfsys, olhs=[Ca], slhs=[Ropen], srhs=[R111])
reac3_b = smodel.SReac('R111_Ropen', surfsys, slhs=[R111], orhs=[Ca], srhs=[Ropen])



# Rate constants and their units
reac1_f.kcst = 12e6       # 1/(Ms)  = 12 1/(mcM x s)
reac1_b.kcst = 8           # 1/s
reac2_f.kcst = 23.4e6         # 1/(Ms) = 24 1/(mcM x s)
reac2_b.kcst = 1.65          # 1/s
reac3_f.kcst = 2.81e6         # 1/(Ms)
reac3_b.kcst = 0.21          # 1/s

    

# GEOMETRY

cell = sgeom.Geom()

# Create the cytosol compartment
cyt = sgeom.Comp('cyt', cell)					
cyt.addVolsys('vsys')
cyt.setVol(0.1e-18)

# Create the endoplasmic reticulum lumen compartment
ER_lumen = sgeom.Comp('ER_lumen', cell)					
ER_lumen.addVolsys('vsys')
ER_lumen.setVol(0.1e-18)

# Create the er membrane (ER_lumen inner, cyt outer)
ER_memb = sgeom.Patch('ER_memb', cell, ER_lumen, cyt)		
ER_memb.addSurfsys('ssys')
ER_memb.setArea(0.21544369e-12)



