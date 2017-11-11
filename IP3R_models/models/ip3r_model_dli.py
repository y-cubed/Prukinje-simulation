# ip3r_model.py
#
# Katri Hituri
# IP3R model by Dawson, Lea, Irvine (2003)
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
RP = smodel.Spec('RP', mdl)

R = smodel.Spec('R', mdl) #naive state 2
RI = smodel.Spec('RI', mdl) 
RI2 = smodel.Spec('RI2', mdl)
RI3 = smodel.Spec('RI3', mdl)

P = smodel.Spec('P', mdl)     #naive state 1
PI = smodel.Spec('PI', mdl)
PI2 = smodel.Spec('PI2', mdl)
PI3 = smodel.Spec('PI3', mdl)

C1 = smodel.Spec('C1', mdl)
O1 = smodel.Spec('O1', mdl) # open state 1
O2 = smodel.Spec('O2', mdl) # open state 1


# REACTIONS


#1 R <=> P
reac1_f = smodel.SReac('R_P', surfsys, slhs=[R], srhs=[P])
reac1_b = smodel.SReac('P_R', surfsys, slhs=[P], srhs=[R])

#2 R + IP3 <=> RI
reac2_f = smodel.SReac('R_RI', surfsys, olhs=[IP3], slhs=[R], srhs=[RI])
reac2_b = smodel.SReac('RI_R', surfsys, slhs=[RI], orhs=[IP3], srhs=[R])

#3 RI + IP3 <=> RI2
reac3_f = smodel.SReac('RI_RI2', surfsys, olhs=[IP3], slhs=[RI], srhs=[RI2])
reac3_b = smodel.SReac('RI2_RI', surfsys, slhs=[RI2], orhs=[IP3], srhs=[RI])

#4 RI2 + IP3 <=> RI3
reac4_f = smodel.SReac('RI2_RI3', surfsys, olhs=[IP3], slhs=[RI2], srhs=[RI3])
reac4_b = smodel.SReac('RI3_RI2', surfsys, slhs=[RI3], orhs=[IP3], srhs=[RI2])

#5 RI3 + IP3 <=> O1
reac5_f = smodel.SReac('RI3_O1', surfsys, olhs=[IP3], slhs=[RI3], srhs=[O1])
reac5_b = smodel.SReac('O1_RI3', surfsys, slhs=[O1], orhs=[IP3], srhs=[RI3])

#6 P + IP3 <=> PI
reac6_f = smodel.SReac('P_PI', surfsys, olhs=[IP3], slhs=[P], srhs=[PI])
reac6_b = smodel.SReac('PI_P', surfsys, slhs=[PI], orhs=[IP3], srhs=[P])

#7 PI + IP3 <=> PI2
reac7_f = smodel.SReac('PI_PI2', surfsys, olhs=[IP3], slhs=[PI], srhs=[PI2])
reac7_b = smodel.SReac('PI2_PI', surfsys, slhs=[PI2], orhs=[IP3], srhs=[PI])

#8 PI2 + IP3 <=> PI3
reac8_f = smodel.SReac('PI2_PI3', surfsys, olhs=[IP3], slhs=[PI3], srhs=[PI3])
reac8_b = smodel.SReac('PI3_PI2', surfsys, slhs=[PI3], orhs=[IP3], srhs=[PI2])

#9 PI3 + IP3 <=> C1
reac9_f = smodel.SReac('PI3_C1', surfsys, olhs=[IP3], slhs=[PI3], srhs=[C1])
reac9_b = smodel.SReac('C1_PI3', surfsys, slhs=[C1], orhs=[IP3], srhs=[PI3])

#10 RI <=> PI
reac10_f = smodel.SReac('RI_PI', surfsys, slhs=[RI], srhs=[PI])
reac10_b = smodel.SReac('PI_RI', surfsys, slhs=[PI], srhs=[RI])

#11 RI2 <=> PI2
reac11_f = smodel.SReac('RI2_PI2', surfsys, slhs=[RI2], srhs=[PI2])
reac11_b = smodel.SReac('PI2_RI2', surfsys, slhs=[PI2], srhs=[RI2])

#12 RI3 <=> PI3
reac12_f = smodel.SReac('RI3_PI3', surfsys, slhs=[RI3], srhs=[PI3])
reac12_b = smodel.SReac('PI3_RI3', surfsys, slhs=[PI3], srhs=[RI3])

#13 O1 <=> C1
reac13_f = smodel.SReac('O1_C1', surfsys, slhs=[O1], srhs=[C1])
reac13_b = smodel.SReac('C1_O1', surfsys, slhs=[C1], srhs=[O1])

#14 (flux thru O1) not included

#15 O1 + Ca <=> O2
reac15_f = smodel.SReac('O1_O2', surfsys, olhs=[Ca], slhs=[O1], srhs=[O2])
reac15_b = smodel.SReac('O2_O1', surfsys, slhs=[O2], orhs=[Ca], srhs=[O1])

#16 (flux thru O2) and 17 (diffucion of Ca away from channel mouth) not included

#18 R + Ca <=> RP
reac18_f = smodel.SReac('R_RP', surfsys, olhs=[Ca], slhs=[R], srhs=[RP])
reac18_b = smodel.SReac('RP_R', surfsys, slhs=[RP], orhs=[Ca], srhs=[R])

#19 P + Ca <=> RP
reac19_f = smodel.SReac('P_RP', surfsys, olhs=[Ca], slhs=[P], srhs=[RP])
reac19_b = smodel.SReac('RP_P', surfsys, slhs=[RP], orhs=[Ca], srhs=[P])



# Rate constants and their units
reac1_f.kcst = 1           # 1/s
reac1_b.kcst = 100         # 1/s

reac2_f.kcst = 4000e6      # 1/(Ms) 
reac2_b.kcst = 1000        # 1/s

reac3_f.kcst = 3000e6      # 1/(Ms)
reac3_b.kcst = 2000        # 1/s

reac4_f.kcst = 2000e6      # 1/Ms
reac4_b.kcst = 3000        # 1/s

reac5_f.kcst = 1000e6      # 1/Ms
reac5_b.kcst = 4000        # 1/s

reac6_f.kcst = 400e6       # 1/Ms
reac6_b.kcst = 10          # 1/s

reac7_f.kcst = 300e6       # 1/Ms
reac7_b.kcst = 20          # 1/s

reac8_f.kcst = 200e6       # 1/Ms
reac8_b.kcst = 30          # 1/s

reac9_f.kcst = 100e6       # 1/Ms
reac9_b.kcst = 40          # 1/s

reac10_f.kcst = 1          # 1/s
reac10_b.kcst = 10         # 1/s

reac11_f.kcst = 1          # 1/s
reac11_b.kcst = 1          # 1/s

reac12_f.kcst = 10         # 1/s
reac12_b.kcst = 1          # 1/s

reac13_f.kcst = 10         # 1/s
reac13_b.kcst = 0.1        # 1/s

reac15_f.kcst = 100e6      # 1/Ms
reac15_b.kcst = 10         # 1/s

reac18_f.kcst = 1e6        # 1/Ms
reac18_b.kcst = 0.1        # 1/s

reac19_f.kcst = 10e6       # 1/Ms
reac19_b.kcst = 0.01       # 1/s

    

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



