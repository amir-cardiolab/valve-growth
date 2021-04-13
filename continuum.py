from dolfin import *
import numpy as np
from numpy import linalg as LA
import time
import sys
import argparse
import glob
import re
import math
from Flags import *

from functions import *
from parameters import *



#Generalized alpha method
# 0.5 #0 annihilates high freq in one step. 1: high freq preserved


######################################
# anisotropic problem
######################################

############# total deformation gradient implementation


I = Identity(3)   # define second order identity
if (Flag_generalized_alpha ): #generalized alpha integrator
 F_t = I + grad( (1.0 - alpha_f)* ub + alpha_f * ub0)    ##### total deformation gradient
else: #Newmark or typical integrators
 F_t = I + grad(ub)

J_t=det(F_t) #jacobian
################### Growth problem implementation #################

loc_flag=Constant(0.0)  ## will turn on after 60 years old

glob_flag=Constant(1.0)

n_dir_tensor=outer(e_normal,e_normal)  ##### making global normal direction tensor 

n_dir_tensor_loc=outer(e_normal_loc,e_normal_loc)   ##### making local normal direction tensor 

growth_constant=Constant(0.0) ## it is zero for healthy case. it will be increased during the main loop

growth_constant_loc=Constant(0.0)  ## it is zero for healthy case. it will be increased during the main loop

F_growth=growth_constant*n_dir_tensor+growth_constant_loc*n_dir_tensor_loc*loc_flag+I  #### total growth part of deformation gradient. includes global and local growth

F_growth_inv=inv(F_growth)


#################### elastic part of deformation gradient

F_e = F_t*F_growth_inv

C_e = dot(F_e.T,F_e)  ### right cauchy deformation tensor

J_e = det(F_e)

E = 0.5*(C_e - I)  # Cauchy strain tensor

F_iso = J_e**(-1.0/3.0)*F_e #isotropic  part of elastic deformation gradient 

C_iso = dot(F_iso.T, F_iso)   

# I1 = tr(C_e)  # first invariant of C
# I2 =  0.5*( tr(C_e)**2 - tr(C_e*C_e) )  # second invariant of C


I1_iso = tr(C_iso)
# I2_iso =  0.5*( tr(C_iso)**2 - tr(C_iso*C_iso) )


########## circumferential direction and I_4 implementation

circ_Dir_tensor = outer(e_circum ,e_circum)

I4_iso = inner(C_iso, circ_Dir_tensor) 


##### calculating strain if neeeded
E_output  = 0.5*(C_e - I)

strain_eps = inner(E_output, circ_Dir_tensor)

############## nearly incompressible material implementation 

Kappa = Function(Vb_scalar)

temp_array_Kappa = Kappa.vector().get_local() 

temp_array_Kappa[:] = Kappa_incomp

Kappa.vector().set_local(temp_array_Kappa)

dU_dJ = Kappa*(J_e-1)


####### definition of S_volume
S_vol = J_e*dU_dJ*inv(C_e)

############### derivative of constitutive equation dpsi/dc #############


 ######## Isotropic terms included wall and leaflets. C0 defined as 0 for wall and 1 for the leaflet. so the second equation will turn zero for the wall

d_psi_bar_d_C_bar = ( C10Y + C0 * C_H * C1_H * exp(C1_H*(I1_iso - 3.0)) ) * I   

######## ansotropic term 
d_psi_bar_d_C_bar = d_psi_bar_d_C_bar + ( C0 * C0_H * C2_H * (1.0 - (1./(I4_iso**.5)) ) * exp(C2_H* ( (I4_iso**.5) - 1.0)**2)  ) *circ_Dir_tensor

####### definition of S_isochoric
S_isc = 2.0*J_e**(-2.0/3.0)*Dev(d_psi_bar_d_C_bar, C_e)

second_PK_stress =  S_vol + S_isc


#### Inititation of velocity and acceleration based on the generalized alpha method

vel = (gamma/(beta2*dt))*(ub - ub0) - (gamma/beta2 - 1.0)*vel0 - dt*(gamma/(2.0*beta2) - 1.0)*a0 #vel and acceleration at t_(n+1)
a  = (1.0/(beta2*dt**2))*(ub - ub0 - dt*vel0) - (1.0/(2.0*beta2) - 1.0) * a0
if (Flag_generalized_alpha ): #generalized alpha integrator
 a = (1.0 - alpha_m) * a + alpha_m * a0
 vel = (1.0 - alpha_f) * vel + alpha_f * vel0


##########################################################
############ Implementation of neuman boundary conditions and forming a weak form

P_trans = Function(Vb_scalar)
P_wall = Function(Vb_scalar)
P_wall1 = Function(Vb_scalar)

########### Sideflag is a flag to see if we want to apply the transvalvular pressure to the ventricular side of the valve or aortic side. 
########### it turns on when the transvalvular pressure is positive. and we apply the transvalvular pressure to the ventricular side of the leaflet
Side_flag = Constant(1.0)
Functional_b_isotropic = inner(dot(F_e,second_PK_stress),grad(vb))*dx() + Side_flag* P_trans*J_e*inner(inv(F_e.T)*n_function,vb)*ds(1)

########### Side_flag2 turns on when the transvalvular pressure is negative. and we apply the transvalvular pressure to the aortic side of the leaflet
Side_flag2 = Constant(0.0)
Functional_b_isotropic = Functional_b_isotropic + Side_flag2 * P_trans*J_e* inner(inv(F_e.T)*n_function,vb)*(ds(2))

########### applying pressure to the aortic root and ascending aorta
Functional_b_isotropic = Functional_b_isotropic +  P_wall*J_e*inner(inv(F_e.T)*n_function,vb)* (ds(4)) 
Functional_b_isotropic = Functional_b_isotropic +  P_wall1*J_e*inner(inv(F_e.T)*n_function,vb)* (ds(3)) 

#### acceleration term applied when we have dynamic problem
Functional_b_isotropic = Functional_b_isotropic + rho * dot(a,vb) * dx() #+damp*dot(vel,vb)*dx()

########## applying symmetry conditions to the symmetry walls
Functional_b_isotropic = Functional_b_isotropic + beta/h * dot(ub,n_function) * dot(vb,n_function)* ds(5) + beta/h * dot(ub,n_function) * dot(vb,n_function)* ds(8)








