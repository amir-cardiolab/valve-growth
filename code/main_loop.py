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
import Pressure_trans_shift as P_trans_shift
import Pres_final_edition as Pressure_BC
from parameters import *
from continuum import *


# import Pres as Pres




#---------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-pd", "--polynomial_degree",
                    help="polynomial degree of the ansatz functions",
                    default=1,
                    type=int)
args = parser.parse_args()


set_log_level(30)
# set_log_level(16)
#---------------------------------
# Optimization options for the form compiler
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['representation'] = "uflacs"
parameters['form_compiler']['quadrature_degree'] = 4
parameters['form_compiler']['optimize'] = True
ffc_options = {'optimize' : True,
               'eliminate_zeros' : True,
               'precompute_basis_const' : True,
               'precompute_ip_const' : True}



''' Assembeling the weak form solver'''

J_b_isotropic = derivative(Functional_b_isotropic, ub, du_b)
if rank == 0:
    print ('Assembling solver...')
problem_isotropic = NonlinearVariationalProblem(Functional_b_isotropic, ub, bcs_b, J_b_isotropic)
solver_isotropic = NonlinearVariationalSolver(problem_isotropic)


###########################################################

###########################################################

'''saving everything to see they are imported correctly'''

file_BC = File(output_dir + 'facet_domains_main.pvd')
# file_BC << facet_domains
file_circum = File(output_dir + 'circum_main.pvd')
# file_circum << e_circum
file_normal_direction_loc = File(output_dir + 'normal_loc.pvd')
# file_normal_direction_loc << e_normal_loc
file_normal_direction = File(output_dir + 'normal_main.pvd')
# file_normal_direction << e_normal
file_leaflet = File(output_dir + 'leaflet_tag.pvd')
# file_leaflet << Tag_leaflet_node
Tags=Tag_leaflet_element.vector().get_local()
# file_material_3 = File(output_dir + 'tags.pvd')
file_material = File(output_dir + 'C10y.pvd')
# file_material << C10Y
file_material_2 = File(output_dir + 'C0.pvd')
# file_material_2 << C0
file_material_3= File(output_dir + 'C0_H.pvd')
file_isotropic_ub = File(output_dir + 'displacement_multi.pvd')
file_isotropic_ub_dynamic=File(output_dir + 'displacement_multi_dynamic.pvd')
file_isotropic_ub_dynamic_Global_growth=File(output_dir + 'displacement_multi_dynamic_global_growth.pvd')
file_sigma = File(output_dir + 'stress_aniso.pvd')
file_sigma1 = File(output_dir + 'stress_cycle.pvd')
Filename_ub = 'disp_anisoHumphYin_'
Filename_strain = 'strain_'

xdmffile_c0 = XDMFFile(output_dir + 'stiffness.xdmf')


########################################
#### tagging the aortic valve calcified elements to increasing local stiffening

# calcif_tags=np.zeros((element_numbers))
# calcif_tags[:]=calcif.vector().get_local()[:]
###########################################################

if rank == 0:
    print ('Entering loop:')
###########################################################
'''manual projection: to calculate strain and I_4 if needed'''

w_m = TrialFunction(V00)
v_m = TestFunction(V00)
a_m = w_m * v_m * dx
L_m = I4_iso * v_m * dx
I_4_proj = Function(V00)
# #manual projection 2:
w_m2 = TrialFunction(V00)
v_m2 = TestFunction(V00)
a_m2 = w_m2 * v_m2 * dx
L_m2 = strain_eps * v_m2 * dx
strain_eps_proj = Function(V00)

###########################################################

'''solver parameter settiing'''

#solver.parameters.update(snes_solver_parameters)
prm = solver_isotropic.parameters
#print prm["newton_solver"]["maximum_iterations"]
prm["newton_solver"]["maximum_iterations"] = 200
#prm["newton_solver"]["relaxation_parameter"] = 0.6 #.7
# prm["newton_solver"]["monitor_convergence"] = True
# prm["newton_solver"]["error_on_nonconvergence"] = False
prm["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = True
#prm["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1e-9 
prm['newton_solver']['absolute_tolerance'] = 1E-6 #1E-5
#prm['newton_solver']['relative_tolerance'] = 1E-1
prm['newton_solver']['linear_solver'] = 'tfqmr' #'gmres' bicgstab
prm['newton_solver']['preconditioner'] = 'petsc_amg'
# prm["newton_solver"]["krylov_solver"]["monitor_convergence"] = True
prm['newton_solver']['krylov_solver']['maximum_iterations'] = 20000


###########################################################
q_degree = 5
dx = dx(metadata={'quadrature_degree': q_degree})
###########################################################

###########################################################


count = 0
t_prev = 0.
count_2=0
count_3=0
N_cycles = 834  ### total number of the simulation
count=0
cycle_number=0   ### starting cycle. 

count_saving_file1=0
count_saving_file2=0

###########################################################
Flag_global_growth=True
Flag_local_growth=False
Flag_cardiac_cycle=False

###########################################################
#### linearly increasing the stiffness of the leaflet calcified elements, and the aortic root
stiffness_factor_aorta=1.003
stiffness_factore_global=1.007
stiffness_factore_local=1.01

dt_growth=0.00001  ### small time to apply nearly zero pressure when we have growth. 

C_0_temp=np.zeros((element_numbers))



###########################################################
#### main loop
###########################################################
while cycle_number <=N_cycles:
  t=0 #### initiate the time for new cycle

  vel0.vector()[:] = vel0.vector() * 0.   ### updating velocity and accelartation to the zero for new cycle.
  a0.vector()[:] = a0.vector() * 0.
  xdmffile_u = XDMFFile(output_dir + 'displacement_'+str(cycle_number)+'.xdmf')     ### saving displacement and stresses for different cycles
  xdmffile_s = XDMFFile(output_dir + 'stress_'+str(cycle_number)+'.xdmf')
  ###########################################################

  ###cardiac cycle loop
  if cycle_number%80==0 or cycle_number==834:    ###### calculating 11 cardiac cycle. 
    Flag_global_growth=False
    Flag_local_growth=False
    Flag_cardiac_cycle=True
    count_3+=1
  else:
    ###growth loop'''
    if cycle_number<165: ### till now each cycle of growth is 56.2 days. so the total cycle will be 165 total cycles - 3 cardiac cycle =162  global cycle till the 60 years old
      Flag_global_growth=True
      Flag_local_growth=False
      Flag_cardiac_cycle=False
      glob_growth_cycle=cycle_number-count3
    else:
##### after year 60 we have combination of global growth and local growth. 500 hundred of local growth cycle divided by 162 remaining will be 3 so each 3 cycle of local growth
### we calculate 1 global growth
      if (cycle_number-165)%4==0:
        Flag_global_growth=True
        Flag_local_growth=False
        Flag_cardiac_cycle=False
        glob_growth_cycle=cycle_number-3*count_2-count3### this value will be multiplied to the global growth rate (0.25). it should increase linearly. so we use this equation
### for example if the cycle_number is 169 we are in cycle of 163 of the global growth cycle(3 cycle of cardiac cycle 0,80,160) and 3 cycle of local growth cycle(166,167,168)) 
## also if we are in in 173 cycle we are in 164 global growth cycle
        count_2+=1
      else:
        Flag_global_growth=False
        Flag_local_growth=True
        Flag_cardiac_cycle=False
        local_growth_cycle=cycle_number -161 - count_2-count3 ## we have the same situation for local growth cycle. 
  ###########################################################
  if (Flag_cardiac_cycle):
        rho.assign(1100)
        while t<=0.6:
          if (1):
            if t<0.0:
              dt_s=dt_coarser
            elif t<0.04:
              dt_s=dt_fine
            elif t<.25:
              dt_s=dt_med
            else:
              dt_s=dt_coarse
          t_prev=t
          t+=dt_s
          dt.assign(dt_s)
          temp_array_P = P_trans.vector().get_local()
          temp_array_P_wall = P_wall.vector().get_local()
          temp_array_P_wall1 = P_wall1.vector().get_local()
          T_alpha = (1.0 - alpha_f) * t + (alpha_f * t_prev)
          P_current_trans = np.interp(T_alpha,Time_array,P_trans_array)
          P_current_wall = np.interp(T_alpha,Time_array,P_aort_array)
          P_current_wall1 = np.interp(T_alpha,Time_array,P_vent_array)
          shift=np.interp(T_alpha,Time_array2,P_shift)
          if (P_current_trans < 0 ):
            Side_flag.assign(0.0)
            Side_flag2.assign(1.0)
            P_current_trans_in = -P_current_trans
            vent_side=False
          else:
            P_current_trans_in=P_current_trans#+cycle_number*shift
            Side_flag.assign(1.0)
            Side_flag2.assign(0.0)
            vent_side=True
          if rank==0:
              print('cardiac_cycle')
              print ('cycle=', cycle_number)
              print ('time=', t)
              print ("tot_pressure", P_current_trans_in)
              if vent_side:
                print ("pressure_in_vent=", P_current_trans)
                print ("shifting_pressure", shift)
              else:
                print ("pressure_in aort=", P_current_trans)
          temp_array_P[:] = P_current_trans_in
          temp_array_P_wall1[:] = P_current_wall1
          temp_array_P_wall[:] = P_current_wall
          P_trans.vector().set_local(temp_array_P)
          P_wall.vector().set_local(temp_array_P_wall)
          P_wall1.vector().set_local(temp_array_P_wall1)
          solver_isotropic.solve()
          ###########################################################
          if count_saving_file1%7==0:
             ub.rename("displacement","Label")
             xdmffile_u.write(ub, t)
             stress_tensor = (1./J_e)* F_e*second_PK_stress*F_e.T
             STRESS_TENSOR_FUNCTION=TensorFunctionSpace(bmesh, 'DG', 0)
             eqSigma = Function(STRESS_TENSOR_FUNCTION)
             eqSigma.assign(project(stress_tensor,STRESS_TENSOR_FUNCTION))
             eqSigma.rename("stress","Label")
             xdmffile_s.write(eqSigma, t)
          update(ub,ub0,vel0,a0,beta2,gamma,dt_s)
          count_saving_file1+=1
  ###Growth_loop'''
  elif (Flag_global_growth):
      ###global growth loop'''
      growth_constant.assign(glob_growth_cycle*.0025)
      rho.assign(0.0)
      C_0_temp[:]=C0.vector().get_local()[:]
      for hh in range (element_numbers):
        C_0_temp[hh]=C_0_temp[hh]*stiffness_factore_global
      C0.vector().set_local(C_0_temp[:])
      while t<0.00001:
        dt_s=dt_growth
        t_prev = t
        t+= dt_s
        dt.assign(dt_s) #need to update like this!
        if rank == 0:
          # print('t',t)
          print ('cycle=', cycle_number)
          print('Global_growth_cycle')
        solver_isotropic.solve()
        update(ub,ub0,vel0,a0,beta2,gamma,dt_s)

  elif (Flag_local_growth):
    ########local growth loop######
      rho.assign(0.0)
      loc_flag.assign(1.0)
      growth_constant_loc.assign(local_growth_cycle*0.002)
      C_0_temp[:]=C0.vector().get_local()[:]
      for hh in range (element_numbers):
        if calcif_tags[hh]==1:
          C_0_temp[hh]=C_0_temp[hh]*stiffness_factore_local
      C0.vector().set_local(C_0_temp[:])
      while t<0.00001:
        dt_s=dt_growth
        t_prev = t
        t+= dt_s
        dt.assign(dt_s) #need to update like this!
        if rank == 0:
          # print('t',t)
          print ('cycle=', cycle_number)
          print('local_growth_cycle')
        solver_isotropic.solve()
        update(ub,ub0,vel0,a0,beta2,gamma,dt_s)

      
  if (cycle_number%300==0):
    xdmffile_c0.write(C0, cycle_number)
    file_material_2<<C0
  cycle_number=cycle_number+1

