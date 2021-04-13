
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
import Pres_final_edition as Pressure_BC
import Pressure_trans_shift as P_trans_shift


###########################################################

''' Loading Geometry (contains directions, material stiffness, and boundaries) '''

mesh_filename='./whol_model_n.h5'
output_dir='./results/'
bmesh = Mesh()
if (rank == 0):
    print ('Loading mesh:', mesh_filename)
hdf = HDF5File(bmesh.mpi_comm(), mesh_filename, "r")
hdf.read(bmesh, "mesh",False)


####### circumferential direction for the fibers in constitutive equation (I_4)  
#### for numerical issues, it is better to import circumferential direction cellwise (error raises in linear solver with nodewise circumferential direction)
e_circum=Function(VectorFunctionSpace(bmesh,'DG',0))
hdf.read(e_circum,'e_circum')

####### Normal direction for global thicknening and growth 
#### ill condition happened when we used cellwise normal direction. it reaches maximum iteration in linear solver
e_normal = Function(VectorFunctionSpace(bmesh,'CG',1))
hdf.read(e_normal, "e_normal")


####### Normal direction for local thicknening and growth (nodewise)
e_normal_loc = Function(VectorFunctionSpace(bmesh,'CG',1))
hdf.read(e_normal_loc, "e_normal_loc_2la")


####### tagging leaflet to use two different strain energy in a single equation 
Tag_leaflet_element = Function(FunctionSpace(bmesh, 'DG', 0))
hdf.read(Tag_leaflet_element,'Tag_leaflet')
element_numbers = bmesh.num_cells()
Tags_element =  np.zeros((element_numbers, 1),dtype="intc")
for gol in range(element_numbers):
  Tags_element[gol] = Tag_leaflet_element.vector()[gol]


####### All surfaces of the geometry for applying boundary conditions 
facet_domains = MeshFunction('size_t', bmesh,2, bmesh.domains())
hdf.read(facet_domains,'facet_domains')
ds = Measure("ds")[facet_domains]


####### calcified elements tags for increasing the local stiffening of these elements
V0=FunctionSpace(bmesh,'CG',1)
V00=FunctionSpace(bmesh,'DG',0)
Vb = VectorFunctionSpace(bmesh, 'Lagrange', 2)
Vb_scalar = FunctionSpace(bmesh, 'DG', 0)



###########################################################
#### ###### Time integration (generalized alpha method) [J Chung et al 1993]
###########################################################



rho_inf = 0.0 
alpha_m =  (2.0* rho_inf - 1.0 ) / (1.0 + rho_inf)
alpha_f = rho_inf / ( 1.0 + rho_inf)
if (Flag_generalized_alpha):
  gamma = 0.5 - alpha_m + alpha_f
  beta2 = 0.25*( ( 1.0 - alpha_m + alpha_f )**2 )
  if rank == 0:
    print ('Using Generalized-alpha!!!!! ')


#### total time of simulation (T) and the timestep sizes. using different timestep size to reduce the computational costs
T = .6 #0.84   #integration time
t = 0.0
dt = Constant(5e-5) #2e-4
dt_fine = 5e-5    #5e-5 fails at the start of the steep acceleration starting at the end of 1st cycle
dt_finer = 1e-6
dt_coarse = 1e-3
dt_med = 5e-4

###########################################################
#### material properties
###########################################################
node_numbers = bmesh.num_vertices()
element_numbers = bmesh.num_cells()
connectivity_matrix = bmesh.cells()

rho = Constant(1100.0)  #material density
C10Y_aorta = 2e6/ 6.  
C0_aorta = 0.0 #Aorta is neo-Hookean 

###Humphrey-alpha model (Auricchio, Ferrara, Morganti Ann. Solid Struct Mech 2012)

C1_H = 5.81  #Auricchio et al 2012 paper
C2_H = 24.97  #Auricchio et al 2012 paper
C0_H_s = 0.062e6 #Auricchio et al 2012 paper
C_H_s = 0.022e6

C_H = Function(V00)
C0_H = Function(V00)
C_H_temp = C_H.vector().get_local()
C_H_temp[:] = C_H_s
C_H.vector().set_local(C_H_temp)
C_H.vector().apply('insert')
C0_H_temp = C0_H.vector().get_local()
C0_H_temp[:] = C0_H_s
C0_H.vector().set_local(C0_H_temp)
C0_H.vector().apply('insert')

C10Y_values = [0 , C10Y_aorta ]
C0_values = [1.0,  C0_aorta ]


'''changing the material property to numpy array'''
C10Y = Function(V00)
C0=Function(V00)
local_range = V0.dofmap().ownership_range()
local_dim = local_range[1] - local_range[0]
normal_numpy_flatten=e_normal.vector().get_local()
normal_numpy=np.reshape(normal_numpy_flatten,(local_dim,3))

sss = element_numbers
C10Y_numpy = np.zeros(sss)
C0_numpy = np.zeros(sss)
damp_values=[4e5,5e9]
damp_regions_numpy=np.zeros(sss)
for cell_no in range(sss):
    if (Tags_element[cell_no] == 1 ):
     subdomain_no = 0
    else:
     subdomain_no = 1
    C10Y_numpy[cell_no] = C10Y_values[subdomain_no]
    C0_numpy[cell_no] = C0_values[subdomain_no]
    damp_regions_numpy[cell_no]=damp_values[subdomain_no]

C10Y.vector().set_local(C10Y_numpy)
C0.vector().set_local(C0_numpy)
C10Y.vector().apply('insert')
C0.vector().apply('insert')


### Penalty term for symmetry boundary condition
h = CellDiameter(bmesh)
beta = 1e7 #1e7  #penalty for Nitsche
n_function = FacetNormal(bmesh)

##### penalty term for the nearly incompressible materials
Kappa_incomp = 1e6 #5e8 #5e6





###########################################################
#### Adding damp coeffiecient
###########################################################

damp_values=[4e5,5e9]
damp_regions_numpy=np.zeros(sss)

damp_co=Function(V00)
damp_co.vector().set_local(damp_regions_numpy)
damp_co.vector().apply('insert')



if rank == 0:
 print ('number of elements:', element_numbers)
 print ('number of nodes:', node_numbers)



###########################################################
####  boundary conditions
###########################################################

'''Drichlet Boundary conditions'''
bc1 = DirichletBC(Vb.sub(1), 0.0 , facet_domains, 6) #wall_bottom  (they cannot move in axial direction)
bc2 = DirichletBC(Vb.sub(1), 0.0 , facet_domains, 7) #wall_top


bc3z_top = DirichletBC(Vb,(0.0,0.0,0.0), Pinpoint([.0108, .043, 0]), 'pointwise')   ##### four point rules. 2 points of each end wall mentioned above, is fixed inn any direction
bc4x_top = DirichletBC(Vb, (0.0,0.0,0.0), Pinpoint([0.0054, 0.043, -0.00935307]), 'pointwise')
bc3z_bot = DirichletBC(Vb, (0.0,0.0,0.0), Pinpoint([0.00599995, -.015, -0.0103922]), 'pointwise')
bc4x_bot = DirichletBC(Vb, (0.0,0.0,0.0), Pinpoint([0.0120001, -0.015, 6.899e-10]), 'pointwise')
bcs_b = [bc1,bc2,bc3z_top,bc4x_top,bc3z_bot,bc4x_bot]


'''neuman Boundary conditions'''
Time_array = Pressure_BC.time_BC()
Time_array2=P_trans_shift.time_BC()
P_shift=P_trans_shift.P_trans_added()
P_trans_array = Pressure_BC.P_trans_BC()
P_aort_array = Pressure_BC.P_aort_BC()
P_vent_array = Pressure_BC.P_vent_BC()
F_interp_trans = np.interp(t,Time_array,P_trans_array)
F_interp_aort = np.interp(t,Time_array,P_aort_array)
F_interp_vent = np.interp(t,Time_array,P_vent_array)

###########################################################


'''Trial and test functions'''

du_b = TrialFunction(Vb)   # displacement
vb = TestFunction(Vb)
ub = Function(Vb) # the most recently computed solution
vel0, a0 , ub0 = Function(Vb), Function(Vb), Function(Vb)

