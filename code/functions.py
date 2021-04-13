
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

''' Quasi-state or unsteady Velocity Update'''
if Flag_quasi:

	def update(ub,ub0):
	  ub_vec, ub0_vec = ub.vector() , ub0.vector()
	  ub0.vector()[:] = ub.vector()
else:
	def update(ub,ub0,vel0,a0,beta2,gamma,dt):
	  ub_vec, ub0_vec = ub.vector(), ub0.vector()
	  vel0_vec, a0_vec = vel0.vector(), a0.vector()
	    #update acceleration and vel
	  a_vec = (1.0/(2.0*beta2))*(  (ub_vec -ub0_vec - vel0_vec*dt) / (0.5*dt*dt) - (1.0 -2.0*beta2)*a0_vec)
	  v_vec = dt*((1.0-gamma) * a0_vec + gamma*a_vec ) + vel0_vec
	  #update
	  vel0.vector()[:] , a0.vector()[:] = v_vec , a_vec
	  ub0.vector()[:] = ub.vector()




def Dev(mytensor, C_e):
    return mytensor - 1.0/3.0*(inner(mytensor, C_e))*inv(C_e)
mesh_edge_size = 5e-4
TOL = mesh_edge_size/3.0
class Pinpoint(SubDomain):
    
    def __init__(self, coords):
        self.coords = np.array(coords)
        SubDomain.__init__(self)
    def move(self, coords):
        self.coords[:] = np.array(coords)
    def inside(self, x, on_boundary):
        return np.linalg.norm(x-self.coords) < TOL





