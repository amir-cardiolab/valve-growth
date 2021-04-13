

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

###########################################################

''' Loading Geometry (contains direction, material stiffness, and boundaries) '''


mesh_filename='./whole_model.h5'
bmesh = Mesh()
if (rank == 0):
    print ('Loading mesh:', mesh_filename)
hdf = HDF5File(bmesh.mpi_comm(), mesh_filename, "r")
hdf.read(bmesh, "mesh",False)





e_circum=Function(VectorFunctionSpace(bmesh,'DG',0))
e_normal = Function(VectorFunctionSpace(bmesh,'CG',1))
e_normal_loc = Function(VectorFunctionSpace(bmesh,'CG',1))
Tag_leaflet_element = Function(FunctionSpace(bmesh, 'DG', 0))
facet_domains = MeshFunction('size_t', bmesh,2, bmesh.domains())
hdf.read(e_normal, "e_normal")
hdf.read(e_normal_loc, "e_normal_loc_2la")
hdf.read(facet_domains,'facet_domains')
hdf.read(e_circum,'e_circum')
hdf.read(Tag_leaflet_element,'Tag_leaflet')
calcif = Function(FunctionSpace(bmesh, 'DG', 0))
hdf.read(calcif,'calcif')
file_calcif = File(output_dir + 'calcif.pvd')
file_calcif << calcif
ds = Measure("ds")[facet_domains]


V0=FunctionSpace(bmesh,'CG',1)
V00=FunctionSpace(bmesh,'DG',0)
Vb = VectorFunctionSpace(bmesh, 'Lagrange', 2)
Vb_scalar = FunctionSpace(bmesh, 'DG', 0)

element_numbers = bmesh.num_cells()
Tags_element =  np.zeros((element_numbers, 1),dtype="intc")
for gol in range(element_numbers):
  Tags_element[gol] = Tag_leaflet_element.vector()[gol]



