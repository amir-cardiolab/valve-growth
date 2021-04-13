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


############# Some flags used for different methods.

Flag_mpi = True  ## if you want to run in parrallel or not. 
Flag_quasi=False ## if you want to run the simulation in quasi static problem, or transient. If it is false, you run the transient motion 
Flag_fenics2017=False  ### the implementation of rank (mpi) is different in Fenics 2017 and newer versions like 2018 and 2019. which version do you want to use
Flag_generalized_alpha = True  # time integration is implemented in 2 different format: NEWMARK and generalized alpha method. Generalized alpha method is more stable to converge for high accerealtion problems
Humphrey_flag = True   ## the constitutive equations are implemented in 2 different formats. Humphry and Yin will be used if Humphrey_flag=True.





if Flag_fenics2017:
	rank = MPI.rank(mpi_comm_world())
#Fenics 2018 or later
else:
	rank = MPI.rank(MPI.comm_world)