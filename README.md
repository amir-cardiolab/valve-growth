# valve-growth
1-Introduction:

This package is an open-source nonlinear transient structural mechanics solver written in Python (FEniCS). The package couples growth and remodeling (G&R) of a biological tissue/structure with the transient dynamic of the structure. Specifically, we developed a solver to couple the transient dynamic of an aortic valve with the growth and remodeling of the aortic valve due to aging and calcification for multiscale modeling of disease growth. 
The solver runs with MPI and interfaces, through FEniCS, to the state-of-the-art linear algebra backends like PETSc. To be able to work with this package, a user should know basic knowledge in how to solve PDEs through FEniCS, and basic programming skills in Python. A good introduction to the FEniCS could be obtained by following the tutorial: 
https://fenicsproject.org/tutorial/

This code is used in the following publication:



2- FEniCS installation:
To run the solver, users need to install FEniCS, an open-source finite element software, written in Python: the best reference to install FEniCS is: 
https://fenicsproject.org/download/
	
This solver is compatible with FEniCS 2017, 2018 and 2019. However, it is recommended to install this version of FEniCS. 

$ conda create -n fenicsproject -c conda-forge fenics=2019.2
$ source activate fenicsproject

3- Nonlinear transient solver installation:
a.	The solver can be installed simply by cloning the GitHub repository to your own computer:


	$ git clone https://github.com/amir-cardiolab/valve-growth.git
	$ cd code
	
	
b.	You can just download the package from the link below:


	$ https://github.com/amir-cardiolab/valve-growth.git
	$ unzip code.zip
	$ cd code

4- Files and folders:
	The solver has a folder named code, which includes six python codes: 
All of the codes are implemented in Python and the main executable code name is main_loop.py. 
To run the solver, you need to preset the assumptions:

a.	Flags.py  includes some of the settings used in the solver e.g., which version of FEniCS you need to use. By turning on or off the flags the settings may be changed. All of the flags have comments in the code and are explained in the code. 
b.	functions.py   includes some of the functions, which are used in the solver. 
c.	Parameters.py  includes all of the parameters used in the solver. The time integration, the material properties, the circumferential and normal directions for growth and remodeling are implemented in this file. You can change the material properties such as the constitutive equation, the constitutive equation constants, time integration method, and parameters. 
d.	continuum.py   includes the G&R equations and Cauchy’s equation of motion. Also, the weak form of Cauchy’s equation is assembled in this code. 
e.	main_loop.py  is the executable code, which is included the main loop where you can couple the transient structural dynamics with G&R. 
To run the solver you need to run the main_loop.py:
$ python main_loop.py

 To download the geometry  you can use the google drive link: https://drive.google.com/file/d/1QVIGgRxDCMKsZHVPIoSGhBz3DxKBMSlv/view?usp=sharing
 the geometry should be in the same folder as the codes are.
 
If you want to run in parallel (recommended), you need to run the code as:

	$ time mpi -n n_processors python main_loop.py

In the above “n_processors” is the number of processors you want to use in parallel (MPI).


5- Post processing: 

You can see the results in the open-source software Paraview. To download the software, you can use the link below: 
https://www.paraview.org/download/
	
To see the displacement and the growth of the aortic valve, you need to use 
“warp by vector” option in Filters->alphabetical->warp by vector:

<img width="468" src="https://user-images.githubusercontent.com/61292399/114954222-ea25ab00-9e0e-11eb-941c-eea6ef1bda13.png">



Figure-2 is an example of the visualization of displacement of the valve which is  the output of the solver:

<img width="468" alt="Picture1" src="https://user-images.githubusercontent.com/61292399/114954318-1e996700-9e0f-11eb-85b4-e18dc97c2d41.png">






