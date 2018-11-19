# LAMMPS
# Author: Jianfeng Huang
# SEBE, Glasgow Caledonian University, Unitied Kingdom, G4 0BA

Content:

1 Voronoi tessellation.

The Voronoi polycrystalline structure integrated in LAMMPS.

File name:
region.h
region.cpp

The LAMMPS official region commond is defined as in LAMMPS site.

	region ID style args keyword arg ...

 ID = user-assigned name for the region
 
 style = delete or block or cone or cylinder or plane or prism or sphere or union or intersect
    zero or more keyword/arg pairs may be appended
    
 keyword = side or units or move or rotate or open
    accelerated styles (with same args) = block/kk

This modification add a keyword voronoi which will be followed by three arguments nx ny nz.

voronoi keywords args = nx ny nz

	nx = the integer number of cells in the x direction of this region
	
	ny = the integer number of cells in the y direction of this region
	
	nz = the integer number of cells in the z direction of this region
	

2 Dislocation recognition

The dislocation class is an object oriented model as a result of the dislocation recognition algorithm. This class will be instanced once a dislocation in MD model is recognised and will be updated if the dislocation has any changes in the current step.

Dislocation file:
dislocation.h
dislocation.cpp

3 Dislocation pair force

The pair of dislocation is used to calculate the dislocation potential based on the dislocation recognition method.

There are two files: 

pair_dislocation.h
pair_dislocation.cpp

4 Compute command

The compute command in LAMMPS has been revised with a new style 'dislocation'. The usage of this command in LAMMPS is:

compute ID group-ID style args

    ID = user-assigned name for the computation
    group-ID = ID of the group of atoms to perform the computation on
    style = one of a list of possible style names (see below)
    args = arguments used by a particular style


The style is 'dislocation' added with a args to determine the cutoff radius in this computing as a parameter.

For example:

compute 1 all dislocation 3

files:
compute_dislocation.h
compute_dislocation.cpp
