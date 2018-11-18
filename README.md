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

