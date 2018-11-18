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

\begin{lstlisting}
	region ID style args keyword arg ...
\end{lstlisting}

\begin{itemize}
\item ID = user-assigned name for the region
\item style = delete or block or cone or cylinder or plane or prism or sphere or union or intersect
\item zero or more keyword/arg pairs may be appended
\item keyword = side or units or move or rotate or open
\item accelerated styles (with same args) = block/kk
\end{itemize}


This modification add a keyword voronoi which will be followed by three arguments nx ny nz.

\begin{lstlisting}
voronoi keywords args = nx ny nz
	nx = the integer number of cells in the x direction of this region
	ny = the integer number of cells in the y direction of this region
	nz = the integer number of cells in the z direction of this region
\end{lstlisting}
