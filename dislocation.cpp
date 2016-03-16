/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   * Jianfeng Huang  20-Nov-2015  Change for dislocation simulation in MD
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "dislocation.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "comm.h"
#include "domain.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "lattice.h"
#include <algorithm>

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 256
#define EPSILON 1.0e-7
#define BIG 1.0e20

#define SINERTIA 0.4            // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

Dislocation::Dislocation(LAMMPS *lmp, int numAtoms, tagint* atomIdsIn) : Pointers(lmp)
{
  me = comm->me;

  // initialize all fields to empty

  initialize();

  // scan file for sizes of all fields and allocate them
  natoms = numAtoms;
  atomIds = new std::vector<tagint>;
  newAtomIds = new std::vector<tagint>;
  for (int i = 0; i < natoms; i++) {
	  atomIds->push_back(atomIdsIn[i]); newAtomIds->push_back(atomIdsIn[i]);
  }

  allocate();

  // read file again to populate all fields

  // stats
  char tempId[255];
  static int number = 1;
  sprintf(tempId, "%d", number++);
  int n = strlen(tempId) + 1;
  id = new char[n];
  sprintf(id, "Dislocation[%s]", tempId);

 // if (me == 0) {
 //   if (screen)
 //     fprintf(screen,"Create dislocation %s:\n"
 //             "  %d atoms with %d types\n  %d bonds with %d types\n"
 //             "  %d angles with %d types\n  %d dihedrals with %d types\n"
 //             "  %d impropers with %d types\n",
 //             id,natoms,ntypes,
 //             nbonds,nbondtypes,nangles,nangletypes,
 //             ndihedrals,ndihedraltypes,nimpropers,nimpropertypes);
	//if (logfile)
	//{
	//	fprintf(logfile, "Create dislocation %s:\n"
	//		"  %d atoms with %d types\n  %d bonds with %d types\n"
	//		"  %d angles with %d types\n  %d dihedrals with %d types\n"
	//		"  %d impropers with %d types\n",
	//		id, natoms, ntypes,
	//		nbonds, nbondtypes, nangles, nangletypes,
	//		ndihedrals, ndihedraltypes, nimpropers, nimpropertypes);
	//	fprintf(logfile, "atoms list:\n");
 //       std::vector<tagint>::iterator it;
 //       for (it=atomIds->begin(); it!=atomIds->end(); it++)
	//	{
	//		fprintf(logfile, " %d ", *it);
	//	}
	//	fprintf(logfile, "\n");
	//}
 // }
  CalcSlipSystem();
}

/* ---------------------------------------------------------------------- */
Dislocation::~Dislocation()
{
  //delete [] id;
  deallocate();
}

/* ----------------------------------------------------------------------
   compute center = geometric center of dislocation
   also compute:
     dx = displacement of each atom from center
     molradius = radius of dislocation from center including finite-size particles
------------------------------------------------------------------------- */

void Dislocation::compute_center()
{
    double ** x = atom->x;
    std::vector<tagint> catoms;
    std::vector<tagint>::iterator disAtomIt;
    double deltXij[3];
    int n = 1;
    int nerror = 0;
    std::vector<tagint> temp;
    center[0] = center[1] = center[2] = 0.0;

    {
        center[0] += x[*disAtomIt][0];
        center[1] += x[*disAtomIt][1];
        center[2] += x[*disAtomIt][2];
    }

    center[0] /= natoms;
    center[1] /= natoms;
    center[2] /= natoms;

    memory->destroy(dx);
    memory->create(dx, natoms, 3, "dislocation:dx");
    n = 0;
    for (disAtomIt = atomIds->begin(); disAtomIt != atomIds->end(); disAtomIt++)
    {
        dx[n][0] = x[*disAtomIt][0] - center[0];
        dx[n][1] = x[*disAtomIt][1] - center[1];
        dx[n][2] = x[*disAtomIt][2] - center[2];
        n++;
    }

    molradius = 0.0;
    for (int i = 0; i < natoms; i++) {
        double rad = MathExtra::len3(dx[i]);
        if (radiusflag) rad += radius[i];
        molradius = MAX(molradius, rad);
    }
}

/* ----------------------------------------------------------------------
   compute masstotal = total mass of dislocation
   could have been set by user, otherwise calculate it
------------------------------------------------------------------------- */

void Dislocation::compute_mass()
{
  if (massflag) return;
  massflag = 1;

  if (!rmassflag) atom->check_mass();

  masstotal = 0.0;
  for (int i = 0; i < natoms; i++) {
    if (rmassflag) masstotal += rmass[i];
    else masstotal += atom->mass[type[i]];
  }
}

/* ----------------------------------------------------------------------
   compute com = center of mass of dislocation
   could have been set by user, otherwise calculate it
   works for finite size particles assuming no overlap
   also compute:
     dxcom = displacement of each atom from COM
     comatom = which atom (1-Natom) is nearest the COM
     maxextent = furthest any atom in dislocation is from comatom (not COM)
------------------------------------------------------------------------- */

void Dislocation::compute_com()
{
  if (!comflag) {
    comflag = 1;

    if (!rmassflag) atom->check_mass();

    double onemass;
    com[0] = com[1] = com[2] = 0.0;
    for (int i = 0; i < natoms; i++) {
      if (rmassflag) onemass = rmass[i];
      else onemass = atom->mass[type[i]];
      com[0] += x[i][0]*onemass;
      com[1] += x[i][1]*onemass;
      com[2] += x[i][2]*onemass;
    }
    com[0] /= masstotal;
    com[1] /= masstotal;
    com[2] /= masstotal;
  }

  memory->destroy(dxcom);
  memory->create(dxcom,natoms,3,"dislocation:dxcom");

  for (int i = 0; i < natoms; i++) {
    dxcom[i][0] = x[i][0] - com[0];
    dxcom[i][1] = x[i][1] - com[1];
    dxcom[i][2] = x[i][2] - com[2];
  }

  double rsqmin = BIG;
  for (int i = 0; i < natoms; i++) {
    double rsq = MathExtra::lensq3(dxcom[i]);
    if (rsq < rsqmin) {
      comatom = i;
      rsqmin = rsq;
    }
  }

  double rsqmax = 0.0;
  for (int i = 0; i < natoms; i++) {
    double dx = x[comatom][0] - x[i][0];
    double dy = x[comatom][1] - x[i][1];
    double dz = x[comatom][2] - x[i][2];
    double rsq = dx*dx + dy*dy + dz*dz;
    rsqmax = MAX(rsqmax,rsq);
  }

  comatom++;
  maxextent = sqrt(rsqmax);
}

/* ----------------------------------------------------------------------
   compute itensor = 6 moments of inertia of dislocation around xyz axes
   could have been set by user, otherwise calculate it
   accounts for finite size spheres, assuming no overlap
   also compute:
     inertia = 3 principal components of inertia
     ex,ey,ez = principal axes in space coords
     quat = quaternion for orientation of dislocation
     dxbody = displacement of each atom from COM in body frame
------------------------------------------------------------------------- */

void Dislocation::compute_inertia()
{
  if (!inertiaflag) {
    inertiaflag = 1;

    if (!rmassflag) atom->check_mass();

    double onemass,dx,dy,dz;
    for (int i = 0; i < 6; i++) itensor[i] = 0.0;
    for (int i = 0; i < natoms; i++) {
      if (rmassflag) onemass = rmass[i];
      else onemass = atom->type[type[i]];
      dx = dxcom[i][0];
      dy = dxcom[i][1];
      dz = dxcom[i][2];
      itensor[0] += onemass * (dy*dy + dz*dz);
      itensor[1] += onemass * (dx*dx + dz*dz);
      itensor[2] += onemass * (dx*dx + dy*dy);
      itensor[3] -= onemass * dy*dz;
      itensor[4] -= onemass * dx*dz;
      itensor[5] -= onemass * dx*dy;
    }

    if (radiusflag) {
      for (int i = 0; i < natoms; i++) {
        if (rmassflag) onemass = rmass[i];
        else onemass = atom->type[type[i]];
        itensor[0] += SINERTIA*onemass * radius[i]*radius[i];
        itensor[1] += SINERTIA*onemass * radius[i]*radius[i];
        itensor[2] += SINERTIA*onemass * radius[i]*radius[i];
      }
    }
  }

  // diagonalize inertia tensor for each body via Jacobi rotations
  // inertia = 3 eigenvalues = principal moments of inertia
  // evectors and exzy = 3 evectors = principal axes of rigid body

  double cross[3];
  double tensor[3][3],evectors[3][3];

  tensor[0][0] = itensor[0];
  tensor[1][1] = itensor[1];
  tensor[2][2] = itensor[2];
  tensor[1][2] = tensor[2][1] = itensor[3];
  tensor[0][2] = tensor[2][0] = itensor[4];
  tensor[0][1] = tensor[1][0] = itensor[5];
  
  if (MathExtra::jacobi(tensor,inertia,evectors))
    error->all(FLERR,"Insufficient Jacobi rotations for rigid dislocation");
  
  ex[0] = evectors[0][0];
  ex[1] = evectors[1][0];
  ex[2] = evectors[2][0];
  ey[0] = evectors[0][1];
  ey[1] = evectors[1][1];
  ey[2] = evectors[2][1];
  ez[0] = evectors[0][2];
  ez[1] = evectors[1][2];
  ez[2] = evectors[2][2];

  // if any principal moment < scaled EPSILON, set to 0.0
  
  double max;
  max = MAX(inertia[0],inertia[1]);
  max = MAX(max,inertia[2]);
  
  if (inertia[0] < EPSILON*max) inertia[0] = 0.0;
  if (inertia[1] < EPSILON*max) inertia[1] = 0.0;
  if (inertia[2] < EPSILON*max) inertia[2] = 0.0;
  
  // enforce 3 evectors as a right-handed coordinate system
  // flip 3rd vector if needed
  
  MathExtra::cross3(ex,ey,cross);
  if (MathExtra::dot3(cross,ez) < 0.0) MathExtra::negate3(ez);
  
  // create quaternion
  
  MathExtra::exyz_to_q(ex,ey,ez,quat);

  // compute displacements in body frame defined by quat

  memory->destroy(dxbody);
  memory->create(dxbody,natoms,3,"dislocation:dxbody");

  for (int i = 0; i < natoms; i++)
    MathExtra::transpose_matvec(ex,ey,ez,dxcom[i],dxbody[i]);
}

/* ----------------------------------------------------------------------
   read dislocation info from file
   flag = 0, just scan for sizes of fields
   flag = 1, read and store fields
------------------------------------------------------------------------- */

void Dislocation::read(int flag)
{
  char line[MAXLINE],keyword[MAXLINE];
  char *eof,*ptr;

  // skip 1st line of file

  if (me == 0) {
    eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of dislocation file");
  }

  // read header lines
  // skip blank lines or lines that start with "#"
  // stop when read an unrecognized line

  while (1) {

    readline(line);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // search line for header keywords and set corresponding variable

    if (strstr(line,"atoms")) sscanf(line,"%d",&natoms);
    else if (strstr(line,"bonds")) sscanf(line,"%d",&nbonds);
    else if (strstr(line,"angles")) sscanf(line,"%d",&nangles);
    else if (strstr(line,"dihedrals")) sscanf(line,"%d",&ndihedrals);
    else if (strstr(line,"impropers")) sscanf(line,"%d",&nimpropers);

    else if (strstr(line,"mass")) {
      massflag = 1;
      sscanf(line,"%lg",&masstotal);
    }
    else if (strstr(line,"com")) {
      comflag = 1;
      sscanf(line,"%lg %lg %lg",&com[0],&com[1],&com[2]);
      if (domain->dimension == 2 && com[2] != 0.0)
        error->all(FLERR,"Dislocation file z center-of-mass must be 0.0 for 2d");
    }
    else if (strstr(line,"inertia")) {
      inertiaflag = 1;
      sscanf(line,"%lg %lg %lg %lg %lg %lg",
             &itensor[0],&itensor[1],&itensor[2],
             &itensor[3],&itensor[4],&itensor[5]);
    }

    else break;
  }

  // error check

  if (flag == 0) {
    if (natoms == 0) error->all(FLERR,"No atom count in dislocation file");
  }

  // count = vector for tallying bonds,angles,etc per atom

  if (flag == 0) memory->create(count,natoms,"dislocation:count");
  else count = NULL;
  
  // grab keyword and skip next line

  parse_keyword(0,line,keyword);
  readline(line);

  // loop over sections of dislocation file

  while (strlen(keyword)) {
    if (strcmp(keyword,"Coords") == 0) {
      xflag = 1;
      if (flag) coords(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Types") == 0) {
      typeflag = 1;
      if (flag) types(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Charges") == 0) {
      qflag = 1;
      if (flag) charges(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Diameters") == 0) {
      radiusflag = 1;
      if (flag) diameters(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Masses") == 0) {
      rmassflag = 1;
      if (flag) masses(line);
      else skip_lines(natoms,line);

    } else if (strcmp(keyword,"Bonds") == 0) {
      if (nbonds == 0)
	error->all(FLERR,"Dislocation file has bonds but no nbonds setting");
      bondflag = tag_require = 1;
      bonds(flag,line);
    } else if (strcmp(keyword,"Angles") == 0) {
      if (nangles == 0) 
	error->all(FLERR,"Dislocation file has angles but no nangles setting");
      angleflag = tag_require = 1;
      angles(flag,line);
    } else if (strcmp(keyword,"Dihedrals") == 0) {
      if (ndihedrals == 0) error->all(FLERR,"Dislocation file has dihedrals "
				      "but no ndihedrals setting");
      dihedralflag = tag_require = 1;
      dihedrals(flag,line);
    } else if (strcmp(keyword,"Impropers") == 0) {
      if (nimpropers == 0) error->all(FLERR,"Dislocation file has impropers "
				      "but no nimpropers setting");
      improperflag = tag_require = 1;
      impropers(flag,line);

    } else if (strcmp(keyword,"Special Bond Counts") == 0) {
      nspecialflag = 1;
      nspecial_read(flag,line);
    } else if (strcmp(keyword,"Special Bonds") == 0) {
      specialflag = tag_require = 1;
      if (flag) special_read(line);
      else skip_lines(natoms,line);

    } else if (strcmp(keyword,"Shake Flags") == 0) {
      shakeflagflag = 1;
      if (flag) shakeflag_read(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Shake Atoms") == 0) {
      shakeatomflag = tag_require = 1;
      if (shaketypeflag) shakeflag = 1;
      if (!shakeflagflag) 
	error->all(FLERR,"Dislocation file shake flags not before shake atoms");
      if (flag) shakeatom_read(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Shake Bond Types") == 0) {
      shaketypeflag = 1;
      if (shakeatomflag) shakeflag = 1;
      if (!shakeflagflag) 
	error->all(FLERR,"Dislocation file shake flags not before shake bonds");
      if (flag) shaketype_read(line);
      else skip_lines(natoms,line);

    } else error->one(FLERR,"Unknown section in dislocation file");
	     
    parse_keyword(1,line,keyword);
  }

  // clean up

  memory->destroy(count);

  // error check

  if (flag == 0) {
    if ((nspecialflag && !specialflag) || (!nspecialflag && specialflag))
      error->all(FLERR,"Dislocation file needs both Special Bond sections");
    if (specialflag && !bondflag) 
      error->all(FLERR,"Dislocation file has special flags but no bonds");
    if (!specialflag && bondflag) 
      error->all(FLERR,"Dislocation file has bonds but no special flags");

    if ((shakeflagflag || shakeatomflag || shaketypeflag) && !shakeflag)
      error->all(FLERR,"Dislocation file shake info is incomplete");
  }
}

/* ----------------------------------------------------------------------
   read coords from file
------------------------------------------------------------------------- */

void Dislocation::coords(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %lg %lg %lg",&tmp,&x[i][0],&x[i][1],&x[i][2]);
  }

  if (domain->dimension == 2) {
    for (int i = 0; i < natoms; i++)
      if (x[i][2] != 0.0) 
        error->all(FLERR,"Dislocation file z coord must be 0.0 for 2d");
  }
}

/* ----------------------------------------------------------------------
   read types from file
   set ntypes = max of any atom type
------------------------------------------------------------------------- */

void Dislocation::types(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %d",&tmp,&type[i]);
  }

  for (int i = 0; i < natoms; i++)
    if (type[i] <= 0)
      error->all(FLERR,"Invalid atom type in dislocation file");

  for (int i = 0; i < natoms; i++)
    ntypes = MAX(ntypes,type[i]);
}

/* ----------------------------------------------------------------------
   read charges from file
------------------------------------------------------------------------- */

void Dislocation::charges(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %lg",&tmp,&q[i]);
  }
}

/* ----------------------------------------------------------------------
   read diameters from file and set radii
------------------------------------------------------------------------- */

void Dislocation::diameters(char *line)
{
  int tmp;
  maxradius = 0.0;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %lg",&tmp,&radius[i]);
    radius[i] *= 0.5;
    maxradius = MAX(maxradius,radius[i]);
  }

  for (int i = 0; i < natoms; i++)
    if (radius[i] < 0.0) 
      error->all(FLERR,"Invalid atom diameter in dislocation file");
}

/* ----------------------------------------------------------------------
   read masses from file
------------------------------------------------------------------------- */

void Dislocation::masses(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %lg",&tmp,&rmass[i]);
  }

  for (int i = 0; i < natoms; i++)
    if (rmass[i] <= 0.0) error->all(FLERR,"Invalid atom mass in dislocation file");
}

/* ----------------------------------------------------------------------
   read bonds from file
   set nbondtypes = max type of any bond
   store each with both atoms if newton_bond = 0
   if flag = 0, just count bonds/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Dislocation::bonds(int flag, char *line)
{
  int tmp,itype;
  tagint m,atom1,atom2;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_bond[i] = 0;

  for (int i = 0; i < nbonds; i++) {
    readline(line);
    sscanf(line,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT,
           &tmp,&itype,&atom1,&atom2);

    if (atom1 <= 0 || atom1 > natoms ||
	atom2 <= 0 || atom2 > natoms)
      error->one(FLERR,"Invalid atom ID in Bonds section of dislocation file");
    if (itype <= 0)
      error->one(FLERR,"Invalid bond type in Bonds section of dislocation file");

    if (flag) {
      m = atom1-1;
      nbondtypes = MAX(nbondtypes,itype);
      bond_type[m][num_bond[m]] = itype;
      bond_atom[m][num_bond[m]] = atom2;
      num_bond[m]++;
      if (newton_bond == 0) {
	m = atom2-1;
	bond_type[m][num_bond[m]] = itype;
	bond_atom[m][num_bond[m]] = atom1;
	num_bond[m]++;
      }
    } else {
      count[atom1-1]++;
      if (newton_bond == 0) count[atom2-1]++;
    }
  }

  // bond_per_atom = max of count vector

  if (flag == 0) {
    bond_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      bond_per_atom = MAX(bond_per_atom,count[i]);
  }
}

/* ----------------------------------------------------------------------
   read angles from file
   store each with all 3 atoms if newton_bond = 0
   if flag = 0, just count angles/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Dislocation::angles(int flag, char *line)
{
  int tmp,itype;
  tagint m,atom1,atom2,atom3;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_angle[i] = 0;

  for (int i = 0; i < nangles; i++) {
    readline(line);
    sscanf(line,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT,
           &tmp,&itype,&atom1,&atom2,&atom3);

    if (atom1 <= 0 || atom1 > natoms ||
        atom2 <= 0 || atom2 > natoms ||
        atom3 <= 0 || atom3 > natoms)
      error->one(FLERR,"Invalid atom ID in Angles section of dislocation file");
    if (itype <= 0)
      error->one(FLERR,"Invalid angle type in Angles section of dislocation file");

    if (flag) {
      m = atom2-1;
      nangletypes = MAX(nangletypes,itype);
      angle_type[m][num_angle[m]] = itype;
      angle_atom1[m][num_angle[m]] = atom1;
      angle_atom2[m][num_angle[m]] = atom2;
      angle_atom3[m][num_angle[m]] = atom3;
      num_angle[m]++;
      if (newton_bond == 0) {
	m = atom1-1;
	angle_type[m][num_angle[m]] = itype;
	angle_atom1[m][num_angle[m]] = atom1;
	angle_atom2[m][num_angle[m]] = atom2;
	angle_atom3[m][num_angle[m]] = atom3;
	num_angle[m]++;
	m = atom3-1;
	angle_type[m][num_angle[m]] = itype;
	angle_atom1[m][num_angle[m]] = atom1;
	angle_atom2[m][num_angle[m]] = atom2;
	angle_atom3[m][num_angle[m]] = atom3;
	num_angle[m]++;
      }
    } else {
      count[atom2-1]++;
      if (newton_bond == 0) {
	count[atom1-1]++;
	count[atom3-1]++;
      }
    }
  }

  // angle_per_atom = max of count vector

  if (flag == 0) {
    angle_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      angle_per_atom = MAX(angle_per_atom,count[i]);
  }
}

/* ----------------------------------------------------------------------
   read dihedrals from file
   store each with all 4 atoms if newton_bond = 0
   if flag = 0, just count dihedrals/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Dislocation::dihedrals(int flag, char *line)
{
  int tmp,itype;
  tagint m,atom1,atom2,atom3,atom4;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_dihedral[i] = 0;

  for (int i = 0; i < ndihedrals; i++) {
    readline(line);
    sscanf(line,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT " " 
           TAGINT_FORMAT " " TAGINT_FORMAT " ",
           &tmp,&itype,&atom1,&atom2,&atom3,&atom4);

    if (atom1 <= 0 || atom1 > natoms ||
        atom2 <= 0 || atom2 > natoms ||
        atom3 <= 0 || atom3 > natoms ||
        atom4 <= 0 || atom4 > natoms)
      error->one(FLERR,
		 "Invalid atom ID in dihedrals section of dislocation file");
    if (itype <= 0)
      error->one(FLERR,
		 "Invalid dihedral type in dihedrals section of dislocation file");

    if (flag) {
      m = atom2-1;
      ndihedraltypes = MAX(ndihedraltypes,itype);
      dihedral_type[m][num_dihedral[m]] = itype;
      dihedral_atom1[m][num_dihedral[m]] = atom1;
      dihedral_atom2[m][num_dihedral[m]] = atom2;
      dihedral_atom3[m][num_dihedral[m]] = atom3;
      dihedral_atom4[m][num_dihedral[m]] = atom4;
      num_dihedral[m]++;
      if (newton_bond == 0) {
	m = atom1-1;
	dihedral_type[m][num_dihedral[m]] = itype;
	dihedral_atom1[m][num_dihedral[m]] = atom1;
	dihedral_atom2[m][num_dihedral[m]] = atom2;
	dihedral_atom3[m][num_dihedral[m]] = atom3;
	dihedral_atom4[m][num_dihedral[m]] = atom4;
	num_dihedral[m]++;
	m = atom3-1;
	dihedral_type[m][num_dihedral[m]] = itype;
	dihedral_atom1[m][num_dihedral[m]] = atom1;
	dihedral_atom2[m][num_dihedral[m]] = atom2;
	dihedral_atom3[m][num_dihedral[m]] = atom3;
	dihedral_atom4[m][num_dihedral[m]] = atom4;
	num_dihedral[m]++;
	m = atom4-1;
	dihedral_type[m][num_dihedral[m]] = itype;
	dihedral_atom1[m][num_dihedral[m]] = atom1;
	dihedral_atom2[m][num_dihedral[m]] = atom2;
	dihedral_atom3[m][num_dihedral[m]] = atom3;
	dihedral_atom4[m][num_dihedral[m]] = atom4;
	num_dihedral[m]++;
      }
    } else {
      count[atom2-1]++;
      if (newton_bond == 0) {
	count[atom1-1]++;
	count[atom3-1]++;
	count[atom4-1]++;
      }
    }
  }

  // dihedral_per_atom = max of count vector

  if (flag == 0) {
    dihedral_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      dihedral_per_atom = MAX(dihedral_per_atom,count[i]);
  }
}

/* ----------------------------------------------------------------------
   read impropers from file
   store each with all 4 atoms if newton_bond = 0
   if flag = 0, just count impropers/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Dislocation::impropers(int flag, char *line)
{
  int tmp,itype;
  tagint m,atom1,atom2,atom3,atom4;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_improper[i] = 0;

  for (int i = 0; i < nimpropers; i++) {
    readline(line);
    sscanf(line,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT " " 
           TAGINT_FORMAT " " TAGINT_FORMAT " ",
           &tmp,&itype,&atom1,&atom2,&atom3,&atom4);

    if (atom1 <= 0 || atom1 > natoms ||
        atom2 <= 0 || atom2 > natoms ||
        atom3 <= 0 || atom3 > natoms ||
        atom4 <= 0 || atom4 > natoms)
      error->one(FLERR,
		 "Invalid atom ID in impropers section of dislocation file");
    if (itype <= 0)
      error->one(FLERR,
		 "Invalid improper type in impropers section of dislocation file");

    if (flag) {
      m = atom2-1;
      nimpropertypes = MAX(nimpropertypes,itype);
      improper_type[m][num_improper[m]] = itype;
      improper_atom1[m][num_improper[m]] = atom1;
      improper_atom2[m][num_improper[m]] = atom2;
      improper_atom3[m][num_improper[m]] = atom3;
      improper_atom4[m][num_improper[m]] = atom4;
      num_improper[m]++;
      if (newton_bond == 0) {
	m = atom1-1;
	improper_type[m][num_improper[m]] = itype;
	improper_atom1[m][num_improper[m]] = atom1;
	improper_atom2[m][num_improper[m]] = atom2;
	improper_atom3[m][num_improper[m]] = atom3;
	improper_atom4[m][num_improper[m]] = atom4;
	num_improper[m]++;
	m = atom3-1;
	improper_type[m][num_improper[m]] = itype;
	improper_atom1[m][num_improper[m]] = atom1;
	improper_atom2[m][num_improper[m]] = atom2;
	improper_atom3[m][num_improper[m]] = atom3;
	improper_atom4[m][num_improper[m]] = atom4;
	num_improper[m]++;
	m = atom4-1;
	improper_type[m][num_improper[m]] = itype;
	improper_atom1[m][num_improper[m]] = atom1;
	improper_atom2[m][num_improper[m]] = atom2;
	improper_atom3[m][num_improper[m]] = atom3;
	improper_atom4[m][num_improper[m]] = atom4;
	num_improper[m]++;
      }
    } else {
      count[atom2-1]++;
      if (newton_bond == 0) {
	count[atom1-1]++;
	count[atom3-1]++;
	count[atom4-1]++;
      }
    }
  }

  // improper_per_atom = max of count vector

  if (flag == 0) {
    improper_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      improper_per_atom = MAX(improper_per_atom,count[i]);
  }
}

/* ----------------------------------------------------------------------
   read 3 special bonds counts from file
   if flag = 0, just tally maxspecial
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Dislocation::nspecial_read(int flag, char *line)
{
  int tmp,c1,c2,c3;

  if (flag == 0) maxspecial = 0;

  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %d %d %d",&tmp,&c1,&c2,&c3);

    if (flag) {
      nspecial[i][0] = c1;
      nspecial[i][1] = c1+c2;
      nspecial[i][2] = c1+c2+c3;
    } else maxspecial = MAX(maxspecial,c1+c2+c3);
  }
}

/* ----------------------------------------------------------------------
   read special bond indices from file
------------------------------------------------------------------------- */

void Dislocation::special_read(char *line)
{
  int m,nwords;
  char **words = new char*[maxspecial+1];

  for (int i = 0; i < natoms; i++) {
    readline(line);
    nwords = parse(line,words,maxspecial+1);
    if (nwords != nspecial[i][2]+1)
      error->all(FLERR,"Dislocation file special list "
		 "does not match special count");

    for (m = 1; m < nwords; m++) {
      special[i][m-1] = ATOTAGINT(words[m]);
      if (special[i][m-1] <= 0 || special[i][m-1] > natoms ||
	  special[i][m-1] == i+1)
	error->all(FLERR,"Invalid special atom index in dislocation file");
    }
  }

  delete [] words;
}

/* ----------------------------------------------------------------------
   read SHAKE flags from file
------------------------------------------------------------------------- */

void Dislocation::shakeflag_read(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %d",&tmp,&shake_flag[i]);
  }

  for (int i = 0; i < natoms; i++)
    if (shake_flag[i] < 0 || shake_flag[i] > 4) 
      error->all(FLERR,"Invalid shake flag in dislocation file");
}

/* ----------------------------------------------------------------------
   read SHAKE atom info from file
------------------------------------------------------------------------- */

void Dislocation::shakeatom_read(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (shake_flag[i] == 1)
      sscanf(line,"%d " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT,
             &tmp,&shake_atom[i][0],&shake_atom[i][1],&shake_atom[i][2]);
    else if (shake_flag[i] == 2)
      sscanf(line,"%d " TAGINT_FORMAT " " TAGINT_FORMAT,
             &tmp,&shake_atom[i][0],&shake_atom[i][1]);
    else if (shake_flag[i] == 3)
      sscanf(line,"%d " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT,
             &tmp,&shake_atom[i][0],&shake_atom[i][1],&shake_atom[i][2]);
    else if (shake_flag[i] == 4)
      sscanf(line,"%d " TAGINT_FORMAT " " TAGINT_FORMAT " " 
             TAGINT_FORMAT " " TAGINT_FORMAT,
             &tmp,&shake_atom[i][0],&shake_atom[i][1],
             &shake_atom[i][2],&shake_atom[i][3]);
  }

  for (int i = 0; i < natoms; i++) {
    int m = shake_flag[i];
    if (m == 1) m = 3;
    for (int j = 0; j < m; j++)
      if (shake_atom[i][j] <= 0 || shake_atom[i][j] > natoms)
        error->all(FLERR,"Invalid shake atom in dislocation file");
  }
}

/* ----------------------------------------------------------------------
   read SHAKE bond type info from file
------------------------------------------------------------------------- */

void Dislocation::shaketype_read(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (shake_flag[i] == 1)
      sscanf(line,"%d %d %d %d",&tmp,
             &shake_type[i][0],&shake_type[i][1],&shake_type[i][2]);
    else if (shake_flag[i] == 2)
      sscanf(line,"%d %d",&tmp,&shake_type[i][0]);
    else if (shake_flag[i] == 3)
      sscanf(line,"%d %d %d",&tmp,&shake_type[i][0],&shake_type[i][1]);
    else if (shake_flag[i] == 4)
      sscanf(line,"%d %d %d %d",&tmp,
             &shake_type[i][0],&shake_type[i][1],&shake_type[i][2]);
  }

  for (int i = 0; i < natoms; i++) {
    int m = shake_flag[i];
    if (m == 1) m = 3;
    for (int j = 0; j < m-1; j++)
      if (shake_type[i][j] <= 0)
        error->all(FLERR,"Invalid shake bond type in dislocation file");
    if (shake_flag[i] == 1)
      if (shake_type[i][2] <= 0)
        error->all(FLERR,"Invalid shake angle type in dislocation file");
  }
}

/* ----------------------------------------------------------------------
   init all data structures to empty
------------------------------------------------------------------------- */

void Dislocation::initialize()
{
  id = NULL;
  natoms = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;
  ntypes = 0;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;
  dislocationLength = 0;

  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
  maxspecial = 0;

  xflag = typeflag = qflag = radiusflag = rmassflag = 0;
  bondflag = angleflag = dihedralflag = improperflag = 0;
  nspecialflag = specialflag = 0;
  shakeflag = shakeflagflag = shakeatomflag = shaketypeflag = 0;

  centerflag = massflag = comflag = inertiaflag = 0;
  tag_require = 0;

  x = NULL;
  type = NULL;
  q = NULL;
  radius = NULL;
  rmass = NULL;

  num_bond = NULL;
  bond_type = NULL;
  bond_atom = NULL;

  num_angle = NULL;
  angle_type = NULL;
  angle_atom1 = angle_atom2 = angle_atom3 = NULL;

  num_dihedral = NULL;
  dihedral_type = NULL;
  dihedral_atom1 = dihedral_atom2 = dihedral_atom3 = dihedral_atom4 = NULL;

  num_improper = NULL;
  improper_type = NULL;
  improper_atom1 = improper_atom2 = improper_atom3 = improper_atom4 = NULL;

  nspecial = NULL;
  special = NULL;

  shake_flag = NULL;
  shake_atom = NULL;
  shake_type = NULL;

  dx = NULL;
  dxcom = NULL;
  dxbody = NULL;
  burgersMap = new std::map<tagint, SlipVector*>;
}

/* ----------------------------------------------------------------------
   allocate all data structures
   also initialize values for data structures that are always allocated
------------------------------------------------------------------------- */

void Dislocation::allocate()
{
  if (xflag) memory->create(x,natoms,3,"dislocation:x");
  if (typeflag) memory->create(type,natoms,"dislocation:type");
  if (qflag) memory->create(q,natoms,"dislocation:q");
  if (radiusflag) memory->create(radius,natoms,"dislocation:radius");
  if (rmassflag) memory->create(rmass,natoms,"dislocation:rmass");

  // always allocate num_bond,num_angle,etc and special+nspecial
  // even if not in dislocation file, initialize to 0
  // this is so methods that use these arrays don't have to check they exist

  memory->create(num_bond,natoms,"dislocation:num_bond");
  for (int i = 0; i < natoms; i++) num_bond[i] = 0;
  memory->create(num_angle,natoms,"dislocation:num_angle");
  for (int i = 0; i < natoms; i++) num_angle[i] = 0;
  memory->create(num_dihedral,natoms,"dislocation:num_dihedral");
  for (int i = 0; i < natoms; i++) num_dihedral[i] = 0;
  memory->create(num_improper,natoms,"dislocation:num_improper");
  for (int i = 0; i < natoms; i++) num_improper[i] = 0;

  memory->create(special,natoms,maxspecial,"dislocation:special");
  memory->create(nspecial,natoms,3,"dislocation:nspecial");
  for (int i = 0; i < natoms; i++)
    nspecial[i][0] = nspecial[i][1] = nspecial[i][2] = 0;

  if (bondflag) {
    memory->create(bond_type,natoms,bond_per_atom,
		   "dislocation:bond_type");
    memory->create(bond_atom,natoms,bond_per_atom,
		   "dislocation:bond_atom");
  }

  if (angleflag) {
    memory->create(angle_type,natoms,angle_per_atom,
		   "dislocation:angle_type");
    memory->create(angle_atom1,natoms,angle_per_atom,
		   "dislocation:angle_atom1");
    memory->create(angle_atom2,natoms,angle_per_atom,
		   "dislocation:angle_atom2");
    memory->create(angle_atom3,natoms,angle_per_atom,
		   "dislocation:angle_atom3");
  }

  if (dihedralflag) {
    memory->create(dihedral_type,natoms,dihedral_per_atom,
		   "dislocation:dihedral_type");
    memory->create(dihedral_atom1,natoms,dihedral_per_atom,
		   "dislocation:dihedral_atom1");
    memory->create(dihedral_atom2,natoms,dihedral_per_atom,
		   "dislocation:dihedral_atom2");
    memory->create(dihedral_atom3,natoms,dihedral_per_atom,
		   "dislocation:dihedral_atom3");
    memory->create(dihedral_atom4,natoms,dihedral_per_atom,
		   "dislocation:dihedral_atom4");
  }

  if (improperflag) {
    memory->create(improper_type,natoms,improper_per_atom,
		   "dislocation:improper_type");
    memory->create(improper_atom1,natoms,improper_per_atom,
		   "dislocation:improper_atom1");
    memory->create(improper_atom2,natoms,improper_per_atom,
		   "dislocation:improper_atom2");
    memory->create(improper_atom3,natoms,improper_per_atom,
		   "dislocation:improper_atom3");
    memory->create(improper_atom4,natoms,improper_per_atom,
		   "dislocation:improper_atom4");
  }

  if (shakeflag) {
    memory->create(shake_flag,natoms,"dislocation:shake_flag");
    memory->create(shake_atom,natoms,4,"dislocation:shake_flag");
    memory->create(shake_type,natoms,3,"dislocation:shake_flag");
  }
}

/* ----------------------------------------------------------------------
   deallocate all data structures
------------------------------------------------------------------------- */

void Dislocation::deallocate()
{
  memory->destroy(x);
  memory->destroy(type);
  memory->destroy(q);
  memory->destroy(radius);
  memory->destroy(rmass);
  
  memory->destroy(num_bond);
  memory->destroy(bond_type);
  memory->destroy(bond_atom);
  
  memory->destroy(num_angle);
  memory->destroy(angle_type);
  memory->destroy(angle_atom1);
  memory->destroy(angle_atom2);
  memory->destroy(angle_atom3);
  
  memory->destroy(num_dihedral);
  memory->destroy(dihedral_type);
  memory->destroy(dihedral_atom1);
  memory->destroy(dihedral_atom2);
  memory->destroy(dihedral_atom3);
  memory->destroy(dihedral_atom4);
  
  memory->destroy(num_improper);
  memory->destroy(improper_type);
  memory->destroy(improper_atom1);
  memory->destroy(improper_atom2);
  memory->destroy(improper_atom3);
  memory->destroy(improper_atom4);

  memory->destroy(nspecial);
  memory->destroy(special);

  memory->destroy(shake_flag);
  memory->destroy(shake_atom);
  memory->destroy(shake_type);

  memory->destroy(dx);
  memory->destroy(dxcom);
  memory->destroy(dxbody);
  std::map<tagint, SlipVector*>::iterator it;
  for (it = burgersMap->begin(); it != burgersMap->end(); it++)
  {
      delete it->second;
  }
  delete burgersMap;
  delete slipSystem;
  delete atomIds;
  delete newAtomIds;

}

/* ----------------------------------------------------------------------
   open dislocation file
------------------------------------------------------------------------- */

void Dislocation::open(char *file)
{
  fp = fopen(file,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open dislocation file %s",file);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   read and bcast a line
------------------------------------------------------------------------- */

void Dislocation::readline(char *line)
{
  int n;
  if (me == 0) {
    if (fgets(line,MAXLINE,fp) == NULL) n = 0;
    else n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n == 0) error->all(FLERR,"Unexpected end of dislocation file");
  MPI_Bcast(line,n,MPI_CHAR,0,world);
}

/* ----------------------------------------------------------------------
   extract keyword from line
   flag = 0, read and bcast line
   flag = 1, line has already been read
------------------------------------------------------------------------- */

void Dislocation::parse_keyword(int flag, char *line, char *keyword)
{
  if (flag) {

    // read upto non-blank line plus 1 following line
    // eof is set to 1 if any read hits end-of-file

    int eof = 0;
    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
	if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      }
      if (fgets(keyword,MAXLINE,fp) == NULL) eof = 1;
    }

    // if eof, set keyword empty and return

    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) {
      keyword[0] = '\0';
      return;
    }

    // bcast keyword line to all procs

    int n;
    if (me == 0) n = strlen(line) + 1;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);
  }

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
         || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

/* ----------------------------------------------------------------------
   skip N lines of file
------------------------------------------------------------------------- */

void Dislocation::skip_lines(int n, char *line)
{
  for (int i = 0; i < n; i++) readline(line);
}

/* ----------------------------------------------------------------------
   parse line into words separated by whitespace
   return # of words
   max = max pointers storable in words
------------------------------------------------------------------------- */

int Dislocation::parse(char *line, char **words, int max)
{
  char *ptr;

  int nwords = 0;
  words[nwords++] = strtok(line," \t\n\r\f");

  while ((ptr = strtok(NULL," \t\n\r\f"))) {
    if (nwords < max) words[nwords] = ptr;
    nwords++;
  }

  return nwords;
}

/***********************************************************************
 Function name: add_atom
 Description:   add number of atoms in the list atomsIn into the dislocation
 Input:         int natomIn
                tagint* atomsIn
 Output:        None
 Return:        None
***********************************************************************/
void Dislocation::add_atom(int natomIn, tagint* atomsIn)
{
	newAtomIds->clear();
    for (int inx = 0; inx < natomIn; inx++)
    {
        atomIds->push_back(atomsIn[inx]);
		newAtomIds->push_back(atomsIn[inx]);
    }
    natoms = atomIds->size();
	//fprintf(logfile, "Add atom to dislocation %s:\n"
	//	"  %d atoms\n", id, natoms);
	//fprintf(logfile, "atoms list:\n");
	//for (int inx = 0; inx < natoms; inx++)
	//{
	//	fprintf(logfile, " %d ", (*atomIds)[inx]);
	//}
	//fprintf(logfile, "\n");
    //Need to recalculate the physical properties of dislocations
    //To be continue:
}

/***********************************************************************
Function name:  remove_atom
Description:    Remove number of atoms in the list atomsIn from the dislocation
Input:          int natomIn, 
                tagint* atomsIn
Output:         None
Return:         None
***********************************************************************/
void Dislocation::remove_atom(int natomIn, tagint* atomsIn)
{
    int n = 0;
    std::vector<tagint> temp;
    for (int inx = 0; inx < natomIn; inx++)
    {
        std::vector<tagint>::iterator it = std::find(atomIds->begin(), atomIds->end(), atomsIn[inx]);
        if (atomIds->end()!=it)
        {
            temp.push_back(atomsIn[inx]);
            atomIds->erase(it);
            n++;
        }
    }

	//if (n > 0){
	//	natoms -= n;
	//	fprintf(logfile, "Removed atom from dislocation %s:\n"
	//		"  %d atoms\n", id, n);
	//	fprintf(logfile, "atoms list:\n");
	//	for (int inx = 0; inx < n; inx++)
	//	{
	//		fprintf(logfile, " %d ", temp[inx]);
	//	}
	//	fprintf(logfile, "\n");
	//}
}

/***********************************************************************
Function name:  find_atom
Description:    Remove number of atoms in the list atomsIn from the dislocation
Input:          tagint atomsIn
Output:         None
Return:         int = 1 - the atom is in the dislocation
                    = 0 - the atom is not in the dislocation
***********************************************************************/
int Dislocation::find_atom(tagint atomsIn)
{
    return std::find(atomIds->begin(), atomIds->end(), atomsIn)!=atomIds->end();
}

/***********************************************************************
Function name:  CalcSlipSystem
Description:    Calculate the slip system of this dislocation
				The slip system consists of slip vector and slip plane
				which stored in the member slipSystem 
Input:          None
Output:         None
Return:         None
***********************************************************************/
void Dislocation::CalcSlipSystem()
{
	double* vx = domain->lattice->a1;
	double* vy = domain->lattice->a2;
	double* vz = domain->lattice->a3;

	double x[] = { FCCSlipVector[0][0], FCCSlipVector[0][1], FCCSlipVector[0][2] };
	//fprintf(logfile, "%f, %f, %f box to lattice :", x[0], x[1], x[2]);
	domain->lattice->box2lattice(x[0],x[1],x[2]);
	//fprintf(logfile, "%f, %f, %f \n", x[0], x[1], x[2]);
    tagint seedAtomId = *(atomIds->begin());
    int n[4] = {0}, index = -1, tempMax = 0;
    double myradius = sqrt(2)*domain->lattice->xlattice/4;
    double **tempX = atom->x;
    double seedAtom[] = { tempX[seedAtomId][0], tempX[seedAtomId][1], tempX[seedAtomId][1] };
    //cycle all predefined slip plane
    //Plane normal: A B C, plane through point X0 Y0 Z0
    //Plane A(X-X0)+B(Y-Y0)+C(Z-Z0)=0
    //      AX+BY+CZ+(-AX0-BY0-CZ0)=0  and D = -AX0-BY0-CZ0
    //Point (X1 Y1 Z1)
    //The distance from point to plane is
    //d = abs(AX1+BY1+CZ1+D)/sqrt(A^2+B^2+C^2)
    //if d < r, then this point (X1 Y1 Z1) is the satified point

    for (int np = 0; np < 4; np++)
    {
        double A[3] = { FCCSlipPlane[np][0], FCCSlipPlane[np][1], FCCSlipPlane[np][2] };
        double org[3] = { tempX[seedAtomId][0], tempX[seedAtomId][1], tempX[seedAtomId][1] };

        domain->lattice->lattice2box(A[0], A[1], A[2]);
        std::vector<tagint>::iterator it;
        std::vector<tagint> atomsInPlane;

        // cycle all atoms in the dislocation to check the distance 
        // to plane through the seed atom
        for (it = atomIds->begin(); it != atomIds->end();it++)
        {
            if ( *it != seedAtomId )
            {
                double point[3] = { tempX[*it][0], tempX[*it][1], tempX[*it][2] };
                double d = ProjecDistanceToPlane(A, org, point);

                if (d < myradius)
                {
                    n[np]++;
                    atomsInPlane.push_back(*it);
                }
            }
        }

        double xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
        double length[3];
        // Calculate the distance in this plane
        for (it = atomsInPlane.begin(); it != atomsInPlane.end();++it)
        {
            //project atom x y z to the plane
            double point[3] = { tempX[*it][0], tempX[*it][1], tempX[*it][2] };
            
            //slip vector
            double v1[3] = { FCCSlipVector[np * 3][0], FCCSlipVector[np * 3][1], FCCSlipVector[np * 3][2] };

            domain->lattice->lattice2box(v1[0], v1[1], v1[2]);

            double v2[3] = {0};
            double * temp = v2;

            productVector(A, v1, temp);

            double d1 = ProjecDistanceToPlane(A, org, point);
            double d2 = ProjecDistanceToPlane(v1, org, point);
            double d3 = ProjecDistanceToPlane(v2, org, point);

            bbox(d1,d2,d3,xmin,ymin,zmin,xmax,ymax,zmax);
            length[0] = xmax - xmin;
            length[1] = ymax - ymin;
            length[2] = zmax - zmin;

			maxLength = MAX(length[0], maxLength);
			maxLength = MAX(length[1], maxLength);
			maxLength = MAX(length[2], maxLength);
            //fprintf(logfile, "%s length is %f, %f, %f, Max length is %f\n", id, length[0], length[1], length[2], maxLength);
        }

        if ( tempMax < n[np] )
        {
            tempMax = n[np];
            index = np;
        }
    }

    if (-1 == index)
    {
        index = 0;
        //fprintf(logfile, "default slip plane is used since no slip plane found for %s\n", id);
    }

    SlipPlane* p = new SlipPlane(lmp, FCCSlipPlane[index][0], FCCSlipPlane[index][1], FCCSlipPlane[index][2]);
    SlipVector* v = new SlipVector(lmp, FCCSlipVector[index*3][0], FCCSlipVector[index*3][1], FCCSlipVector[index*3][2]);
    slipSystem = new SlipSystem(lmp, p, v);
    //fprintf(logfile, "%s slip plane is %f, %f, %f \n", id, p->u,p->v,p->w);
}

/***********************************************************************
Function name:  ProjecDistanceToPlane
Description:    Calculate project distance to a plane
Input:          double norm[3], The norm of the plane
                double org[3],  The point which the plane through
                double point[3], The point which will project to the plane
Output:         None
Return:         double
***********************************************************************/
double Dislocation::ProjecDistanceToPlane(double norm[3], double org[3], double point[3])
{
    double root = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
    return abs(norm[0] * (org[0] - point[0]) + norm[1] * (org[1] - point[1]) + norm[2] * (org[2] - point[2])) / root;
}

/***********************************************************************
Function name:  productVector
Description:    Compute the product of two vector
Input:          double v1[3], input vector 1
                double v2[3]��input vector 2
Output:         double *&v3,  output vector 3
Return:         None
***********************************************************************/
void Dislocation::productVector(double v1[3], double v2[3], double *&v3)
{
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

/***********************************************************************
Function name:  bbox
Description:    Calculate the min bound box
Input:          None
Output:         None
Return:         None
***********************************************************************/
void Dislocation::bbox(double x, double y, double z,
    double &xmin, double &ymin, double &zmin,
    double &xmax, double &ymax, double &zmax)
{
    xmin = MIN(x, xmin);  ymin = MIN(y, ymin);  zmin = MIN(z, zmin);
    xmax = MAX(x, xmax);  ymax = MAX(y, ymax);  zmax = MAX(z, zmax);
}

/***********************************************************************
Function name:  Length
Description:    Calculate the dislocation length
Input:          None
Output:         None
Return:         None
***********************************************************************/
void Dislocation::Length()
{
	double* vx = domain->lattice->a1;
	double* vy = domain->lattice->a2;
	double* vz = domain->lattice->a3;

	double x[] = { slipSystem->vector->u, slipSystem->vector->v, slipSystem->vector->w };
	//fprintf(logfile, "%f, %f, %f box to lattice :", x[0], x[1], x[2]);
	domain->lattice->box2lattice(x[0], x[1], x[2]);
	//fprintf(logfile, "%f, %f, %f \n", x[0], x[1], x[2]);
	tagint seedAtomId = *(newAtomIds->begin());
	int n[4] = { 0 }, index = -1, tempMax = 0;
	double myradius = sqrt(2)*domain->lattice->xlattice / 4;
	double **atomDisp = atom->x;
	double seedAtom[] = { atomDisp[seedAtomId][0], atomDisp[seedAtomId][1], atomDisp[seedAtomId][1] };
	//cycle all predefined slip plane
	//Plane normal: A B C, plane through point X0 Y0 Z0
	//Plane A(X-X0)+B(Y-Y0)+C(Z-Z0)=0
	//      AX+BY+CZ+(-AX0-BY0-CZ0)=0  and D = -AX0-BY0-CZ0
	//Point (X1 Y1 Z1)
	//The distance from point to plane is
	//d = abs(AX1+BY1+CZ1+D)/sqrt(A^2+B^2+C^2)
	//if d < r, then this point (X1 Y1 Z1) is the satified point

	double A[3] = { slipSystem->plane->u, slipSystem->plane->v, slipSystem->plane->w };
	double org[3] = { atomDisp[seedAtomId][0], atomDisp[seedAtomId][1], atomDisp[seedAtomId][1] };

	domain->lattice->lattice2box(A[0], A[1], A[2]);
	std::vector<tagint>::iterator it;
	std::vector<tagint> atomsInPlane;
	atomsInPlane.push_back(seedAtomId);
	// cycle all atoms in the dislocation to check the distance 
	// to plane through the seed atom
	for (it = newAtomIds->begin(); it != newAtomIds->end(); it++)
	{
		if (*it != seedAtomId)
		{
			double point[3] = { atomDisp[*it][0], atomDisp[*it][1], atomDisp[*it][2] };
			double d = ProjecDistanceToPlane(A, org, point);

			if (d < myradius)
			{
				atomsInPlane.push_back(*it);
			}
		}
	}

	double xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
	double length[3];
	// Calculate the distance in this plane
	for (it = atomsInPlane.begin(); it != atomsInPlane.end(); ++it)
	{
		double xx = atomDisp[*it][0];
		double yy = atomDisp[*it][1];
		double zz = atomDisp[*it][2];

		double deltx = org[0] - xx;
		double delty = org[1] - yy;
		double deltz = org[2] - zz;

		double distance = sqrt(deltx*deltx+delty*delty+deltz*deltz);
		dislocationLength = MAX(distance, dislocationLength);
	}
	//fprintf(logfile, "%s length is %f, %f, %f, Max length is %f\n", id, length[0], length[1], length[2], dislocationLength);
}

