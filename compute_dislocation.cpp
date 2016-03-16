/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Jianfeng Huang Glasgow Caledonian University
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "compute_dislocation.h"
#include "dislocation.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include <algorithm>
#include "timer.h"
#include "domain.h"
#include "lattice.h"

using namespace LAMMPS_NS;

#define MAXNEAR 16
#define MAXCOMMON 8

enum{UNKNOWN,FCC,HCP,BCC,ICOS,OTHER};
enum{NCOMMON,NBOND,MAXBOND,MINBOND};

/* ---------------------------------------------------------------------- */

ComputeDislocation::ComputeDislocation(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute dislocation command");

  peratom_flag = 1;
  size_peratom_cols = 5;

  double cutoff = force->numeric(FLERR,arg[3]);
  if (cutoff < 0.0) error->all(FLERR,"Illegal compute dislocation command");
  cutsq = cutoff*cutoff;

  numDislocation = 0;
  nmax = 0;
  nearest = NULL;
  nnearest = NULL;
  pattern = NULL;
  dislocations = NULL;
  displacements = NULL;
  atomSlipVector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDislocation::~ComputeDislocation()
{
  memory->destroy(nearest);
  memory->destroy(nnearest);
  memory->destroy(pattern);
  memory->destroy(displacements);
  memory->destroy(atomSlipVector);
  
  for (int inx = 0; inx < numDislocation; inx++)
  {
	  delete dislocations[inx];
	  dislocations[inx] = NULL;
  }
  memory->destroy(dislocations);
  dislocations = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeDislocation::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute dislocation requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute dislocation cutoff is longer than pairwise cutoff");

  // cannot use neighbor->cutneighmax b/c neighbor has not yet been init

  if (2.0*sqrt(cutsq) > force->pair->cutforce + neighbor->skin &&
      comm->me == 0)
    error->warning(FLERR,"Compute dislocation cutoff may be too large to find "
                   "ghost atom neighbors");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"dislocation") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute dislocation defined");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeDislocation::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeDislocation::compute_peratom()
{
  int i,j,k,ii,jj,kk,m,n,inum,jnum,inear,jnear;
  int firstflag,ncommon,nbonds,maxbonds,minbonds;
  int nfcc,nhcp,nbcc4,nbcc6,nico,cj,ck,cl,cm;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int cna[MAXNEAR][4],onenearest[MAXNEAR];
  int common[MAXCOMMON],bonds[MAXCOMMON];
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  std::vector<tagint> candidateAtom;

  invoked_peratom = update->ntimestep;

  // grow arrays if necessary
  int flag = 0;
  if (atom->nlocal > nmax) {
    memory->destroy(nearest);
    memory->destroy(nnearest);
    memory->destroy(pattern);
    memory->destroy(displacements);
    memory->destroy(atomSlipVector);
    nmax = atom->nmax;

    memory->create(nearest,nmax,MAXNEAR,"dislocation:nearest");
    memory->create(nnearest,nmax,"dislocation:nnearest");
    memory->create(pattern,nmax,"dislocation:dislocation_pattern");
	memory->create(atomSlipVector, nmax, 5, "dislocation:dislocation_slip_vector");
    array_atom = atomSlipVector;
    vector_atom = pattern;

	for (int inx = 0; inx < nmax;inx++)
	{
        atomSlipVector[inx][0] = 0.0;
        atomSlipVector[inx][1] = 0.0;
        atomSlipVector[inx][2] = 0.0;
        atomSlipVector[inx][3] = 0.0;
		atomSlipVector[inx][4] = 0.0;
	}
	memory->create(displacements, nmax, 3, "Dislocation::Atom::Displacement");
    flag = 1;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // find the neigbours of each atom within cutoff using full neighbor list
  // nearest[] = atom indices of nearest neighbors, up to MAXNEAR
  // do this for all atoms, not just compute group
  // since dislocation calculation requires neighbors of neighbors

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int nerror = 0;
  for (ii = 0; ii < inum; ii++) {
    //std::vector<tagint>* neighbor = new std::vector<tagint>();

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    if (flag)
    {
        displacements[i][0] = x[i][0];
        displacements[i][1] = x[i][1];
        displacements[i][2] = x[i][2];
    }

    jlist = firstneigh[i];
    jnum = numneigh[i];

    n = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        if (n < MAXNEAR) {
              //neighbor->push_back(j);
              nearest[i][n++] = j;
        }
        else {
          nerror++;
          break;
        }
      }
    }

    nnearest[i] = n;
  }

  // warning message

  int nerrorall;
  MPI_Allreduce(&nerror,&nerrorall,1,MPI_INT,MPI_SUM,world);
  if (nerrorall && comm->me == 0) {
    char str[128];
    sprintf(str,"Too many neighbors in dislocation for %d atoms",nerrorall);
    error->warning(FLERR,str,0);
  }

  // compute dislocation for each atom in group
  // only performed if # of nearest neighbors = 12 or 14 (fcc,hcp)

  nerror = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (!(mask[i] & groupbit)) {
      pattern[i] = UNKNOWN;
      continue;
    }

    if (nnearest[i] != 12 && nnearest[i] != 14) {
      pattern[i] = OTHER;
      continue;
    }

    // loop over near neighbors of I to build dislocation data structure
    // cna[k][NCOMMON] = # of common neighbors of I with each of its neighs
    // cna[k][NBONDS] = # of bonds between those common neighbors
    // cna[k][MAXBOND] = max # of bonds of any common neighbor
    // cna[k][MINBOND] = min # of bonds of any common neighbor

    for (m = 0; m < nnearest[i]; m++) {
      j = nearest[i][m];

      // common = list of neighbors common to atom I and atom J
      // if J is an owned atom, use its near neighbor list to find them
      // if J is a ghost atom, use full neighbor list of I to find them
      // in latter case, must exclude J from I's neighbor list

      if (j < nlocal) {
        firstflag = 1;
        ncommon = 0;
        for (inear = 0; inear < nnearest[i]; inear++)
          for (jnear = 0; jnear < nnearest[j]; jnear++)
            if (nearest[i][inear] == nearest[j][jnear]) {
              if (ncommon < MAXCOMMON) common[ncommon++] = nearest[i][inear];
              else if (firstflag) {
                nerror++;
                firstflag = 0;
              }
            }

      } else {
        xtmp = x[j][0];
        ytmp = x[j][1];
        ztmp = x[j][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        n = 0;
        for (kk = 0; kk < jnum; kk++) {
          k = jlist[kk];
          k &= NEIGHMASK;
          if (k == j) continue;

          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            if (n < MAXNEAR) onenearest[n++] = k;
            else break;
          }
        }

        firstflag = 1;
        ncommon = 0;
        for (inear = 0; inear < nnearest[i]; inear++)
          for (jnear = 0; jnear < n; jnear++)
            if (nearest[i][inear] == onenearest[jnear]) {
              if (ncommon < MAXCOMMON) common[ncommon++] = nearest[i][inear];
              else if (firstflag) {
                nerror++;
                firstflag = 0;
              }
            }
      }

      cna[m][NCOMMON] = ncommon;

      // calculate total # of bonds between common neighbor atoms
      // also max and min # of common atoms any common atom is bonded to
      // bond = pair of atoms within cutoff

      for (n = 0; n < ncommon; n++) bonds[n] = 0;

      nbonds = 0;
      for (jj = 0; jj < ncommon; jj++) {
        j = common[jj];
        xtmp = x[j][0];
        ytmp = x[j][1];
        ztmp = x[j][2];
        for (kk = jj+1; kk < ncommon; kk++) {
          k = common[kk];
          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            nbonds++;
            bonds[jj]++;
            bonds[kk]++;
          }
        }
      }

      cna[m][NBOND] = nbonds;

      maxbonds = 0;
      minbonds = MAXCOMMON;
      for (n = 0; n < ncommon; n++) {
        maxbonds = MAX(bonds[n],maxbonds);
        minbonds = MIN(bonds[n],minbonds);
      }
      cna[m][MAXBOND] = maxbonds;
      cna[m][MINBOND] = minbonds;
    }

    // detect CNA pattern of the atom

    nfcc = nhcp = nbcc4 = nbcc6 = nico = 0;
    pattern[i] = OTHER;
	atomSlipVector[i][0] = OTHER;
    if (nnearest[i] == 12) {
      for (inear = 0; inear < 12; inear++) {
        cj = cna[inear][NCOMMON];
        ck = cna[inear][NBOND];
        cl = cna[inear][MAXBOND];
        cm = cna[inear][MINBOND];
        if (cj == 4 && ck == 2 && cl == 1 && cm == 1) nfcc++;
        else if (cj == 4 && ck == 2 && cl == 2 && cm == 0) nhcp++;
        else if (cj == 5 && ck == 5 && cl == 2 && cm == 2) nico++;
      }
      if (nfcc == 12) pattern[i] = FCC;
      else if (nfcc == 6 && nhcp == 6){
          // add this atom to the candidate list
          candidateAtom.push_back((tagint)i);
          pattern[i] = HCP;
      }
      else if (nico == 12) pattern[i] = ICOS;

    } else if (nnearest[i] == 14) {
      for (inear = 0; inear < 14; inear++) {
        cj = cna[inear][NCOMMON];
        ck = cna[inear][NBOND];
        cl = cna[inear][MAXBOND];
        cm = cna[inear][MINBOND];
        if (cj == 4 && ck == 4 && cl == 2 && cm == 2) nbcc4++;
        else if (cj == 6 && ck == 6 && cl == 2 && cm == 2) nbcc6++;
      }
      if (nbcc4 == 6 && nbcc6 == 8) pattern[i] = BCC;
    }
  }


  std::vector<tagint>::iterator seedIt;

  // Flat all dislocation atoms in one vector;
  // match the dislocation atoms in this vector to the candidate atoms.
  // Search candidate atoms which not in any dislocation and form a new candidate list
  // cycle the new candidate and if it adjacent to one of dislocation,
  // then add it to this dislocation, otherwise try to create a new dislocation
  // with the seed atom.
  std::vector<tagint> existDislocationAtoms;
  flat_dislocation_atoms(existDislocationAtoms);

  std::sort(existDislocationAtoms.begin(),existDislocationAtoms.end());
  std::sort(candidateAtom.begin(), candidateAtom.end());
  std::vector<tagint> newCandidate(candidateAtom.size());
  std::vector<tagint> removeCandidate(existDislocationAtoms.size());

  // Now get the difference in candidate atom list and current exist dislocation list.
  // Generate the new candidate and remove candidate lists.
  // For example, if atoms (1,2,3,4,5) are in exist dislocation list, (4,5,6,7,8) are
  // the candidate atoms. Then, the new candidate are (6 7 8), remove candidate are
  // atoms (1,2,3)
  if ( existDislocationAtoms.size() > 2 )
  {
	  std::set_difference(candidateAtom.begin(), candidateAtom.end(), existDislocationAtoms.begin(), existDislocationAtoms.end(), newCandidate.begin());
      std::set_difference(existDislocationAtoms.begin(), existDislocationAtoms.end(), candidateAtom.begin(), candidateAtom.end(), removeCandidate.begin());
  }
  else
  {
	  newCandidate = candidateAtom;
  }

  std::vector<tagint> temp=newCandidate;
  std::vector<tagint>::iterator tempIt;
  //cycle all new candidate atoms
  for (seedIt = newCandidate.begin(); seedIt != newCandidate.end() && *seedIt != 0; seedIt++)
  {
      flag = 0;

      if (temp.size() <= 1)
          break;

      std::vector<tagint> dislocationAtom;
      
      if (std::find(temp.begin(), temp.end(), *seedIt) == temp.end())
          continue;

      find_neighbor_dislocation_atom(*seedIt, newCandidate, dislocationAtom);
      if (dislocationAtom.size() >= 3)
      {
          std::vector<tagint>::iterator disIt;
          if (numDislocation > 0)
          {
              // Check if this atom already in any dislocation
              for (int jj = 0; jj < numDislocation; jj++)
              {
                  for (disIt = dislocationAtom.begin(); disIt != dislocationAtom.end(); disIt++)
                  {
                      if (is_atom_next_to_dislocation(*disIt, dislocations[jj]) == 1)
                      {
                          flag = 1;
                          break;
                      }
                  }

                  if (flag == 1)
                  {
                      add_atoms_to_dislocation(dislocationAtom.size(), dislocationAtom.data(), dislocations[jj]);
                      for (disIt = dislocationAtom.begin(); disIt != dislocationAtom.end(); disIt++)
                      {
                          std::vector<tagint>::iterator eraseIt = 
                              std::find(temp.begin(), temp.end(), *disIt);
                          if (eraseIt != temp.end())
                          {
                              temp.erase(eraseIt);
                          }
                      }
                      temp.swap(temp);
                      break;
                  }
              }
          }

          if (flag == 0)
          {
              if (numDislocation == 0)
              {
                  numDislocation++;
                  memory->create(dislocations, numDislocation, sizeof(Dislocation*), "dislocation");
              }
              else
              {
                  dislocations = (Dislocation **)
                      memory->srealloc(dislocations, ++numDislocation*sizeof(Dislocation *), "dislocation");
              }

              for (disIt = dislocationAtom.begin(); disIt != dislocationAtom.end(); disIt++)
              {
                  std::vector<tagint>::iterator eraseIt = 
                      std::find(temp.begin(), temp.end(), *disIt);
                  if (eraseIt != temp.end())
                  {
                      temp.erase(eraseIt);
                  }
                  temp.swap(temp);
              }
			  
              dislocations[numDislocation - 1] = new Dislocation(lmp, dislocationAtom.size(), dislocationAtom.data());

			  atom->add_dislocation(dislocations[numDislocation - 1]);
          }
          else
              flag = 0;
      }
  }

  for (int inx = 0; inx < numDislocation; inx++)
  {
	  dislocations[inx]->remove_atom(removeCandidate.size(), removeCandidate.data());
	  if (dislocations[inx]->natoms < 3)
	  {
		  fprintf(logfile, "Deleted dislocation %s:\n", dislocations[inx]->id);

		  atom->remove_dislocation(dislocations[inx]);

		  delete dislocations[inx];

		  for (int jnx = inx; jnx < numDislocation; jnx++)
		  {
			  if (jnx == numDislocation - 1)
			  {
				  dislocations[jnx] = NULL;
			  }
			  else
				  dislocations[jnx] = dislocations[jnx + 1];
		  }
		  numDislocation--;
		  inx--;
	  }

	  Burgers(dislocations[inx]);
	  Length(dislocations[inx]);
	  int nAtom = dislocations[inx]->natoms;
	  std::vector<tagint>::iterator temp;
	  for (temp = dislocations[inx]->atomIds->begin(); temp < dislocations[inx]->atomIds->end(); temp++)
	  {
		  atomSlipVector[*temp][0] = OTHER + inx + 1;
	  }
  }

  if ( numDislocation > 1)
	atom->dislocation = 1;

  // warning message

  MPI_Allreduce(&nerror,&nerrorall,1,MPI_INT,MPI_SUM,world);
  if (nerrorall && comm->me == 0) {
    char str[128];
    sprintf(str,"Too many common neighbors in CNA %d times",nerrorall);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDislocation::memory_usage()
{
  double bytes = nmax * sizeof(int);
  bytes += nmax * MAXNEAR * sizeof(int);
  bytes += nmax * sizeof(double);
  return bytes;
}

/***********************************************************************
Function Name: add_atom_to_dislocation
Description:   Add any atom to a dislocation. This function can tell any
               atom whether should add to an existing dislocation or not.
               It will check the atom is near the dislocation or not. If
               true, it will add the atom to the dislocation. otherwise
               leave it alone.
Input:         tagint seedAtom - the atom id.
               Dislocation *dislocation - the dislocation
output:        None
return:        int 0 = No, not add the seed atom into the dislocation
                   1 = Yes, added the atom to the dislocation
                   2 = atom already in this dislocation
***********************************************************************/
int ComputeDislocation::add_atom_to_dislocation(tagint seedAtom, Dislocation * dislocation)
{
    int ret = is_atom_next_to_dislocation(seedAtom, dislocation);
    if (ret == 1)
	{
		dislocation->add_atom(1, &seedAtom);
	}
    return ret;
}

/***********************************************************************
Function Name: is_atom_next_to_dislocation
Description:   function to determinate whether the atom is near a dislocation 
Input:         tagint seedAtom - the atom id.
               Dislocation *dislocation - the dislocation
output:        None
return:        int 0 = No, the atom is not near the dislocation
                   1 = Yes
                   2 = atom already in this dislocation
***********************************************************************/
int ComputeDislocation::is_atom_next_to_dislocation(tagint seedAtom, Dislocation * dislocation)
{
	std::vector<tagint> nearest;
	std::vector<tagint>::iterator nearestIt;

    find_nearest_neighbor(seedAtom, nearest);
    if (1==dislocation->find_atom(seedAtom))
        return 2;

    for (nearestIt = nearest.begin(); nearestIt != nearest.end(); nearestIt++)
    {
        if (1==dislocation->find_atom(*nearestIt))
        {
            return 1;
        }
    }
    return 0;
}

/***********************************************************************
Function Name: add_atoms_to_dislocation
Description:   Add atoms to a dislocation
Input:         int numAtom - the number of atoms.
               tagint* atomIdIn - atom Ids list will be added into a dislocation
               Dislocation *dislocation - the dislocation
output:        None
Return:        None
***********************************************************************/
void ComputeDislocation::add_atoms_to_dislocation(int numAtom, tagint* atomIdIn, Dislocation *dislocation)
{
    dislocation->add_atom(numAtom, atomIdIn);
}

/***********************************************************************
Function Name: find_neighbor_dislocation_atom
Description:   Find all of the nearest neighbor atoms of the input seed 
               atom. All of those atoms will be connected and are the
               core part of a dislocation
Input:         int seedAtom - the seed atom id.
               std::vector<tagint>& candidate - the candidate atom list
output:        std::vector<tagint>& dislocationAtom - all the dislocation
               atoms
Return:        None
***********************************************************************/
void ComputeDislocation::find_neighbor_dislocation_atom(int seedAtom,
    std::vector<tagint>& candidate,
    std::vector<tagint>& dislocationAtom)
{
	std::vector<tagint> nearest;
	std::vector<tagint>::iterator nearestIt;
	find_nearest_neighbor(seedAtom, nearest);
	if (1==candidate.size())
		return;

	std::vector<tagint> existDislocationAtoms;
	flat_dislocation_atoms(existDislocationAtoms);

	std::vector<tagint>::iterator seedIt = find(candidate.begin(),candidate.end(), seedAtom );
	for (nearestIt = nearest.begin(); nearestIt != nearest.end();nearestIt++)
    {
		if (1==candidate.size())
			break;

		std::vector<tagint>::iterator it = find(candidate.begin(), candidate.end(), *nearestIt);
		std::vector<tagint>::iterator aIt = existDislocationAtoms.end();
		if ( existDislocationAtoms.size() > 0 )
		{
			aIt = find(existDislocationAtoms.begin(), existDislocationAtoms.end(), *nearestIt);
		}
		if (it != candidate.end() && aIt == existDislocationAtoms.end())
        {
            if (0==dislocationAtom.size())
            {
				dislocationAtom.push_back(seedAtom);
				dislocationAtom.push_back(*nearestIt);
				//candidate.erase(seedIt);
				find_neighbor_dislocation_atom(*nearestIt, candidate, dislocationAtom);
            }
			else
			{
				std::vector<tagint>::iterator dIt = find(dislocationAtom.begin(), dislocationAtom.end(), *nearestIt);
                if (dislocationAtom.end() == dIt )
				{
					dislocationAtom.push_back(*nearestIt);
					//candidate.erase(seedIt);
					find_neighbor_dislocation_atom(*nearestIt, candidate, dislocationAtom);
				}
			}
        }
    }
}

/***********************************************************************
Function Name: flat_dislocation_atoms
Description:   Get atoms which is in current existing dislocations and 
               store them in an vector
Input:         None
output:        std::vector<tagint>& atomsOut - the dislocation atoms
Return:        None
***********************************************************************/
void ComputeDislocation::flat_dislocation_atoms(std::vector<tagint>& atomsOut)
{
    std::vector<tagint>::iterator it;

    for (int inx = 1; inx < numDislocation; inx++)
    {
        for (it = dislocations[inx]->atomIds->begin(); it != dislocations[inx]->atomIds->end(); ++it)
        {
            atomsOut.push_back(*it);
        }
	}
}

/***********************************************************************
Function Name: find_nearest_neighbor
Description:   Find the nearest neighbors of an atom
Input:         tagint seedAtom - seed atom
               std::vector<tagint>& nearest - the nearest atoms list
output:        None
Return:        None
***********************************************************************/
void ComputeDislocation::find_nearest_neighbor(tagint seedAtom, std::vector<tagint>& nearestVec)
{
	if (nearest != NULL)
	{
		for (int inx = 0; inx < nnearest[seedAtom]; inx++)
		{
			nearestVec.push_back(nearest[seedAtom][inx]);
		}
		return;
	}

	int i, j, jj, n, jnum, jnear;
	int *jlist, *numneigh, **firstneigh;
	double xtmp, ytmp, ztmp, delx, dely, delz, rsq;

	invoked_peratom = update->ntimestep;

	// invoke full neighbor list (will copy or build if necessary)

	neighbor->build_one(list);

	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

	// find the neigbours of each atom within cutoff using full neighbor list
	// nearest[] = atom indices of nearest neighbors, up to MAXNEAR
	// do this for all atoms, not just compute group
	// since CNA calculation requires neighbors of neighbors

	double **x = atom->x;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;

	int nerror = 0;
	i = seedAtom;
	xtmp = x[i][0];
	ytmp = x[i][1];
	ztmp = x[i][2];
	jlist = firstneigh[i];
	jnum = numneigh[i];
	n = 0;
	for (jj = 0; jj < jnum; jj++) {
		j = jlist[jj];
		j &= NEIGHMASK;

		delx = xtmp - x[j][0];
		dely = ytmp - x[j][1];
		delz = ztmp - x[j][2];
		rsq = delx*delx + dely*dely + delz*delz;
		if (rsq < cutsq) {
			if (n < MAXNEAR) {
				n++;
				nearestVec.push_back(j);
			}
			else {
				nerror++;
				break;
			}
		}
	}
}

/***********************************************************************
Function Name: burgers
Description:   Calculate the burgers vector of dislocation
               the burgers vector stored in burgers vector map
Input:         Dislocation* dislocationIn
output:        None
Return:        None
***********************************************************************/
void ComputeDislocation::Burgers(Dislocation* dislocationIn)
{
    int num = dislocationIn->natoms;
    double ** x = atom->x;
    std::vector<tagint>::iterator disAtomIt;
    double deltXij[3];
    int n = 1;
	int nerror = 0;
	std::vector<tagint> temp1, temp2;

    //Get the dislocation atoms. namely atom A.
    for (disAtomIt = dislocationIn->atomIds->begin(); disAtomIt != dislocationIn->atomIds->end(); disAtomIt++)
    {
		n = 1;
		int flag = 0;
		nerror = 0;
        deltXij[0] = 0.0;
        deltXij[1] = 0.0;
        deltXij[2] = 0.0;
        
        //Get the atoms ID
        tagint atomId = *disAtomIt;

        if ( atomId == 0 )
            continue;

		std::vector<tagint> catoms;
        //Get the nearest neighbors of this atom
        find_nearest_neighbor(atomId, catoms);

		double xiRef = displacements[atomId][0];
		double yiRef = displacements[atomId][1];
		double ziRef = displacements[atomId][2];

        double xi = x[atomId][0];
        double yi = x[atomId][1];
        double zi = x[atomId][2];

        //Cycle all the nearest neighbors
        std::vector<tagint>::iterator catomsIt;
        for (catomsIt = catoms.begin(); catomsIt != catoms.end(); catomsIt++)
        {
            //if the neighbor is in the dislocation lists. catch this atom and calculate
			if (0 != dislocationIn->find_atom(*catomsIt) && *disAtomIt != *catomsIt )
            {
                //Now get the reference vector of the two atoms in the first step
                double xjRef = displacements[*catomsIt][0];
                double yjRef = displacements[*catomsIt][1];
                double zjRef = displacements[*catomsIt][2];

                double xj = x[*catomsIt][0];
                double yj = x[*catomsIt][1];
                double zj = x[*catomsIt][2];

				double delx = xi - xj;
				double dely = yi - yj;
				double delz = zi - zj;
				double rsq = delx*delx + dely*dely + delz*delz;

				if (rsq > cutsq)
				{
					nerror++;
				}
				else
				{
                    // The distance between two atoms is smaller than cutoff distance
                    // then catch the atom A.
					std::vector<tagint>::iterator it;
					if (0==flag)
					{
						it = std::find(temp2.begin(), temp2.end(), *disAtomIt);
						if (it == temp2.end())
						{
							temp2.push_back(*disAtomIt);
							flag = 1;
						}
					}

                    // catch the second atom. namely atom B.
					it = std::find(temp2.begin(), temp2.end(), *catomsIt);
					if ( it == temp2.end() )
					{
						temp2.push_back(*catomsIt);
					}

                    // calculate the accumulate different of the vector of 
                    // atom AB at current step between the reference step
                    // without deformation.

					//validate delx dely delz and reference delx dely delz 
					double rx, ry, rz;
					rx = xiRef - xjRef;
					ry = yiRef - yjRef;
					rz = ziRef - zjRef;

					if (fabs(delx) < domain->lattice->xlattice &&
						fabs(rx) < domain->lattice->xlattice &&
						fabs(dely) < domain->lattice->ylattice&&
						fabs(ry) < domain->lattice->ylattice&&
						fabs(delz) < domain->lattice->zlattice&&
						fabs(rz) < domain->lattice->zlattice)
					{
						deltXij[0] += delx - (xiRef - xjRef);
						deltXij[1] += dely - (yiRef - yjRef);
						deltXij[2] += delz - (ziRef - zjRef);

						n++;
					}
				}
            }
        }

        // output value (will dump out of lammps)
        atomSlipVector[atomId][1] = -deltXij[0] / n;
        atomSlipVector[atomId][2] = -deltXij[1] / n;
        atomSlipVector[atomId][3] = -deltXij[2] / n;

		double mag = sqrt(atomSlipVector[atomId][1] * atomSlipVector[atomId][1] + atomSlipVector[atomId][2] * atomSlipVector[atomId][2] + atomSlipVector[atomId][3] * atomSlipVector[atomId][3]);

		atomSlipVector[atomId][4] = mag;

        // store to burgers map
        std::map<tagint, SlipVector*>::iterator mapIt = dislocationIn->burgersMap->find(atomId);
        if (dislocationIn->burgersMap->end() == mapIt )
        {
            // create the slip vector if not exist
            SlipVector* slipVector = new SlipVector(lmp, deltXij[0]/n, deltXij[1]/n, deltXij[2]/n);
            dislocationIn->burgersMap->insert(std::pair<tagint, SlipVector*>(atomId, slipVector));
        }
        else
        {
            // update the slip vector if it already exist.
            mapIt->second->u = -deltXij[0]/n;
            mapIt->second->v = -deltXij[1]/n;
            mapIt->second->w = -deltXij[2]/n;
        }
    }
}

/***********************************************************************
Function Name: Length
Description:   Calculate the length of dislocation
Input:         Dislocation* dislocationIn
output:        None
Return:        None
***********************************************************************/
void ComputeDislocation::Length(Dislocation* dislocationIn)
{
	double* vx = domain->lattice->a1;
	double* vy = domain->lattice->a2;
	double* vz = domain->lattice->a3;

	double x[] = { dislocationIn->slipSystem->vector->u, dislocationIn->slipSystem->vector->v, dislocationIn->slipSystem->vector->w };
	fprintf(logfile, "%f, %f, %f box to lattice :", x[0], x[1], x[2]);
	domain->lattice->box2lattice(x[0], x[1], x[2]);
	fprintf(logfile, "%f, %f, %f \n", x[0], x[1], x[2]);
	tagint seedAtomId = *(dislocationIn->newAtomIds->begin());
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

	double A[3] = { dislocationIn->slipSystem->plane->u, dislocationIn->slipSystem->plane->v, dislocationIn->slipSystem->plane->w };

	domain->lattice->lattice2box(A[0], A[1], A[2]);
	std::vector<tagint>::iterator it;
	int numNeighbor = 2;

	std::vector<tagint> dislocationAtoms = *(dislocationIn->newAtomIds);
	std::sort(dislocationAtoms.begin(), dislocationAtoms.end());
	for (it = dislocationAtoms.begin(); it != dislocationAtoms.end(); ++it)
	{
		std::vector<tagint> nearestAtoms;
		find_nearest_neighbor(*it, nearestAtoms);
		std::sort(nearestAtoms.begin(), nearestAtoms.end());

		std::vector<tagint> tempV;
		std::set_intersection(nearestAtoms.begin(), nearestAtoms.end(), dislocationAtoms.begin(), dislocationAtoms.end(), back_inserter(tempV));

		if (tempV.size() < numNeighbor)
		{
			numNeighbor = tempV.size();
			seedAtomId = *it;
		}
	}

	std::vector<tagint> atomsInLine;
	double org[3] = { atomDisp[seedAtomId][0], atomDisp[seedAtomId][1], atomDisp[seedAtomId][1] };

	// cycle all atoms in the dislocation to check the distance 
	// to plane through the seed atom

	CalcLengthAtoms(seedAtomId, dislocationAtoms, atomsInLine, dislocationIn);

	double xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;

	// Calculate the distance in this plane
	for (it = atomsInLine.begin(); it != atomsInLine.end(); ++it)
	{
		double xx = atomDisp[*it][0];
		double yy = atomDisp[*it][1];
		double zz = atomDisp[*it][2];

		double deltx = org[0] - xx;
		double delty = org[1] - yy;
		double deltz = org[2] - zz;
		org[0] = xx;
		org[1] = yy;
		org[2] = zz;
		double distance = sqrt(deltx*deltx + delty*delty + deltz*deltz);

		if (distance < cutsq)
		{
			dislocationIn->addLength(distance);
		}
	}
	fprintf(logfile, "%s Max length is %f\n", dislocationIn->id, dislocationIn->dislocationLength);
}

/***********************************************************************
Function Name: CalcLengthAtoms
Description:   Calculate the length atoms of dislocation
Input:         Dislocation* dislocationIn
output:        None
Return:        None
***********************************************************************/
void ComputeDislocation::CalcLengthAtoms(tagint seedAtomId, 
	std::vector<tagint> dislocationAtoms, 
	std::vector<tagint>& nearestInPlaneVec, 
	Dislocation* dislocationIn)
{
	std::vector<tagint> nearestAtoms;
	find_nearest_neighbor(seedAtomId, nearestAtoms);
	std::vector<tagint>::iterator it;

	double myradius = sqrt(2)*domain->lattice->xlattice / 4;
	double **atomDisp = atom->x;
	double A[3] = { dislocationIn->slipSystem->plane->u, dislocationIn->slipSystem->plane->v, dislocationIn->slipSystem->plane->w };
	domain->lattice->lattice2box(A[0], A[1], A[2]);
	double org[3] = { atomDisp[seedAtomId][0], atomDisp[seedAtomId][1], atomDisp[seedAtomId][1] };
	double V[3] = { dislocationIn->slipSystem->vector->u, dislocationIn->slipSystem->vector->v, dislocationIn->slipSystem->vector->w };
	domain->lattice->lattice2box(V[0], V[1], V[2]);

	std::sort(nearestAtoms.begin(), nearestAtoms.end());
	std::sort(dislocationAtoms.begin(), dislocationAtoms.end());

	std::vector<tagint> nearestAtomIn;
	std::set_intersection(nearestAtoms.begin(), nearestAtoms.end(), dislocationAtoms.begin(), dislocationAtoms.end(), back_inserter(nearestAtomIn));

	for (it = nearestAtomIn.begin(); it != nearestAtomIn.end(); it++)
	{
		if (*it != seedAtomId)
		{
			double point[3] = { atomDisp[*it][0], atomDisp[*it][1], atomDisp[*it][2] };
			double d = dislocationIn->ProjecDistanceToPlane(A, org, point);

			if (d < myradius && std::find(nearestInPlaneVec.begin(), nearestInPlaneVec.end(), *it) == nearestInPlaneVec.end())
			{
				double xx = atomDisp[*it][0];
				double yy = atomDisp[*it][1];
				double zz = atomDisp[*it][2];

				double deltx = org[0] - xx;
				double delty = org[1] - yy;
				double deltz = org[2] - zz;

				double interProduct = deltx*V[0] + delty*V[1] + deltz*V[2];
				//if (interProduct == 0)
				{
					nearestInPlaneVec.push_back(*it);
					CalcLengthAtoms(*it, dislocationAtoms, nearestInPlaneVec, dislocationIn);
				}
			}
		}
	}
}
