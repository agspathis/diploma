#ifndef LP_GRID_H
#define LP_GRID_H

#include <btBulletDynamicsCommon.h>
#include <vector>

#include "terrain.h"
#include "fluid.h"

/* Structure LP_GRID: Data structure to store pointers to the simulation particles,
   preserving locality between simulation and memory space.
   Member documentation:
   ORIGIN	   : Coordinates of grid origin (minimum-coordinate point of AABB).
   STEP		   : Grid consists of adjacent cubic cells with side length equal to STEP.
   X, Y, Z	   : Number of cells along axes x, y anx z, respectively.
   CELL_COUNT	   : Total number of cells (X*Y*Z).
   PARTICLE_COUNT  : Total number of particles.
   MAP             : Array map of linearized addresses (3D->1D by function LINEARIZE_ID) to
                     linear indices (on ANCHOR_MAP), thus implementing the locality-preserving
		     3D->1D mapping.
   ANCHORS	   : Array containing the cell ANCHORs, i.e. the linearized address of the
		     cell in MAP (CELL_LADDRESS) and the index of the first of cell's
                     particles in PARTICLES (CELL_START) . ANCHORs correspond to cells in
		     ascending linear order, so the relation CS[i-1] <= CS[i] (CS=CELL_START)
		     holds for all elements of this array. When CS[i-1]=CS[i], the cell 
		     corresponding to CS[i-1] contains no particles. This array has
		     CELL_COUNT+1 elements, where the last anchor has CELL_START=PARTICLE_COUNT+1
		     and LADDRESS=-1 for loop termination conditions.
   PARTICLES       : Array containing pointers to actual simulation particles
*/

struct particle_range {
    particle* start;
    particle* end;
};

struct anchor {
    unsigned int cell_laddress;
    unsigned int cell_start;
};

struct lp_grid {
    btVector3 origin;
    float step;
    unsigned int x, y, z, cell_count, particle_count;
    unsigned int* map;
    anchor* anchors;
    particle** particles;
};

lp_grid make_lp_grid (btVector3 origin, float step,
		      unsigned int x, unsigned int y, unsigned int z,
		      std::vector<particle> particles);

int update_lp_grid (lp_grid lpg);

particle_range get_cell(lp_grid lpg, unsigned int i, unsigned int j, unsigned int k);

#endif
