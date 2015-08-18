#ifndef LP_GRID_H
#define LP_GRID_H

#include <btBulletDynamicsCommon.h>
#include <vector>

#include "terrain.h"
#include "fluid.h"

struct anchor {
    unsigned int linear_id;
    index3 grid_id;
};

/* Structure LP_GRID: Spatial hash table to store pointers to the simulation particles.
   Member documentation:
   ORIGIN	: Coordinates of grid origin (minimum coordinate point of AABB).
   STEP		: Grid consists of adjacent cubic cells with side length equal to STEP.
   X, Y, Z	: Number of cells along axes x, y anx z, respectively.
   LINEAR_MAP   : Array mapping linearized (3D->1D, by function LINEAR_ID) indices
                  to linear indices (on ANCHOR_MAP), thus implementing the locality
		  preserving 3D->1D map.
   ANCHORS      : Array containing the index (anchor) to the first particle of the
		  corresponding cell in LINEAR_GRID. Anchors correspond to cells
		  in ascending linear order, so the relation A[i-1] <= A[i] holds
		  for all elements of this array. When A[i-1]=A[i], the cell
		  corresponding to A[i-1] contains no particles.
   LINEAR_GRID  : Array containing pointers to particles
*/
struct lp_grid {
    unsigned int* linear_map;
    anchor* anchors;
    particle** linear_grid;
    btVector3 origin;
    unsigned int x, y, z, cell_count, particle_count;
    float step;
};

int linear_id (lp_grid lpg, int i, int j, int k);

index3 inv_linear_id (lp_grid lpg, int linear_id);

lp_grid make_lp_grid (btVector3 origin, float step, int x, int y, int z, std::vector<particle> particles);

int populate_lp_grid (lp_grid &lpg, std::vector<particle> particles);

/* std::vector<particle*> lp_grid_ref (lp_grid g, int x, int y, int z); */

#endif
