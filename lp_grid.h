#ifndef LP_GRID_H
#define LP_GRID_H

#include <btBulletDynamicsCommon.h>
#include <vector>

#include "terrain.h"
#include "fluid.h"

typedef particle** anchor;

/*
  Structure LP_GRID: Data structure to store pointers to the simulation particles,
  preserving locality between simulation and memory space.
  Member documentation:
  ORIGIN	  : Coordinates of grid origin (minimum-coordinate point of AABB).
  STEP		  : Grid consists of adjacent cubic cells with side length equal to STEP.
  X, Y, Z	  : Number of cells along axes x, y and z, respectively.
  XSL, YSL, ZSL   : Segment lengths along each axis for load distribution among threads
  CELL_COUNT      : Total number of cells (X*Y*Z).
  PARTICLE_COUNT  : Total number of particles.
  MAP		  : Array map of linearized addresses (3D->1D by function LINEARIZE_ADDRESS)
                    to linear indices (on ANCHOR_MAP), thus implementing the locality 
		    preserving 3D->1D mapping.
  ANCHORS         : Array containing the cell ANCHORs, i.e. index of the first of cell
		    particles in PARTICLES. ANCHORs correspond to cells in ascending linear
		    order, so the relation A[i-1] <= A[i] holds for all elements of this 
		    array. When A[i-1]=A[i], the cell corresponding to A[i-1] contains
		    no particles. This array has CELL_COUNT+1 elements, where the last
		    anchor refers to a spare cell containing particles which are out of
		    grid bounds (also serves for loop termination conditions).
  PARTICLES       : Array containing pointers to actual simulation particles. It has
                    PARTICLE_COUNT+1 elements. The extra slot is used for spare storage
		    to gradually initialize the array, for loop termination conditions
		    and valid PARTICLE_RANGE for the last cell.
  
*/
struct lp_grid {
    btVector3 origin;
    float step;
    long x, y, z;
    long xsl, ysl, zsl;
    long cell_count, particle_count;
    long* map;
    long* anchors;
    particle** particles;
};

long linearize_address (lp_grid lpg, long i, long j, long k);

lp_grid make_lp_grid (aabb domain, float step, std::vector<particle> particles);

int update_lp_grid (lp_grid lpg);

/*
  Structure CELL represents a cell as a range of consecutive pointers to its
  particles, from START (included) to END (not included -- first particle of
  next cell in spatial sort order).
*/
struct cell {
    anchor start;
    anchor end;
};

cell get_cell(lp_grid lpg, long i, long j, long k);

#endif
