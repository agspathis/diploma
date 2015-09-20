#ifndef LP_GRID_H
#define LP_GRID_H

#include "fluid.h"

#define DEFAULT_CF_SDF 4

typedef particle** anchor;

/*
  Structure LP_GRID: Data structure to store pointers to the simulation particles,
  preserving locality between simulation and memory space.
  Member documentation:
  ORIGIN	  : Coordinates of grid origin (minimum-coordinate point of AABB).
  STEP		  : Grid consists of adjacent cubic cells with side length equal to STEP.
  X, Y, Z	  : Number of cells along axes x, y and z, respectively.
  XSS, YSS, ZSS   : Segment sizes along each axis for load distribution among threads
  CELL_COUNT      : Total number of cells (X*Y*Z).
  PARTICLE_COUNT  : Total number of particles.
  MAP		  : Array containing pointers to ANCHORS. Each pointer is stored
		    at the index computed by the function LINEARIZE from the 3D
		    address of the respective cell in the grid. The pointer value
		    is determined by the cell's position in the spatial ordering
		    along a space-filling curve, thus implementing the locality
		    preserving 3D->1D mapping. LINEARIZE returns CELL_COUNT if
		    given address out of grid bounds (X, Y, Z), which maps to the
		    last cell on ANCHORS.
  ANCHORS         : Array containing the anchor for each cell in spatial order,
		    i.e. pointer to the first of cell particles in PARTICLES, so
		    the relation A[i-1] <= A[i] holds for all anchors of this
		    array. When A[i-1]=A[i], the cell corresponding to A[i-1]
		    contains no particles. This array has CELL_COUNT+1 elements,
		    where the last anchor refers to a virtual cell containing
		    particles which are out of grid bounds (also used for loop
		    termination conditions).
  PARTICLES       : Array containing pointers to actual simulation particles. It
		    has PARTICLE_COUNT+1 elements. The extra last slot is used
		    for loop termination conditions and valid PARTICLE_RANGE for
		    the last cell (containing out-of-grid particles).
  CF_SDF          : Color field subdivision factor is the number of subdivisions
                    inside each cell along each axis (same for x, y, z), for
		    addressing the higher resolution color field subgrid aligned
		    with the master grid.
  COLOR_FIELD     : Array containing the samples of the scalar color field
                    (color = dimensionless density), which are computed at the
		    minimum vertex of each cell (cube). The isosurface at an
		    appropriate value is a good estimation of fluid surface.
*/
struct lp_grid {
    btVector3 origin;
    float step;
    long x, y, z;
    long xss, yss, zss;
    long cell_count, particle_count;
    anchor** map;
    anchor* anchors;
    particle** particles;
    int cf_sdf;
    float* color_field;
};

lp_grid make_lp_grid(aabb domain, fluid fluid);

void update_lp_grid(lp_grid lpg);

void delete_lp_grid(lp_grid lpg);

/*
  Structure CELL represents a cell as the range between a pair of pointers to
  LPG.PARTICLES, from START (included -- anchor of this cell) to END (not
  included -- anchor of next cell in spatial sort order).
*/
struct cell {
    anchor start;
    anchor end;
};

long linearize(lp_grid lpg, long i, long j, long k, int f);

cell get_cell(lp_grid lpg, long i, long j, long k);

#endif
