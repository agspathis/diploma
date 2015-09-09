#ifndef SPH_H
#define SPH_H

#include "lp_grid.h"

struct interaction {
    particle* p0;
    particle* p1;
    btScalar distance;
    btVector3 direction;
};

struct fluid_sim {
    fluid f;
    lp_grid lpg;
};

void apply_sph(fluid_sim fsim);

void adjust_fluid(fluid* fluid, lp_grid lpg, aabb fluid_aabb, aabb terrain_aabb);

#endif
