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
    terrain t;
    fluid f;
    lp_grid lpg;
    std::vector<terrain_impulse> tis;
};

void apply_sph(fluid_sim* fsimp);

void adjust_fluid(fluid* fluid, lp_grid lpg, aabb fluid_aabb, aabb terrain_aabb);

void compute_cf(lp_grid lpg);

void accumulate_if(lp_grid lpg, std::vector<terrain_impulse> tis);

#endif
