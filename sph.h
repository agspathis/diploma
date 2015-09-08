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
    fluid* f;
    lp_grid* lpg;
};

void apply_sph(fluid_sim* fsim);

#endif
