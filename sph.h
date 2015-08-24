#ifndef SPH_H
#define SPH_H

#include <vector>
#include <btBulletDynamicsCommon.h>

#include "fluid.h"
#include "terrain.h"
#include "lp_grid.h"

struct interaction {
    particle* p0;
    particle* p1;
    btScalar distance;
    btVector3 direction;
    btScalar force;
};

float kernel(float r, float h);

int apply_sph_forces(lp_grid lpg, btScalar particle_mass);

#endif
