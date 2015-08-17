#ifndef SPH_H
#define SPH_H

#include <vector>
#include <btBulletDynamicsCommon.h>

#include "fluid.h"
#include "terrain.h"

struct interaction {
    int i;
    int j;
    btScalar distance;
    btVector3 direction;
    btScalar force;
};

float kernel(float r, float h);

int apply_sph_forces(std::vector<particle> particles,
		     float smoothing_length,
		     btScalar particle_mass);

#endif
