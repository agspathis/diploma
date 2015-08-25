#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include <btBulletDynamicsCommon.h>

#include "terrain.h"

struct particle {
    int samples;
    float density;
    float pressure;
    btRigidBody* rigid_body;
};

std::vector<particle*> fluid_fill(aabb aabb, btScalar pMass, btScalar pRadius,
				  btDiscreteDynamicsWorld* dynamics_world);

#endif
