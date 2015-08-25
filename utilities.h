#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <btBulletDynamicsCommon.h>

#include "fluid.h"

// Random float with upper bound b
float rand_fb(float b);

btVector3 particle_position (particle* pp);

// Export vector of particles to .vtk file
int vtk_export (const char* filename, std::vector<particle*> particles);

int print_aabb(aabb aabb);

int print_long_array(long* array, long count);

#endif
