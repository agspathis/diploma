#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <btBulletDynamicsCommon.h>

#include "fluid.h"

// Random float with upper bound b
float rand_fb(float b);

btVector3 particle_position (particle p);

// Export vector of particles to .vtk file
int vtk_export(const char* filename, std::vector<particle> particles);

#endif
