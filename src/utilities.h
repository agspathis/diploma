#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <btBulletDynamicsCommon.h>

#include "fluid.h"

// Random float with upper bound b
float rand_fb(float b);

int print_aabb(aabb aabb);

int print_int_array(long* array, long count);

#endif
