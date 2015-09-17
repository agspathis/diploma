#ifndef OUTPUT_H
#define OUTPUT_H

#include <unistd.h>

#include "fluid.h"
#include "lp_grid.h"

void particles_to_vtk(const char* output_dir, fluid fluid, long frame);

void terrain_to_obj(const char* output_dir, terrain t);

void color_field_to_vtk(const char* output_dir, lp_grid lpg, long frame);

#endif
