#ifndef OUTPUT_H
#define OUTPUT_H

#include <unistd.h>

#include "fluid.h"

void vtk_export_particles(const char* output_dir, fluid fluid, long frame);

#endif
