#ifndef OUTPUT_H
#define OUTPUT_H

#include "fluid.h"

// Export vector of particles to .vtk file
int vtk_export (const char* filename, std::vector<particle*> particles);

#endif
