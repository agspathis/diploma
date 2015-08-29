#include <fstream>

#include "output.h"

int vtk_export_particles (char* filename, fluid fluid)
{
    int size = fluid.particle_count;
    std::ofstream vtk;
    vtk.open(filename);

    // standard header
    vtk << "# vtk DataFile Version 3.1\n";
    vtk << "Simulation frame\n";
    vtk << "ASCII\n";
    vtk << "DATASET UNSTRUCTURED_GRID\n";

    // point location
    vtk << "POINTS " << size << " FLOAT\n";
    for (long i=0; i<size; i++) {
	btVector3 position = particle_position(fluid.particles+i);
	vtk << position.getX() << " " << position.getY() << " " << position.getZ() << "\n";
    }
    vtk << "\n";

    // dummy cell data
    vtk << "CELLS 1 " << size+1 << "\n";
    vtk << size;
    for (int i=0; i < size; i++) vtk << " " << i;
    vtk << "\n\n";

    vtk << "CELL_TYPES 1\n2\n\n";

    vtk.close();
    return 0;
}
