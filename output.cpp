#include "output.h"

void vtk_export_particles (char* output_dir, fluid fluid, long frame)
{
    chdir(output_dir);
    char filename[sizeof "particles_000000.vtk"];
    sprintf(filename, "particles_%06d.vtk", frame);
    FILE* vtk = fopen(filename, "w");
    int size = fluid.particle_count;

    // standard header
    fprintf(vtk,"# vtk DataFile Version 3.1\n");
    fprintf(vtk,"Simulation frame\n");
    fprintf(vtk,"ASCII\n");
    fprintf(vtk,"DATASET UNSTRUCTURED_GRID\n");

    // point location
    fprintf(vtk, "POINTS %lu FLOAT\n", size);
    for (long i=0; i<size; i++) {
	btVector3 position = particle_position(fluid.particles+i);
	fprintf(vtk, "%f %f %f\n", position.getX(), position.getY(), position.getZ());
    }
    fprintf(vtk, "\n");

    // dummy cell data
    fprintf(vtk, "CELLS 1 %lu\n", size+1);
    fprintf(vtk, "%lu", size);
    for (long i=0; i<size; i++) fprintf(vtk, " %lu", i);
    fprintf(vtk, "\n\n");

    fprintf(vtk, "CELL_TYPES 1\n2\n\n");

    fclose(vtk);
}
