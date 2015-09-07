#include "output.h"

void vtk_export_particles(char* output_dir, fluid fluid, long frame)
{
    chdir(output_dir);
    char filename[sizeof "particles_000000.vtk"];
    sprintf(filename, "particles_%06d.vtk", frame);
    FILE* vtk = fopen(filename, "w");
    long size = fluid.particle_count;

    // standard header
    fprintf(vtk,"# vtk DataFile Version 3.1\n");
    fprintf(vtk,"Simulation particles\n");
    fprintf(vtk,"ASCII\n");
    fprintf(vtk,"DATASET UNSTRUCTURED_GRID\n");

    // point location
    fprintf(vtk, "POINTS %lu FLOAT\n", size);
    for (particle* pp=fluid.particles; pp<fluid.particles+size; pp++) {
	btVector3 position = particle_position(pp);
	fprintf(vtk, "%f %f %f\n", position.getX(), position.getY(), position.getZ());
    }
    fprintf(vtk, "\n");

    // dummy cell data
    fprintf(vtk, "CELLS 1 %lu\n", size+1);
    fprintf(vtk, "%lu", size);
    for (long i=0; i<size; i++) fprintf(vtk, " %lu", i);
    fprintf(vtk, "\n\n");

    fprintf(vtk, "CELL_TYPES 1\n2\n\n");

    // point data
    fprintf(vtk, "POINT_DATA %lu\n", size);
    fprintf(vtk, "VECTORS pressure DOUBLE\n");
    // fprintf(vtk, "LOOKUP_TABLE default\n");
    for (particle* pp=fluid.particles; pp<fluid.particles+size; pp++) {
	btVector3 pforce = pp->pforce;
	fprintf(vtk, "%f %f %f\n", pforce.getX(), pforce.getY(), pforce.getZ());
    }
    fprintf(vtk, "\n");

    fprintf(vtk, "VECTORS viscosity DOUBLE\n");
    // fprintf(vtk, "LOOKUP_TABLE default\n");
    for (particle* pp=fluid.particles; pp<fluid.particles+size; pp++) {
	btVector3 vforce = pp->vforce;
	fprintf(vtk, "%f %f %f\n", vforce.getX(), vforce.getY(), vforce.getZ());
    }
    fprintf(vtk, "\n");

    fclose(vtk);
}
