#include <fstream>

#include "utilities.h"

// Random float with upper bound b
float rand_fb (float b)
{
    return static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/b));
}

btVector3 particle_position (particle p)
{
    btTransform tf;
    p.rigid_body->getMotionState()->getWorldTransform(tf);
    return btVector3(tf.getOrigin().getX(), tf.getOrigin().getY(), tf.getOrigin().getZ());
}


// Export vector of particles to .vtk file
int vtk_export (const char* filename, std::vector<particle> particles)
{
    int size = particles.size();
    std::ofstream vtk;
    vtk.open(filename);

    // standard header
    vtk << "# vtk DataFile Version 3.1\n";
    vtk << "Simulation frame\n";
    vtk << "ASCII\n";
    vtk << "DATASET UNSTRUCTURED_GRID\n";

    // point location
    vtk << "POINTS " << size << " FLOAT\n";
    btTransform tf;
    for (int i=0; i < size; i++) {
	particles[i].rigid_body->getMotionState()->getWorldTransform(tf);
	vtk << tf.getOrigin().getX() << " "
	    << tf.getOrigin().getY() << " "
	    << tf.getOrigin().getZ() << "\n";
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

int print_long_array(long* array, long count)
{
    for (long index=0; index<count; index++) {
	printf("%d ", array[index]);
	if (index % 16) continue;
	printf("\n");
    }
    printf("\n");
    return 0;
}
