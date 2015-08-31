#include <iostream>

#include "utilities.h"
#include "output.h"
#include "terrain.h"
#include "fluid.h"
#include "lp_grid.h"
#include "sph.h"

// Constants
#define STEPS 500
#define PARTICLES 1000

// Global parameters
char obj_filename[] = "/home/agspathis/diplom/models/obj/box.obj";
aabb fluid_aabb = { btVector3(-50, -50, -50), btVector3(50, 50, 50) };
char output_dir[] = "/home/agspathis/diplom/frames";

int main (void)
{
    // dynamics world construction
    btDefaultCollisionConfiguration* collision_configuration = new btDefaultCollisionConfiguration();
    btCollisionDispatcher* dispatcher = new btCollisionDispatcher(collision_configuration);
    btBroadphaseInterface* broadphase = new btDbvtBroadphase();
    btSequentialImpulseConstraintSolver* solver = new btSequentialImpulseConstraintSolver;
    btDiscreteDynamicsWorld* dynamics_world =
	new btDiscreteDynamicsWorld(dispatcher, broadphase, solver, collision_configuration);
    dynamics_world->setGravity(btVector3(0, -9.81, 0));

    // terrain, fluid, lp grid construction
    terrain terrain = make_terrain_obj(obj_filename);
    dynamics_world->addRigidBody(terrain.rigid_body);

    fluid fluid = make_fluid(fluid_aabb, PARTICLES);
    for (long pi=0; pi<fluid.particle_count; pi++)
	dynamics_world->addRigidBody(fluid.particles[pi].rigid_body);

    btVector3 origin = btVector3(0, 0, 0);
    lp_grid lpg = make_lp_grid (terrain.taabb, fluid);

    // change to and clear output directory
    chdir(output_dir);
    system("exec rm *");

    // simulation
    for (long step=0; step<STEPS; step++) {
	dynamics_world->stepSimulation(1/60.f, 10, 1/100.f);
	printf("Frame %lu\n", step);
	// apply_sph(lpg, fluid);
	vtk_export_particles(output_dir, fluid, step);
	update_lp_grid(lpg);
    }

    // cleanup
    dynamics_world->removeRigidBody(terrain.rigid_body);
    delete_terrain(terrain);

    for (long pi=0; pi<fluid.particle_count; pi++)
	dynamics_world->removeRigidBody(fluid.particles[pi].rigid_body);
    delete_fluid(fluid);

    delete_lp_grid(lpg);

    delete dynamics_world;
    delete solver;
    delete collision_configuration;
    delete dispatcher;
    delete broadphase;

    return 0;
}
