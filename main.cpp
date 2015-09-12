#include <iostream>

#include "utilities.h"
#include "output.h"
#include "terrain.h"
#include "fluid.h"
#include "lp_grid.h"
#include "sph.h"

// Constants
#define STEPS 1
#define FRAME_DT 0.05
#define PARTICLES 10000
#define TERRAIN_SCALING_FACTOR 0.04

// Global parameters
const char* output_dir = "/home/agspathis/diplom/frames";
const char* obj_filename = "/home/agspathis/diplom/models/obj/the-city.obj";
aabb fluid_aabb = { btVector3(0, 0, 10), btVector3(5, 10, 80) };

void tick_callback(btDynamicsWorld* dynamics_world, btScalar timeStep) {
    fluid_sim fsim = *((fluid_sim*) dynamics_world->getWorldUserInfo());
    update_lp_grid(fsim.lpg);
    apply_sph(fsim);
}

int main (void)
{
    // dynamics world construction
    btDefaultCollisionConfiguration* collision_configuration =
	new btDefaultCollisionConfiguration();
    btCollisionDispatcher* dispatcher =
	new btCollisionDispatcher(collision_configuration);
    btBroadphaseInterface* broadphase =
	new btDbvtBroadphase();
    btSequentialImpulseConstraintSolver* solver =
	new btSequentialImpulseConstraintSolver;
    btDiscreteDynamicsWorld* dynamics_world =
	new btDiscreteDynamicsWorld(dispatcher, broadphase,
				    solver, collision_configuration);
    dynamics_world->setGravity(btVector3(0, -G, 0));

    // terrain, fluid, lp_grid and fluid_sim construction
    terrain terrain = make_terrain_obj(obj_filename, TERRAIN_SCALING_FACTOR);
    dynamics_world->addRigidBody(terrain.rigid_body);

    fluid fluid = make_fluid(fluid_aabb, PARTICLES);
    for (long pi=0; pi<fluid.particle_count; pi++)
	dynamics_world->addRigidBody(fluid.particles[pi].rigid_body);

    lp_grid lp_grid = make_lp_grid(terrain.taabb, fluid);

    adjust_fluid(&fluid, lp_grid, fluid_aabb, terrain.taabb);
    fluid_sim fluid_sim;
    fluid_sim.f = fluid;
    fluid_sim.lpg = lp_grid;
    dynamics_world->setInternalTickCallback(tick_callback, (void*) &fluid_sim);

    // change to and clear output directory
    chdir(output_dir);
    system("exec rm *");

    // export docked/scaled terrain
    obj_export_terrain(output_dir, terrain);

    // simulation
    for (int step=0; step<STEPS; step++) {
	dynamics_world->stepSimulation(FRAME_DT, ceil(FRAME_DT/fluid.dt), fluid.dt);
	printf("Frame %d / %d\n", step, STEPS-1);
	vtk_export_particles(output_dir, fluid, step);
    }

    // cleanup
    dynamics_world->removeRigidBody(terrain.rigid_body);
    delete_terrain(terrain);

    for (long pi=0; pi<fluid.particle_count; pi++)
	dynamics_world->removeRigidBody(fluid.particles[pi].rigid_body);
    delete_fluid(fluid);

    delete_lp_grid(lp_grid);

    delete dynamics_world;
    delete solver;
    delete collision_configuration;
    delete dispatcher;
    delete broadphase;

    return 0;
}
