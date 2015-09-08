#include <iostream>

#include "utilities.h"
#include "output.h"
#include "terrain.h"
#include "fluid.h"
#include "lp_grid.h"
#include "sph.h"

// Constants
#define STEPS 300
#define PARTICLES 1000

// Global parameters
char output_dir[] = "/home/agspathis/diplom/frames";
char obj_filename[] = "/home/agspathis/diplom/models/obj/box-small.obj";
aabb fluid_aabb = { btVector3(2, 2, 2), btVector3(8, 8, 8) };

void tick_callback(btDynamicsWorld* dynamics_world, btScalar timeStep) {
    fluid_sim* f_simp = (fluid_sim*) dynamics_world->getWorldUserInfo();
    update_lp_grid(*(f_simp->lpg));
    // apply_sph(f_simp);
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
    dynamics_world->setGravity(btVector3(0, -9.81, 0));

    // terrain, fluid, lp_grid and fluid_sim construction
    terrain terrain = make_terrain_obj(obj_filename);
    dynamics_world->addRigidBody(terrain.rigid_body);

    fluid fluid = make_fluid(fluid_aabb, PARTICLES);
    for (long pi=0; pi<fluid.particle_count; pi++)
	dynamics_world->addRigidBody(fluid.particles[pi].rigid_body);

    lp_grid lpg = make_lp_grid (terrain.taabb, fluid);
    
    fluid_sim fluid_sim;
    fluid_sim.f = &fluid;
    fluid_sim.lpg = &lpg;
    dynamics_world->setInternalTickCallback(tick_callback, (void*) &fluid_sim);

    // change to and clear output directory
    chdir(output_dir);
    system("exec rm *");

    // simulation
    for (long step=0; step<STEPS; step++) {
	dynamics_world->stepSimulation(0.1, 100, 0.001);
	printf("Frame %lu\n", step);
	vtk_export_particles(output_dir, fluid, step);
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
