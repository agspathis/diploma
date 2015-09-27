#include <sys/time.h>

#include "utilities.h"
#include "output.h"
#include "terrain.h"
#include "fluid.h"
#include "lp_grid.h"
#include "sph.h"

// Constants
#define FRAMES 100
#define SAMPLES 50
#define FRAME_DT 0.05
#define PARTICLES 4000
#define TERRAIN_SCALING_FACTOR 0.04

// Collision groups
enum collisiontypes { TCOL = 1, PCOL = 2 };

// Global parameters
const char* output_dir = "/home/agspathis/diplom/frames";
const char* terrain_filename = "/home/agspathis/diplom/models/obj/the-city.obj";
aabb fluid_aabb = { btVector3(0, 2, 0), btVector3(6, 5, 84) };

void tick_callback(btDynamicsWorld* dynamics_world, btScalar timeStep) {
    fluid_sim* fsimp = (fluid_sim*) dynamics_world->getWorldUserInfo();
    collect_terrain_impulses(dynamics_world, fsimp->tis);
    update_lp_grid(fsimp->lpg);
    apply_sph(fsimp);
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
    terrain terrain = make_terrain_obj(terrain_filename, TERRAIN_SCALING_FACTOR);
    // dynamics_world->addRigidBody(terrain.rigid_body, TCOL, PCOL);
    dynamics_world->addRigidBody(terrain.rigid_body);
    std::vector<terrain_impulse> terrain_impulses;

    fluid fluid = make_fluid(fluid_aabb, PARTICLES, SAMPLES);
    for (long pi=0; pi<fluid.particle_count; pi++)
	// dynamics_world->addRigidBody(fluid.particles[pi].rigid_body, PCOL, TCOL);
	dynamics_world->addRigidBody(fluid.particles[pi].rigid_body);

    lp_grid lp_grid = make_lp_grid(terrain.taabb, fluid);

    adjust_fluid(&fluid, lp_grid, fluid_aabb, terrain.taabb);
    fluid_sim fluid_sim;
    fluid_sim.f = fluid;
    fluid_sim.lpg = lp_grid;
    fluid_sim.tis = terrain_impulses;
    dynamics_world->setInternalTickCallback(tick_callback, (void*) &fluid_sim);

    // change to and clear output directory
    chdir(output_dir);
    system("exec rm *");

    // export docked/scaled terrain
    terrain_to_obj(output_dir, terrain);

    // simulation
    printf("SIMULATION START\n");
    struct timeval start, end, diff;
    for (int frame=0; frame<FRAMES; frame++) {
	gettimeofday(&start, NULL);

	// sim step and data export
	dynamics_world->stepSimulation(FRAME_DT, ceil(FRAME_DT/fluid.dt), fluid.dt);
	compute_cf(lp_grid);
	particles_to_vtk(output_dir, fluid, frame);
	color_field_to_vtk(output_dir, lp_grid, frame);
	terrain_impulses_to_vtk(output_dir, fluid_sim.tis, frame);

	// clear impulse data accumulated over last step
	fluid_sim.tis.clear();

	// frame log
	gettimeofday(&end, NULL);
	timersub(&end, &start, &diff);
	printf("Frame %03d/%03d, %ld.%06ld seconds\n",
	       frame, FRAMES-1, (int) diff.tv_sec, (long) diff.tv_usec);
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
