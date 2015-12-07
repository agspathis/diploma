#include <sys/time.h>

#include "utilities.h"
#include "output.h"
#include "terrain.h"
#include "fluid.h"
#include "lp_grid.h"
#include "sph.h"

// Constants and global parameters
#define FRAMES 100
#define SAMPLES 50
#define FRAME_DT 0.05
#define PARTICLES 5000
#define TERRAIN_SCALING_FACTOR 0.04
const char* output_dir = "../frames";
const char* coast_filename = "../models/city_0.obj";
aabb sea_aabb = { btVector3(0, 2, 0), btVector3(6, 6, 84) };

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
    terrain coast = make_terrain_obj(coast_filename, TERRAIN_SCALING_FACTOR);
    terrain boundary = terrain_boundary(coast);
    dynamics_world->addRigidBody(coast.rigid_body);
    dynamics_world->addRigidBody(boundary.rigid_body);
    std::vector<terrain_impulse> terrain_impulses;

    fluid sea = make_fluid(sea_aabb, PARTICLES, SAMPLES);
    for (long pi=0; pi<sea.particle_count; pi++)
	dynamics_world->addRigidBody(sea.particles[pi].rigid_body);

    lp_grid lpg = make_lp_grid(coast.taabb, sea);

    adjust_fluid(&sea, lpg, sea_aabb, coast.taabb);
    fluid_sim fsim;
    fsim.f = sea;
    fsim.lpg = lpg;
    fsim.tis = terrain_impulses;
    dynamics_world->setInternalTickCallback(tick_callback, (void*) &fsim);

    // change to and clear output directory
    chdir(output_dir);
    system("exec rm *");

    // export scaled/docked terrain
    terrain_to_obj(output_dir, coast);

    // simulation
    printf("SIMULATION START\n");
    printf("%d physics ticks per frame\n", (int) ceil(FRAME_DT/sea.dt));
    struct timeval start, end, diff;
    for (int frame=0; frame<FRAMES; frame++) {
	gettimeofday(&start, NULL);

	// sim step and data export
	dynamics_world->stepSimulation(FRAME_DT, ceil(FRAME_DT/sea.dt), sea.dt);
	compute_cf(lpg);
	particles_to_vtk(output_dir, sea, frame);
	color_field_to_vtk(output_dir, lpg, frame);
	terrain_impulses_to_vtk(output_dir, fsim.tis, frame);

	// clear impulse data accumulated over last step
	fsim.tis.clear();

	// frame log
	gettimeofday(&end, NULL);
	timersub(&end, &start, &diff);
	printf("Frame %03d/%03d, %ld.%06ld seconds\n",
	       frame, FRAMES-1, (int) diff.tv_sec, (long) diff.tv_usec);
    }

    // cleanup
    dynamics_world->removeRigidBody(coast.rigid_body);
    dynamics_world->removeRigidBody(boundary.rigid_body);
    delete_terrain(coast);
    delete_terrain(boundary);

    for (long pi=0; pi<sea.particle_count; pi++)
	dynamics_world->removeRigidBody(sea.particles[pi].rigid_body);
    delete_fluid(sea);

    delete_lp_grid(lpg);

    delete dynamics_world;
    delete solver;
    delete collision_configuration;
    delete dispatcher;
    delete broadphase;

    return 0;
}
