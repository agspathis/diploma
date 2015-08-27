#include <iostream>

#include "utilities.h"
#include "output.h"
#include "terrain.h"
#include "fluid.h"
#include "lp_grid.h"
#include "sph.h"

// Constants
#define STEPS 20

// Global parameters
char obj_filename[] = "/home/agspathis/diplom/models/obj/box.obj";

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

    // terrain construction
    terrain terrain = import_obj(obj_filename);
    dynamics_world->addRigidBody(terrain.rigid_body);

    // fluid construction
    aabb fluid_aabb;
    fluid_aabb.min = btVector3(-50, -50, -50);
    fluid_aabb.max = btVector3(50, 50, 50);
    btScalar particle_mass = 1.0;
    btScalar particle_radius = 2.0;
    float smoothing_length = 4 * particle_radius;
    std::vector<particle*> particles =
	fluid_fill(fluid_aabb, particle_mass, particle_radius, dynamics_world);

    // grid construction
    btVector3 origin = btVector3(0, 0, 0);
    lp_grid lpg = make_lp_grid (terrain.terrain_aabb, smoothing_length, particles);

    // simulation
    for (int i=0; i<STEPS; i++) {
	// stepping
	dynamics_world->stepSimulation(1/60.f, 10, 1/100.f);
	printf("Frame %d\n", i);
	// sph
	apply_sph(lpg, particle_mass);
	// export to vtk
	std::string filepath = "/home/agspathis/diplom/frames/frame"+std::to_string(i)+".vtk";
	vtk_export((const char*) filepath.c_str(), particles);
    }

    // cleanup
    for (int i=0; i<particles.size(); i++) {
    	dynamics_world->removeRigidBody(particles[i]->rigid_body);
    	delete particles[i]->rigid_body->getMotionState();
    	delete particles[i]->rigid_body;
    }

    dynamics_world->removeRigidBody(terrain.rigid_body);
    delete_terrain(terrain);

    delete dynamics_world;
    delete solver;
    delete collision_configuration;
    delete dispatcher;
    delete broadphase;

    return 0;
}
