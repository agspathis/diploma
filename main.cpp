#include <string>
#include <iostream>

#include "utilities.h"
#include "terrain.h"
#include "fluid.h"
#include "lp_grid.h"
#include "sph.h"

// Constants
#define STEPS 20

// Global parameters
const char* obj_filename = "/home/agspathis/diplom/models/obj/box.obj";

int main (void)
{
    // Construct dynamics world
    btDefaultCollisionConfiguration* collision_configuration = new btDefaultCollisionConfiguration();
    btCollisionDispatcher* dispatcher = new btCollisionDispatcher(collision_configuration);
    btBroadphaseInterface* broadphase = new btDbvtBroadphase();
    btSequentialImpulseConstraintSolver* solver = new btSequentialImpulseConstraintSolver;
    btDiscreteDynamicsWorld* dynamics_world =
	new btDiscreteDynamicsWorld(dispatcher, broadphase, solver, collision_configuration);

    // Set dynamics world gravity
    dynamics_world->setGravity(btVector3(0, -9.81, 0));

    // Terrain construction
    aabb terrain_aabb;		// terrain AABB
    btTriangleMesh* triangle_mesh = import_obj(obj_filename, terrain_aabb);
    btCollisionShape* t_shape = new btBvhTriangleMeshShape(triangle_mesh,true);
    btDefaultMotionState* t_motion_state =
	new btDefaultMotionState(btTransform(btQuaternion(0, 0, 0, 1), btVector3(0, 0, 0)));
    btRigidBody::btRigidBodyConstructionInfo t_ci(0, t_motion_state, t_shape, btVector3(0, 0, 0));
    t_ci.m_restitution = 1.0;
    btRigidBody* terrain = new btRigidBody(t_ci);
    dynamics_world->addRigidBody(terrain);

    printf("%f %f %f\n", terrain_aabb.min.getX(), terrain_aabb.min.getY(), terrain_aabb.min.getZ());

    // Particle construction
    aabb fluid_aabb;
    fluid_aabb.min = btVector3(-50, -50, -50);
    fluid_aabb.max = btVector3(50, 50, 50);
    btScalar particle_mass = 1.0;
    btScalar particle_radius = 4.0;
    float smoothing_length = 4 * particle_radius;
    std::vector<particle> particles = fluid_fill(fluid_aabb, particle_mass,
						 particle_radius, dynamics_world);
    printf("Particle count: %d\n", particles.size());

    // Grid construction
    btVector3 origin = btVector3(0, 0, 0);
    make_lp_grid (terrain_aabb, 10.0, particles);

    // Simulation
    for (int i=0; i<STEPS; i++) {
	// stepping
	dynamics_world->stepSimulation(1/60.f, 10, 1/100.f);
	printf("%d\n", i);
	// sph forces
	apply_sph_forces(particles, smoothing_length, particle_mass);
	// export to vtk
	std::string filepath = "frames/frame"+std::to_string(i)+".vtk";
	vtk_export((const char*) filepath.c_str(), particles);
    }

    // cleanup
    for (int i=0; i<particles.size(); i++) {
    	dynamics_world->removeRigidBody(particles[i].rigid_body);
    	delete particles[i].rigid_body->getMotionState();
    	delete particles[i].rigid_body;
    }

    dynamics_world->removeRigidBody(terrain);
    delete t_shape;
    delete triangle_mesh;
    delete terrain->getMotionState();
    delete terrain;

    delete dynamics_world;
    delete solver;
    delete collision_configuration;
    delete dispatcher;
    delete broadphase;

    return 0;
}
