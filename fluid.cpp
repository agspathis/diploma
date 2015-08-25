#include "fluid.h"

std::vector<particle*> fluid_fill(aabb aabb, btScalar particle_mass, btScalar particle_radius,
				  btDiscreteDynamicsWorld* dynamics_world)
{
    btCollisionShape* rb_shape;
    btVector3 rb_inertia(0, 0, 0);
    rb_shape = new btSphereShape(particle_radius);
    rb_shape->calculateLocalInertia(particle_mass, rb_inertia);
    btDefaultMotionState* rb_motion_state;
    std::vector<particle*> particles;

    float x_min = aabb.min.getX();
    float y_min = aabb.min.getY();
    float z_min = aabb.min.getZ();
    float x_max = aabb.max.getX();
    float y_max = aabb.max.getY();
    float z_max = aabb.max.getZ();

    float dx = x_max - x_min;
    float dy = y_max - y_min;
    float dz = z_max - z_min;

    float x, y, z;

    for (int k=0; (z = 1.633 * k * particle_radius) < dz; k++) {
	for (int j=0; (y = (j + (k%2)/3.0) * 1.732 * particle_radius) < dy; j++) {
	    for (int i=0; (x = (2*i + ((j+k)%2)) * particle_radius) < dx; i++) {
		rb_motion_state =
		    new btDefaultMotionState (btTransform(btQuaternion(0, 0, 0, 1),
							  btVector3(x_min+x, y_min+y, z_min+z)));
		btRigidBody::btRigidBodyConstructionInfo
		    rb_ci(particle_mass, rb_motion_state, rb_shape, rb_inertia);
		rb_ci.m_restitution = 1.0;
		particle* par = (particle*) malloc(sizeof(particle*));
		par->rigid_body = new btRigidBody(rb_ci);
		par->rigid_body->setLinearVelocity(btVector3(0, 0, 0));
		par->rigid_body->setAngularFactor(btVector3(0, 0, 0));
		dynamics_world->addRigidBody(par->rigid_body);
		particles.push_back(par);
	    }
	}
    }
    return particles;
}
