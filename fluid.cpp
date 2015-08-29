#include "fluid.h"

fluid make_fluid(aabb aabb, long desired_particle_count)
{
    fluid f;
    // AABB volume extraction and loop limit setting
    float x_min = aabb.min.getX();
    float y_min = aabb.min.getY();
    float z_min = aabb.min.getZ();
    float x_max = aabb.max.getX();
    float y_max = aabb.max.getY();
    float z_max = aabb.max.getZ();
    float dx = x_max - x_min;
    float dy = y_max - y_min;
    float dz = z_max - z_min;
    float volume = dx*dy*dz;
    // F.PARTICLE_RADIUS computed based on the fraction of space that equal
    // spheres in closest packing (HPC lattice) occupy is pi/(3*sqrt(2))
    f.particle_radius = cbrt(volume/(5.6568542*desired_particle_count));

    // values for water
    f.density = 1000;
    f.dynamic_viscosity = 0.001;
    f.tait_b = 10000;

    // particle counting
    float x, y, z;
    long i, j, k, pi;
    f.particle_count = 0;
    for (k=0; (z = 1.633 * k * f.particle_radius) < dz; k++)
	for (j=0; (y = (j + (k%2)/3.0) * 1.732 * f.particle_radius) < dy; j++)
	    for (i=0; (x = (2*i + ((j+k)%2)) * f.particle_radius) < dx; i++)
		f.particle_count++;
    f.particle_mass = (volume * f.density)/f.particle_count;
    
    // F.SMOOTHING_RADIUS for about 50 smoothing samples
    f.smoothing_radius = 4.2 * f.particle_radius;
    f.sample_density = f.density / f.particle_mass;
    f.particles = (particle*) malloc(f.particle_count * sizeof(particle));

    // particle construction
    pi = 0;			// pi = particle index
    btVector3 fp_inertia(0, 0, 0);
    f.fp_shape = new btSphereShape(f.particle_radius);
    f.fp_shape->calculateLocalInertia(f.particle_mass, fp_inertia);
    for (k=0; (z = 1.633 * k * f.particle_radius) < dz; k++)
	for (j=0; (y = (j + (k%2)/3.0) * 1.732 * f.particle_radius) < dy; j++)
	    for (i=0; (x = (2*i + ((j+k)%2)) * f.particle_radius) < dx; i++, pi++) {
		btDefaultMotionState* fp_motion_state = new btDefaultMotionState
		    (btTransform(btQuaternion(0, 0, 0, 1), btVector3(x_min+x, y_min+y, z_min+z)));
		btRigidBody::btRigidBodyConstructionInfo fp_ci
		    (f.particle_mass, fp_motion_state, f.fp_shape, fp_inertia);
		fp_ci.m_restitution = 1.0;
		f.particles[pi].rigid_body = new btRigidBody(fp_ci);
		f.particles[pi].rigid_body->setLinearVelocity(btVector3(0, 0, 0));
		f.particles[pi].rigid_body->setAngularFactor(btVector3(0, 0, 0));
	    }
    return f;
}


int delete_particle(particle* pp)
{
    delete pp->rigid_body->getMotionState();
    delete pp->rigid_body;
    return 0;
}

int delete_fluid(fluid fluid)
{
    delete fluid.fp_shape;
    free(fluid.particles);
    return 0;
}

btVector3 particle_position (particle* pp)
{
    btTransform tf;
    pp->rigid_body->getMotionState()->getWorldTransform(tf);
    return btVector3(tf.getOrigin().getX(), tf.getOrigin().getY(), tf.getOrigin().getZ());
}

btVector3 particle_velocity (particle* pp)
{
    return pp->rigid_body->getLinearVelocity();
}
