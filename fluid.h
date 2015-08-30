#ifndef FLUID_H
#define FLUID_H

#include "terrain.h"

struct particle {
    long id;
    int samples;
    float density;
    float pressure;
    float p_d2;
    btVector3 force;
    btRigidBody* rigid_body;
};

/*
  SAMPLE_DENSITY is the local density as estimated by summation over the samples
  inside SMOOTHING_RADIUS when the fluid is at rest, but divided by PARTICLE_MASS
 */
struct fluid {
    float particle_mass;
    float particle_radius;
    long particle_count;
    float density;
    float sample_density;
    float tait_b;
    float smoothing_radius;
    btCollisionShape* fp_shape;
    particle* particles;
    float dynamic_viscosity;
};

fluid make_fluid(aabb aabb, long desired_particle_count);

int delete_particle(particle* pp);

int delete_fluid(fluid fluid);

btVector3 particle_position (particle* pp);

#endif
