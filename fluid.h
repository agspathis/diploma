#ifndef FLUID_H
#define FLUID_H

#include "terrain.h"

#define G 9.81
#define WATER_DENSITY 1000
#define WATER_DV 0.001
#define MAX_DENSITY_FLUCTUATION 0.01

struct particle {
    long id;
    int samples;
    float density;
    float pressure;
    float p_d2;
    btVector3 pforce;
    btVector3 vforce;
    btRigidBody* rigid_body;
};

struct fluid {
    /* set in construction */
    float density;
    long particle_count;
    float particle_mass;
    float particle_radius;
    float dynamic_viscosity;
    float smoothing_radius;
    /* set in calibration */
    float dt;
    float tait_b;
    float ideal_k;
    float max_samples;
    /* DENSITY_FACTOR = F.DENSITY / (density kernel estimation) */
    float density_factor;
    /* fluid data */
    btCollisionShape* fp_shape;
    particle* particles;
};

fluid make_fluid(aabb aabb, long desired_particle_count, long desired_sample_count);

int delete_particle(particle* pp);

int delete_fluid(fluid fluid);

btVector3 particle_position (particle* pp);

btVector3 particle_velocity (particle* pp);

#endif
