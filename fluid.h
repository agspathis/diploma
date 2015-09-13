#ifndef FLUID_H
#define FLUID_H

#include "terrain.h"

#define G 9.81
#define WATER_DENSITY 1000
#define WATER_DV 0.001
#define MAX_DENSITY_FLUCTUATION 1

struct particle {
    long id;
    int samples;
    float density;
    float pressure;
    float p_d2;
    btVector3 pforce;
    btVector3 vforce;
    btVector3 neighbour_location;
    btRigidBody* rigid_body;
};

/*
  DENSITY_FACTOR is determined from the fluid at rest and scales the density
  kernel estimation to match the actual density.
 */
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
    float ideal_k;
    int max_samples;
    float density_factor;
    float avg_density_fraction;
    /* fluid data */
    btCollisionShape* fp_shape;
    particle* particles;
};

fluid make_fluid(aabb aabb, long desired_particle_count);

int delete_particle(particle* pp);

int delete_fluid(fluid fluid);

btVector3 particle_position (particle* pp);

btVector3 particle_velocity (particle* pp);

#endif
