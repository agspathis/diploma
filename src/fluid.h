#ifndef FLUID_H
#define FLUID_H

#include "terrain.h"

#define G 9.81
#define WATER_DENSITY 1000
#define WATER_DYNAMIC_VISCOSITY 0.001
#define WATER_IDEAL_K 200
#define WATER_TAIT_B 100
#define COURANT 0.5
#define MAX_DENSITY_FLUCTUATION 0.01
#define FLUID_INIT_VEL btVector3(0, 0, -10)

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
