#include <cmath>
#include <thread>

#include "terrain.h"
#include "sph.h"
#include "utilities.h"

#define THREADS 100

float poly_6 (float r, float h)
{
    float result = 0;
    if (r < h) result = (315/(64*M_PI*pow(h, 9)))*pow(h*h-r*r, 3);
    return result;
}

float spiky (float r, float h)
{
    float result = 0;
    if (r < h) result = (15/(M_PI*pow(h, 6)))*pow(h-r, 3);
    return result;
}

float viscosity (float r, float h)
{
    float result = 0;
    if (r < h) result = (15/(2*M_PI*pow(h, 3)))*(-(pow(r, 3)/(2*pow(h, 3)))
						 +(pow(r, 2)/(pow(h, 2)))
						 +(h/(2*r))
						 - 1);
    return result;
}

// std::vector<interaction> forward_cell_interactions (grid g, int i, int j, int k,
// 						    float smoothing_length
						    

int interactions_in_segment (int particles_per_thread, int segment, std::vector<particle>particles,
			     float smoothing_length, std::vector<std::vector<interaction> > &interactions)
{
    btVector3 posi;
    btVector3 posj;
    btVector3 direction;
    btScalar distance;
    std::vector<interaction> segment_interactions;

    int start_ind = particles_per_thread * segment;
    int end_ind = ((particles_per_thread*(segment+1) < particles.size()) ?
		   start_ind+particles_per_thread :
		   particles.size());
    for (int i=start_ind; i<end_ind; i++) {
	posi = particle_position(particles[i]);
	for (int j=i+1; j<particles.size(); j++) {
	    posj = particle_position(particles[j]);
	    direction = posj-posi;
	    distance = direction.length();
	    if (distance < smoothing_length) {
		interaction interaction;
		interaction.i = i;
		interaction.j = j;
		interaction.distance = distance;
		interaction.direction = direction.normalize();
		segment_interactions.push_back(interaction);
	    }
	}
    }

    interactions[segment] = segment_interactions;

    return 0;
}

std::vector< std::vector<interaction> > collect_interactions_parallel (std::vector<particle> particles, float smoothing_length)
{
    int particle_count = particles.size();
    int particles_per_thread = ceil (particle_count/ (float) THREADS);
    std::vector< std::vector<interaction> > interactions (THREADS);
    std::vector< std::thread > threads (THREADS);
    for (int i=0; i<threads.size(); i++) {
	threads[i] = std::thread (interactions_in_segment, particles_per_thread, i,
				  particles, smoothing_length, std::ref(interactions));
    }
    for (int i=0; i<threads.size(); i++) threads[i].join();
    // DEBUG START
    int interaction_count=0;
    for (int icc=0; icc<interactions.size(); icc++) {
	interaction_count += interactions[icc].size();
    }
    // DEBUG END
    return interactions;
}

int clear_particle_data(std::vector<particle> &particles)
{
    for (int i=0; i<particles.size(); i++) {
	particles[i].samples = 0;
	particles[i].density = 0;
	particles[i].pressure = 0;
    }
    return 0;
}

int apply_sph_forces(std::vector<particle> particles,
		     float smoothing_length,
		     btScalar particle_mass)
{
    clear_particle_data(particles);
    std::vector< std::vector<interaction> > interactions =
	collect_interactions_parallel(particles, smoothing_length);
    return 0;
}

int apply_interactions (std::vector< std::vector<interaction> > interactions,
			std::vector<particle> particles,
			float smoothing_length,
			btScalar particle_mass)
{
    return 0;
}
