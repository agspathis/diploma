#include <cmath>
#include <thread>

#include "terrain.h"
#include "sph.h"
#include "lp_grid.h"
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

int segment_interactions (lp_grid lpg, long xsi, long ysi, long zsi,
			  std::vector<std::vector<interaction> > &interactions)
{
    btVector3 pos0;
    btVector3 pos1;
    btVector3 direction;
    btScalar distance;
    std::vector<interaction> segment_interactions;
    for (int i=xsi*lpg.xsl; i<(xsi+1)*lpg.xsl; i++) {
	for (int j=ysi*lpg.ysl; i<(ysi+1)*lpg.ysl; i++) {
	    for (int i=zsi*lpg.zsl; i<(zsi+1)*lpg.zsl; i++) {
		particle_range cell = get_cell(lpg, i, j, k);
	    }
	}
    }
    
    for (int p0i=start_ind; i<end_ind; i++) {
	pos0 = particle_position(particles[i]);
	for (int p1i=i+1; j<particles.size(); j++) {
	    posj = particle_position(particles[j]);
	    direction = posj-posi;
	    distance = direction.length();
	    if (distance < lpg.step) {
		interaction interaction;
		interaction.p0 = i;
		interaction.p1 = j;
		interaction.distance = distance;
		interaction.direction = direction.normalize();
		segment_interactions.push_back(interaction);
	    }
	}
    }

    interactions[segment] = segment_interactions;

    return 0;
}

std::vector< std::vector<interaction> > collect_interactions_parallel(lp_grid lpg)
{
    int particles_per_thread = ceil (lpg.particle_count / (float) THREADS);
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

int clear_particle_data(lp_grid lpg)
{
    for (long pi=0; pi<lpg.particle_count; pi++) {
	lpg.particles[i]->samples = 0;
	lpg.particles[i]->density = 0;
	lpg.particles[i]->pressure = 0;
    }
    return 0;
}


int apply_interactions (std::vector< std::vector<interaction> > interactions,
			float smoothing_length,
			btScalar particle_mass)
{
    return 0;
}

int apply_sph_forces(lp_grid lpg, btScalar particle_mass)
{
    clear_particle_data(lpg);
    std::vector< std::vector<interaction> > interactions =
	collect_interactions_parallel(lpg);
    return 0;
}
