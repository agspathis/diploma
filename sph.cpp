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

int get_neighbour_cells(lp_grid lpg, long xi, long yi, long zi, std::vector<particle*> &neighbours)
{
    return 0;
}

int segment_interactions (lp_grid lpg, long xsi, long ysi, long zsi,
			  std::vector<std::vector<interaction> > &interactions)
{
    std::vector<particle> neighbours;
    std::vector<interaction> segment_interactions;
    for (int xi=xsi*lpg.xsl; xi<(xsi+1)*lpg.xsl; xi++) {
	for (int yi=ysi*lpg.ysl; yi<(ysi+1)*lpg.ysl; yi++) {
	    for (int zi=zsi*lpg.zsl; zi<(zsi+1)*lpg.zsl; zi++) {
		cell ccell = get_cell(lpg, xi, yi, zi);
		get_neighbour_cells(lpg, xi, yi, zi, neighbours);
		particle* cpp=ccpp.start; // cpp = center particle pointer
		while (cpp < cell.end;) {
		    btVector3 pos0 = particle_position(*cpp);
		    // npi = neighbour particle index
		    for (int npi=0; npi<neighbours.size(); npi++) {
			btVector3 pos1 = particle_position(neighbours[npi]);
			btVector3 direction = pos1-pos0;
			btVector3 distance = direction.lenght();
			if (distance < lpg.step) {
			    interaction interaction;
			    interaction.p0 = cpp;
			    interaction.p1 = &neighbours[npi];
			    interaction.distance = distance;
			    interaction.direction = direction.normalize();
			    segment_interactions.push_back(interaction);
			}
		    }
		    cpp++;
		}
		neighbours.clear();
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
	threads[i] = std::thread (segment_interactions, particles_per_thread, i,
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
