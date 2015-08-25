#include <cmath>
#include <thread>

#include "terrain.h"
#include "sph.h"
#include "utilities.h"

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

int clear_particle_data(lp_grid lpg)
{
    for (long pi=0; pi<lpg.particle_count; pi++) {
	lpg.particles[pi]->samples = 0;
	lpg.particles[pi]->density = 0;
	lpg.particles[pi]->pressure = 0;
    }
    return 0;
}

int get_neighbour_cells(lp_grid lpg, long xi, long yi, long zi, std::vector<particle*> &neighbours)
{
    return 0;
}

int segment_interactions (lp_grid lpg, long xsi, long ysi, long zsi, long sii,
			  std::vector<std::vector<interaction> > &interactions)
{
    long ixi = xsi*lpg.xss;
    long iyi = ysi*lpg.yss;
    long izi = zsi*lpg.zss;
    long txi = std::min((xsi+1)*lpg.xss, lpg.x);
    long tyi = std::min((ysi+1)*lpg.yss, lpg.y);
    long tzi = std::min((zsi+1)*lpg.zss, lpg.z);
    std::vector<particle*> neighbours;
    std::vector<interaction> segment_interactions;
    for (long xi=ixi; xi<txi; xi++) 
	for (long yi=iyi; yi<tyi; yi++) 
	    for (long zi=izi; zi<tzi; zi++) {
		// printf("%lu, %lu, %lu\n", xi, yi, zi);
		cell ccell = get_cell(lpg, xi, yi, zi);
		get_neighbour_cells(lpg, xi, yi, zi, neighbours);
		anchor cca = ccell.start; // cca = center cell anchor
		while (cca < ccell.end) {
		    btVector3 pos0 = particle_position(*cca);
		    // npi = neighbour particle index
		    for (long npi=0; npi<neighbours.size(); npi++) {
			btVector3 pos1 = particle_position(neighbours[npi]);
			btVector3 direction = pos1-pos0;
			btScalar distance = direction.length();
			if (distance < lpg.step) {
			    interaction interaction;
			    interaction.p0 = *cca;
			    interaction.p1 = neighbours[npi];
			    interaction.distance = distance;
			    interaction.direction = direction.normalize();
			    segment_interactions.push_back(interaction);
			}
		    }
		    cca++;
		}
		neighbours.clear();
	    }
    interactions[sii] = segment_interactions;
    return 0;
}

std::vector< std::vector<interaction> > collect_interactions_mt(lp_grid lpg)
{
    long xs_count = 1+lpg.x/lpg.xss;
    long ys_count = 1+lpg.y/lpg.yss;
    long zs_count = 1+lpg.z/lpg.zss;
    long thread_count = xs_count*ys_count*zs_count;
    std::vector< std::thread > threads (thread_count);
    std::vector< std::vector<interaction> > interactions (thread_count);
    long sii=0;			// sii = segment interactions index
    for (long xsi=0; xsi<xs_count; xsi++)
	for (long ysi=0; ysi<ys_count; ysi++)
	    for (long zsi=0; zsi<zs_count; zsi++, sii++)
		threads[sii] = std::thread (segment_interactions, lpg, xsi, ysi, zsi, sii,
					    std::ref(interactions));
    for (long ti=0; ti<thread_count; ti++) threads[ti].join();
    return interactions;
}

std::vector<interaction> collect_interactions_st (lp_grid lpg)
{
    std::vector<particle*> neighbours;
    std::vector<interaction> interactions;
    for (long xi=0; xi<lpg.x; xi++)
	for (long yi=0; yi<lpg.y; yi++)
	    for (long zi=0; zi<lpg.z; zi++) {
		// printf("%lu, %lu, %lu\n", xi, yi, zi);
		cell ccell = get_cell(lpg, xi, yi, zi);
		get_neighbour_cells(lpg, xi, yi, zi, neighbours);
		anchor cca = ccell.start; // cca = center cell anchor
		while (cca < ccell.end) {
		    btVector3 pos0 = particle_position(*cca);
		    // npi = neighbour particle index
		    for (long npi=0; npi<neighbours.size(); npi++) {
			btVector3 pos1 = particle_position(neighbours[npi]);
			btVector3 direction = pos1-pos0;
			btScalar distance = direction.length();
			if (distance < lpg.step) {
			    interaction interaction;
			    interaction.p0 = *cca;
			    interaction.p1 = neighbours[npi];
			    interaction.distance = distance;
			    interaction.direction = direction.normalize();
			    interactions.push_back(interaction);
			}
		    }
		    cca++;
		}
		neighbours.clear();
	    }
    return interactions;
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
    std::vector<std::vector<interaction> > interactions = collect_interactions_mt(lpg);
    //std::vector<interaction> interactions = collect_interactions_st(lpg);
    return 0;
}
