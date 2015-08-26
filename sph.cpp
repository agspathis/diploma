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
    // scan only half of surrounding cells to avoid duplication of work
    std::vector<cell> neighbour_cells(13);
    neighbour_cells[ 0] = get_cell(lpg, xi-1, yi-1, zi-1);
    neighbour_cells[ 1] = get_cell(lpg, xi  , yi-1, zi-1);
    neighbour_cells[ 2] = get_cell(lpg, xi+1, yi-1, zi-1);
    neighbour_cells[ 3] = get_cell(lpg, xi-1, yi  , zi-1);
    neighbour_cells[ 4] = get_cell(lpg, xi  , yi  , zi-1);
    neighbour_cells[ 5] = get_cell(lpg, xi+1, yi  , zi-1);
    neighbour_cells[ 6] = get_cell(lpg, xi-1, yi+1, zi-1);
    neighbour_cells[ 7] = get_cell(lpg, xi  , yi+1, zi-1);
    neighbour_cells[ 8] = get_cell(lpg, xi+1, yi+1, zi-1);
    neighbour_cells[ 9] = get_cell(lpg, xi-1, yi-1, zi  );
    neighbour_cells[10] = get_cell(lpg, xi  , yi-1, zi  );
    neighbour_cells[11] = get_cell(lpg, xi+1, yi-1, zi  );
    neighbour_cells[12] = get_cell(lpg, xi-1, yi  , zi  );
    
    for(int nci=0; nci<13; nci++)
	for(anchor a=neighbour_cells[nci].start; a<neighbour_cells[nci].end; a++)
	    neighbours.push_back(*a);

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

std::vector< std::vector<interaction> > collect_interactions(lp_grid lpg)
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
    // DEBUG START
    int interaction_count=0;
    for (long icc=0; icc<interactions.size(); icc++)
	interaction_count += interactions[icc].size();
    printf("Found %lu interactions\n", interaction_count);
    // DEBUG END
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
    std::vector<std::vector<interaction> > interactions = collect_interactions(lpg);
    apply_interactions(interactions, lpg.step, particle_mass);
    return 0;
}
