#include <thread>
#include <vector>

#include "sph.h"

// density kernel
float poly_6(float r, float h)
{
    float result = 0;
    // if (r<h) result = (315/(64*M_PI*pow(h, 9)))*pow(h*h-r*r, 3);
    if (r<h) result = pow(h*h-r*r, 3);
    return result;
}

// pressure kernel
float spiky(float r, float h)
{
    float result = 0;
    if (r<h) result = (15/(M_PI*pow(h, 6)))*pow(h-r, 3);
    return result;
}

float spiky_grad(float r, float h)
{
    float result = 0;
    if (r<h) result = (-45/(M_PI*pow(h, 6))) * pow((h-r), 2);
    // if (r<h) result = pow((h-r), 2);
    return result;
}

// viscosity kernel
float viscy(float r, float h)
{
    float result = 0;
    if (r<h) result = (15/(2*M_PI*pow(h, 3)))*(-(pow(r, 3)/(2*pow(h, 3)))
					       +(pow(r, 2)/(pow(h, 2)))
					       +(h/(2*r))
					       - 1);
    return result;
}

float viscy_lapl(float r, float h)
{
    float result = 0;
    // if (r<h) result = (45/(M_PI*pow(h, 6))) * (h-r);
    if (r<h) result = (h-r);
    return result;
}

// equation of state
float tait(float density, fluid f)
{
    return f.tait_b * (pow((density*f.particle_mass) / f.density,
			   TAIT_GAMMA)
		       -1);
}

float ideal(float density, fluid f)
{
    return f.ideal_k * (density - f.density);
}

void clear_particle_data(lp_grid lpg)
{
    for (long pi=0; pi<lpg.particle_count; pi++) {
	particle* pp = lpg.particles[pi];
	pp->samples = 0;
	pp->density = 0;
	pp->pressure = 0;
	pp->p_d2 = 0;
	pp->pforce = btVector3(0, 0, 0);
	pp->vforce = btVector3(0, 0, 0);
	pp->rigid_body->clearForces();
    }
}

void get_neighbour_cells(lp_grid lpg, long xi, long yi, long zi,
			 std::vector<particle*> &neighbours)
{
    // scan only half of surrounding cells to avoid duplication of work
    cell neighbour_cells[13];
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
}

void collect_segment_interactions (lp_grid lpg, long xsi, long ysi, long zsi, long sii,
				   std::vector< std::vector<interaction> > &interactions)
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
}

// The grid gets divided into segments (of dimensions LPG.XSS * LPG.YSS *
// LPG.ZSS cells) each of which is assigned to a separate thread for particle
// interaction scan. Each slot of vector SEGMENT_INTERACTIONS contains a vector
// with the interactions scanned in a segment. This vector is then flattened and
// returned.
std::vector<interaction> collect_interactions(lp_grid lpg)
{
    long xs_count = 1+lpg.x/lpg.xss;
    long ys_count = 1+lpg.y/lpg.yss;
    long zs_count = 1+lpg.z/lpg.zss;
    long thread_count = xs_count*ys_count*zs_count;
    std::vector< std::thread > threads (thread_count);
    std::vector< std::vector<interaction> > segment_interactions (thread_count);
    long sii=0;			// sii = segment interactions index
    for (long xsi=0; xsi<xs_count; xsi++)
	for (long ysi=0; ysi<ys_count; ysi++)
	    for (long zsi=0; zsi<zs_count; zsi++, sii++)
		threads[sii] = std::thread (collect_segment_interactions, lpg, xsi, ysi, zsi, sii,
					    std::ref(segment_interactions));
    for (long ti=0; ti<thread_count; ti++) threads[ti].join();
    // flatten interaction vector
    long ic=0;			// interaction count
    for(sii=0; sii<segment_interactions.size(); sii++)
	ic += segment_interactions[sii].size();
    std::vector<interaction> interactions (ic);
    long fii=0;			// flat interactions index
    for(sii=0; sii<segment_interactions.size(); sii++)
	for(long ii=0; ii<segment_interactions[sii].size(); ii++, fii++)
	    interactions[fii] = segment_interactions[sii][ii];
    return interactions;
}

void compute_densities(std::vector<interaction> interactions, float smoothing_radius)
{
    for(long ii=0; ii<interactions.size(); ii++) {
	float density_fraction = poly_6(interactions[ii].distance, smoothing_radius);
	interactions[ii].p0->density += density_fraction;
	interactions[ii].p0->samples++;
	interactions[ii].p1->density += density_fraction;
	interactions[ii].p1->samples++;
    }
}

// compute pressures from EOS and pressure/density^2 for each particle
void compute_pressures(fluid f)
{
    for(particle* pp=f.particles; pp<f.particles+f.particle_count; pp++) {
	pp->density += (f.max_samples - pp->samples) * f.avg_density_fraction;
	pp->density *= f.particle_mass * f.density_factor;
	pp->pressure = ideal(pp->density, f);
	pp->p_d2 = pp->pressure / pow(pp->density, 2);
    }
}

// forces arise from pressure difference and viscosity
void apply_forces(std::vector<interaction> interactions, fluid f)
{
    btVector3 dir;
    for (long ii=0; ii<interactions.size(); ii++) {
	interaction i = interactions[ii];
	float pf_fraction = 
	    (i.p0->p_d2 + i.p1->p_d2)
	    * spiky_grad(i.distance, f.smoothing_radius);
	float vf_fraction = 
	    f.dynamic_viscosity
	    * (particle_velocity(i.p0)- particle_velocity(i.p1)).length()
	    * viscy_lapl(i.distance, f.smoothing_radius);
	dir = i.direction;
	i.p0->pforce += (dir *= pf_fraction);
	i.p1->pforce -= dir;
	dir = i.direction;
	i.p0->vforce += (dir *= vf_fraction);
	i.p1->vforce -= dir;
    }
    for (particle* pp=f.particles; pp<f.particles + f.particle_count; pp++) {
	pp->rigid_body->applyCentralForce(pp->pforce + pp->vforce);
	// printf("Samples = %lu, density = %f\n", pp->samples, pp->density/f.density);
    }
}

void apply_sph(fluid_sim fsim)
{
    clear_particle_data(fsim.lpg);
    std::vector<interaction> interactions = collect_interactions(fsim.lpg);
    compute_densities(interactions, fsim.lpg.step);
    compute_pressures(fsim.f);
    apply_forces(interactions, fsim.f);
}

void adjust_fluid(fluid* f, lp_grid lpg, aabb faabb, aabb taabb)
{
    float height = faabb.max.getY() - taabb.min.getY();
    float v_max = sqrt(2 * G * height);
    float cs = v_max / sqrt(MAX_DENSITY_FLUCTUATION);
    f->tait_b = (f->density) * pow(cs, 2) / TAIT_GAMMA;
    f->ideal_k = 2000;

    // dt to match MAX_DENSITY_FLUCTUATION between ticks
    f->dt = MAX_DENSITY_FLUCTUATION * f->particle_radius / v_max;
    
    // measure and store max density and samples in the initial fluid
    // configuration in order to ajust later sph computations
    clear_particle_data(lpg);
    std::vector<interaction> interactions = collect_interactions(lpg);
    compute_densities(interactions, lpg.step);

    int max_samples = 0;
    float max_density = 0;
    for(particle* pp=f->particles; pp<f->particles + f->particle_count; pp++) {
	if (pp->samples > max_samples) {
	    max_samples = pp->samples;
	    max_density = pp->density;
	}
    }
    f->max_samples = max_samples;
    f->avg_density_fraction = max_density / max_samples;
    f->density_factor = f->density / (max_density * f->particle_mass);

    printf("FLUID ADJUSTMENT INFO:\n");
    printf("dt = %f\n", f->dt);
    printf("tait_b = %f\n", f->tait_b);
    printf("max_density = %f\n", max_density);
    printf("max_samples = %d\n", max_samples);
    printf("density_factor = %f\n\n", f->density_factor);
}
