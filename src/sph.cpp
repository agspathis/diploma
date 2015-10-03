#include <thread>
#include <vector>

#include "sph.h"

// density kernel (poly 6)
float dens_k(float r, float h)
{
    return (315/(64*M_PI*pow(h, 9)))*pow(h*h-r*r, 3);
}
// pressure kernel (spiky) gradient
float press_kg(float r, float h)
{
    return (-45/(M_PI*pow(h, 6))) * pow(h-r, 2);
}
// viscosity kernel laplacian
float visc_kl(float r, float h)
{
    return (45/(M_PI*pow(h, 6))) * (h-r);
}
// equation of state
float eos(float density, fluid f)
{
    // return f.tait_b * (pow((density * f.particle_mass) / f.density, 7) - 1);
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
    }
}

void neighbour_particles(lp_grid lpg, long xi, long yi, long zi,
			 std::vector<particle*> &neighbours)
{
    neighbours.clear();
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
    // segment looping limits
    long ixi = xsi*lpg.xss; long txi = std::min((xsi+1)*lpg.xss, lpg.x);
    long iyi = ysi*lpg.yss; long tyi = std::min((ysi+1)*lpg.yss, lpg.y);
    long izi = zsi*lpg.zss; long tzi = std::min((zsi+1)*lpg.zss, lpg.z);
    // interaction vars
    cell ccell;			 // center cell
    anchor cca0, cca1;		 // anchors to p0, p1
    btVector3 pos0, pos1;	 // positions of p0, p1
    btVector3 direction;
    btScalar distance;
    interaction inter;
    std::vector<particle*> neighbours;
    std::vector<interaction> segment_interactions;

    for (long zi=izi; zi<tzi; zi++)
	for (long yi=iyi; yi<tyi; yi++)
	    for (long xi=ixi; xi<txi; xi++) {
		ccell = get_cell(lpg, xi, yi, zi); // ccell = center cell
		// collect interactions with particles inside center cell
		cca0 = ccell.start;
		while (cca0 < ccell.end-1) {
		    pos0 = particle_position(*cca0);
		    for (cca1 = cca0+1; cca1 < ccell.end; cca1++) {
			pos1 = particle_position(*cca1);
			direction = pos1-pos0;
			distance = direction.length();
			if (distance < lpg.step) {
			    inter.p0 = *cca0;
			    inter.p1 = *cca1;
			    inter.distance = distance;
			    inter.direction = direction.normalize();
			    segment_interactions.push_back(inter);
			}
		    }
		    cca0++;
		}
		// collect interactions with particles in neighboring cells
		neighbour_particles(lpg, xi, yi, zi, neighbours);
		cca0 = ccell.start;
		while (cca0 < ccell.end) {
		    pos0 = particle_position(*cca0);
		    // npi = neighbour particle index
		    for (long npi=0; npi<neighbours.size(); npi++) {
			pos1 = particle_position(neighbours[npi]);
			direction = pos1-pos0;
			distance = direction.length();
			if (distance < lpg.step) {
			    inter.p0 = *cca0;
			    inter.p1 = neighbours[npi];
			    inter.distance = distance;
			    inter.direction = direction.normalize();
			    segment_interactions.push_back(inter);
			}
		    }
		    cca0++;
		}
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
	float density_fraction = dens_k(interactions[ii].distance, smoothing_radius);
	interactions[ii].p0->density += density_fraction;
	interactions[ii].p0->samples++;
	interactions[ii].p1->density += density_fraction;
	interactions[ii].p1->samples++;
    }
}

// compute pressures from equation of state and pressure/density^2
void compute_pressures(fluid f)
{
    for(particle* pp=f.particles; pp<f.particles+f.particle_count; pp++) {
	pp->density *= f.particle_mass * f.density_factor;
	// correction of undersampled density
	if (pp->samples < 0.9 * f.max_samples) pp->density = f.density;
	pp->pressure = eos(pp->density, f);
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
	    - f.particle_mass
	    * (i.p0->p_d2 + i.p1->p_d2)
	    * press_kg(i.distance, f.smoothing_radius);
	float vf_fraction =
	    f.dynamic_viscosity
	    * (particle_velocity(i.p0)- particle_velocity(i.p1)).length()
	    * visc_kl(i.distance, f.smoothing_radius);
	dir = i.direction;
	i.p0->pforce += (dir *= pf_fraction);
	i.p1->pforce -= dir;
	dir = i.direction;
	i.p0->vforce += (dir *= vf_fraction);
	i.p1->vforce -= dir;
    }
    for (particle* pp=f.particles; pp<f.particles + f.particle_count; pp++) {
	pp->rigid_body->applyCentralImpulse((pp->pforce + pp->vforce) * f.dt);
	// printf("Samples = %lu, density = %f\n", pp->samples, pp->density/f.density);
    }
}

void apply_sph(fluid_sim* fsimp)
{
    clear_particle_data(fsimp->lpg);
    std::vector<interaction> interactions = collect_interactions(fsimp->lpg);
    compute_densities(interactions, fsimp->lpg.step);
    compute_pressures(fsimp->f);
    apply_forces(interactions, fsimp->f);
}

// adjust fluid to terrain, compute sph parameters
void adjust_fluid(fluid* f, lp_grid lpg, aabb faabb, aabb taabb)
{
    float height = faabb.max.getY() - taabb.min.getY();
    float v_max = sqrt(2 * G * height);
    float cs = v_max / sqrt(MAX_DENSITY_FLUCTUATION);
    f->ideal_k = WATER_IDEAL_K;
    f->tait_b = WATER_TAIT_B;

    // dt to match MAX_DENSITY_FLUCTUATION between ticks
    f->dt = 0.1 * f->particle_radius / v_max;

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
    f->max_samples = (float) max_samples;
    f->density_factor = f->density / (max_density * f->particle_mass);

    printf("FLUID ADJUSTMENT INFO:\n");
    printf("dt = %f\n", f->dt);
    printf("max_density = %f\n", max_density);
    printf("max_samples = %d\n", max_samples);
    printf("density_factor = %f\n\n", f->density_factor);
}

// color field value at position CPOS, smoothing radius H
float sum_cf(std::vector<particle*> neighbours, btVector3 cpos, float h)
{
    float cf = 0;
    for (int pi=0; pi<neighbours.size(); pi++) {
	float distance = (particle_position(neighbours[pi]) - cpos).length();
 	if (distance < h) cf += dens_k(distance, h);
    }
    return cf;
}

void compute_cf_segment(lp_grid lpg, long xsi, long ysi, long zsi)
{
    cell ccell;
    btVector3 cell_origin, cpos;
    std::vector<particle*> neighbours;
    float cf_step = lpg.step / lpg.cf_sdf;
    long ixi = xsi*lpg.xss; long txi = std::min((xsi+1)*lpg.xss, lpg.x);
    long iyi = ysi*lpg.yss; long tyi = std::min((ysi+1)*lpg.yss, lpg.y);
    long izi = zsi*lpg.zss; long tzi = std::min((zsi+1)*lpg.zss, lpg.z);
    for (long zi=izi; zi<tzi; zi++)
	for (long yi=iyi; yi<tyi; yi++)
	    for (long xi=ixi; xi<txi; xi++) {
		cell_origin = btVector3(xi*lpg.step, yi*lpg.step, zi*lpg.step);
		neighbour_particles(lpg, xi, yi, zi, neighbours);
		// add center cell particles to NEIGHBOURS
		ccell = get_cell(lpg, xi, yi, zi);
		for (anchor a=ccell.start; a<ccell.end; a++)
		    neighbours.push_back(*a);
		// loop over cf sampling points inside ccell
		for (int z_sd=0; z_sd<lpg.cf_sdf; z_sd++)
		    for (int y_sd=0; y_sd<lpg.cf_sdf; y_sd++)
			for (int x_sd=0; x_sd<lpg.cf_sdf; x_sd++) {
			    cpos = btVector3(x_sd, y_sd, z_sd);
			    cpos *= cf_step;
			    cpos += cell_origin;
			    lpg.color_field[linearize(lpg,
						      xi*lpg.cf_sdf + x_sd,
						      yi*lpg.cf_sdf + y_sd,
						      zi*lpg.cf_sdf + z_sd,
						      lpg.cf_sdf)]
				= sum_cf(neighbours, cpos, lpg.step);
			}
	    }
}

void compute_cf(lp_grid lpg)
{
    long ti = 0;
    long xs_count = 1+lpg.x/lpg.xss;
    long ys_count = 1+lpg.y/lpg.yss;
    long zs_count = 1+lpg.z/lpg.zss;
    long thread_count = xs_count*ys_count*zs_count;
    std::thread threads[thread_count];
    for (long xsi=0; xsi<xs_count; xsi++)
	for (long ysi=0; ysi<ys_count; ysi++)
	    for (long zsi=0; zsi<zs_count; zsi++, ti++)
		threads[ti] = std::thread(compute_cf_segment, lpg, xsi, ysi, zsi);
    for (long ti=0; ti<thread_count; ti++) threads[ti].join();
}
