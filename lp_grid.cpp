#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/spatial_sort.h>
#include <math.h>

#include "lp_grid.h"
#include "utilities.h"

using CGAL::Cartesian;
using CGAL::Point_3;

typedef Cartesian<int> Kernel;
typedef Point_3<Kernel> Point;

uint linearize_address (lp_grid lpg, uint i, uint j, uint k)
{
    return (i*lpg.y*lpg.z + j*lpg.z + k);
}

uint particle_laddress (particle p, lp_grid lpg)
{
    btVector3 position = particle_position(p);
    uint i = floor(position.getX()/lpg.step);
    uint j = floor(position.getY()/lpg.step);
    uint k = floor(position.getZ()/lpg.step);
    return linearize_address(lpg, i, j, k);
}

// Move particle when ORIGIN < DESTINATION (= lpg.anchors[target]-1).
int move_particle_right(lp_grid lpg, uint origin, uint destination)
{
    uint i;
    particle* tmp_storage = lpg.particles[origin];
    for (i=origin+1; i<=destination; i++) {
	lpg.particles[i-1] = lpg.particles[i];
    }
    lpg.particles[destination] = tmp_storage;

    // Adjustment of anchors inbetween
    for (i=0; i<lpg.cell_count; i++) {
	if (lpg.anchors[i] >= (origin+1))
	    break;
    }
    while (lpg.anchors[i] <= destination) {
	lpg.anchors[i++]--;
    }
    return 0;
}

// Move particle when DESTINATION (= lpg.anchors[target+1]) < ORIGIN
int move_particle_left(lp_grid lpg, uint origin, uint destination)
{
    uint i;
    particle* tmp_storage = lpg.particles[origin];
    for (i=origin-1; i>=destination; i--) {
	lpg.particles[i+1] = lpg.particles[i];
    }
    lpg.particles[destination] = tmp_storage;

    // Adjustment of anchors inbetween
    for (i=0; i<lpg.cell_count; i++) {
	if (lpg.anchors[i] >= destination)
	    break;
    }
    while (lpg.anchors[i] <= (origin-1)) {
	lpg.anchors[i++]++;
    }
    return 0;
}

int allocate_lp_grid (lp_grid* lpg, std::vector<particle> particles, uint desired_cell_count)
{
    lpg->particle_count = particles.size();
    // Extract PARTICLES aabb
    btVector3 position;
    float x_min, y_min, z_min, x_max, y_max, z_max;
    for (uint pind=0; pind<lpg->particle_count; pind++) {
	position = particle_position(particles[pind]);
	if (position.getX() < x_min) x_min = position.getX();
	if (position.getY() < y_min) y_min = position.getY();
	if (position.getZ() < z_min) z_min = position.getZ();
	if (position.getX() > x_max) x_max = position.getX();
	if (position.getY() > y_max) y_max = position.getY();
	if (position.getZ() > z_max) z_max = position.getZ();
    }
    
    float x_diff = x_max - x_min;
    float y_diff = y_max - y_min;
    float z_diff = z_max - z_min;
    
    float aabb_volume = x_diff * y_diff * z_diff;
    float cell_volume = aabb_volume / desired_cell_count;
    lpg->step = cbrtf(cell_volume);
    lpg->x = ceil(x_diff / lpg->step)+1;
    lpg->y = ceil(y_diff / lpg->step)+1;
    lpg->z = ceil(z_diff / lpg->step)+1;
    lpg->cell_count = lpg->x * lpg->y * lpg->z;
    lpg->origin = btVector3(x_min - lpg->step, y_min - lpg->step, z_min - lpg->step);
    
    // Allocate memory for the 3 arrays
    lpg->map	    = (uint*) malloc(lpg->cell_count * sizeof(uint));
    lpg->anchors    = (uint*) malloc((lpg->cell_count+1) * sizeof(uint));
    lpg->particles  = (particle**) malloc ((lpg->particle_count+1) * sizeof(particle*));

    return 0;
}

lp_grid make_lp_grid (std::vector<particle> particles, uint desired_cell_count)
{
    // Grid parameter initialization/allocation
    lp_grid lpg; uint i,j,k;
    allocate_lp_grid(&lpg, particles, desired_cell_count);

    // Make points in centers of cells and spatially sort them
    std::vector<Point> points;
    for (i=0; i<lpg.x; i++) {
	for (j=0; j<lpg.y; j++) {
	    for (k=0; k<lpg.z; k++) {
		points.push_back(Point(i, j, k));
	    }
	}
    }
    CGAL::hilbert_sort(points.begin(), points.end());

    // Populate MAP
    Point p;
    uint cind;		// cind = cell index
    for (cind=0; cind<points.size(); cind++) {
	p = points[cind];
	lpg.map[linearize_address(lpg, p.x(), p.y(), p.z())] = cind;
    }

    // ANCHOR initialization
    for (i=0; i<lpg.x; i++) {
	for (j=0; j<lpg.y; j++) {
	    for (k=0; k<lpg.z; k++) {
		cind = lpg.map[linearize_address(lpg, i, j, k)];
		lpg.anchors[cind] = 0;
	    }
	}
    }
    lpg.anchors[lpg.cell_count] = lpg.particle_count;

    // Populate PARTICLES
    uint pind;		// pind = particle index
    for (pind=0; pind<lpg.particle_count; pind++) {
	// destination = lpg.anchors[target+1]-1 (last particle of target cell)
	cind = lpg.map[particle_laddress(particles[pind], lpg)];
	uint destination = lpg.anchors[cind+1]-1;
	// origin = the extra last slot in PARTICLES array
	lpg.particles[lpg.particle_count] = &particles[pind];
	move_particle_left(lpg, lpg.particle_count-1, destination);
    }

    return lpg;
}

particle_range get_cell(lp_grid lpg, uint i, uint j, uint k)
{
    particle_range pr;
    uint cind = lpg.map[linearize_address(lpg, i, j, k)];
    pr.start = lpg.particles[lpg.anchors[cind]];
    pr.end = lpg.particles[lpg.anchors[cind+1]];
    return pr;
}

int update_lp_grid (lp_grid lpg)
{
    for (uint origin_cind=0; origin_cind<lpg.cell_count; origin_cind++) {
	for (uint pind=lpg.anchors[origin_cind];
	     pind<lpg.anchors[origin_cind+1];
	     pind++) {
	    // destination_cind is the index of the cell to which the current particle
	    // actully belongs according to its position in space
	    uint destination_cind = lpg.map[particle_laddress(*lpg.particles[pind], lpg)];
	    if (origin_cind != destination_cind) {
		if (origin_cind < destination_cind)
		    move_particle_right(lpg, pind, lpg.anchors[destination_cind]);
		else move_particle_left(lpg, pind, lpg.anchors[destination_cind+1]-1);
	    }
	}
    }
    return 0;
}
