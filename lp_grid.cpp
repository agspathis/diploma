#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/spatial_sort.h>

#include "lp_grid.h"
#include "utilities.h"

using CGAL::Cartesian;
using CGAL::Point_3;

typedef Cartesian<ulong> Kernel;
typedef Point_3<Kernel> Point;

ulong linearize_address (lp_grid lpg, ulong i, ulong j, ulong k)
{
    return (i*lpg.y*lpg.z + j*lpg.z + k);
}

ulong particle_laddress (particle p, lp_grid lpg)
{
    btVector3 position = particle_position(p) - lpg.origin;
    ulong i = floor(position.getX()/lpg.step);
    ulong j = floor(position.getY()/lpg.step);
    ulong k = floor(position.getZ()/lpg.step);
    return linearize_address(lpg, i, j, k);
}

// Move particle when ORIGIN < DESTINATION (= lpg.anchors[target]-1).
int move_particle_right(lp_grid lpg, ulong origin, ulong destination)
{
    ulong i;
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
int move_particle_left(lp_grid lpg, ulong origin, ulong destination)
{
    ulong i;
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

int allocate_lp_grid (lp_grid* lpg, aabb domain, float step, std::vector<particle> particles)
{
    lpg->step = step;
    lpg->particle_count = particles.size();
    lpg->x = ceil((domain.max.getX() - domain.min.getX()) / step)+1;
    lpg->y = ceil((domain.max.getY() - domain.min.getY()) / step)+1;
    lpg->z = ceil((domain.max.getZ() - domain.min.getZ()) / step)+1;
    lpg->cell_count = lpg->x * lpg->y * lpg->z;
    lpg->origin = btVector3(domain.min.getX() - step/2,
			    domain.min.getY() - step/2,
			    domain.min.getZ() - step/2);

    // Allocate memory for the 3 arrays
    lpg->map	    = (ulong*) malloc(lpg->cell_count * sizeof(ulong));
    lpg->anchors    = (ulong*) malloc((lpg->cell_count+1) * sizeof(ulong));
    lpg->particles  = (particle**) malloc ((lpg->particle_count+1) * sizeof(particle*));

    return 0;
}

lp_grid make_lp_grid (aabb domain, float step, std::vector<particle> particles)
{
    // Grid parameter initialization/allocation
    lp_grid lpg; ulong i,j,k;
    allocate_lp_grid(&lpg, domain, step, particles);

    printf("x=%lu, y=%lu, z=%lu\n", lpg.x, lpg.y, lpg.z);
    printf("cell_count=%lu\n", lpg.cell_count);

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

    // Initialize MAP
    ulong cind;			// cind = cell index
    for (cind=0; cind<points.size(); cind++) {
	lpg.map[linearize_address(lpg, points[cind].x(), points[cind].y(), points[cind].z())] = cind;
    }

    // Initialize ANCHORS to 0
    for (cind=0; cind<=lpg.cell_count; cind++) lpg.anchors[cind] = 0;

    // Populate PARTICLES
    ulong pind;		// pind = particle index
    for (pind=0; pind<lpg.particle_count; pind++) {
    	// destination = lpg.anchors[target+1]-1 (last particle of target cell)
    	cind = lpg.map[particle_laddress(particles[pind], lpg)];
    	ulong destination = lpg.anchors[cind+1]-1;
    	// origin = the extra last slot in PARTICLES array
    	lpg.particles[lpg.particle_count] = &particles[pind];
	move_particle_left(lpg, lpg.particle_count-1, destination);
    }

    return lpg;
}

particle_range get_cell(lp_grid lpg, ulong i, ulong j, ulong k)
{
    particle_range pr;
    ulong cind = lpg.map[linearize_address(lpg, i, j, k)];
    pr.start = lpg.particles[lpg.anchors[cind]];
    pr.end = lpg.particles[lpg.anchors[cind+1]];
    return pr;
}

int update_lp_grid (lp_grid lpg)
{
    for (ulong origin_cind=0; origin_cind<lpg.cell_count; origin_cind++) {
	for (ulong pind=lpg.anchors[origin_cind];
	     pind<lpg.anchors[origin_cind+1];
	     pind++) {
	    // destination_cind is the index of the cell to which the current particle
	    // actully belongs according to its position in space
	    ulong destination_cind = lpg.map[particle_laddress(*lpg.particles[pind], lpg)];
	    if (origin_cind != destination_cind) {
		if (origin_cind < destination_cind)
		    move_particle_right(lpg, pind, lpg.anchors[destination_cind]);
		else move_particle_left(lpg, pind, lpg.anchors[destination_cind+1]-1);
	    }
	}
    }
    return 0;
}
