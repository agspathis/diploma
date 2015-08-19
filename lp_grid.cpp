#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/spatial_sort.h>

#include "lp_grid.h"
#include "utilities.h"

using CGAL::Cartesian;
using CGAL::Point_3;

typedef Cartesian<int> Kernel;
typedef Point_3<Kernel> Point;

unsigned int linearize_address (lp_grid lpg, unsigned int i, unsigned int j, unsigned int k)
{
    return (i*lpg.y*lpg.z + j*lpg.z + k);
}

unsigned int particle_laddress (particle p, lp_grid lpg)
{
    btVector3 position = particle_position(p);
    unsigned int i = floor(position.getX()/lpg.step);
    unsigned int j = floor(position.getY()/lpg.step);
    unsigned int k = floor(position.getZ()/lpg.step);
    return linearize_address(lpg, i, j, k);
}

// Move particle when ORIGIN < DESTINATION (= lpg.anchors[target]-1).
int move_particle_right(lp_grid lpg, unsigned int origin, unsigned int destination)
{
    unsigned int i;
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
int move_particle_left(lp_grid lpg, unsigned int origin, unsigned int destination)
{
    unsigned int i;
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


lp_grid make_lp_grid (btVector3 origin, float step,
		      unsigned int x, unsigned int y, unsigned int z,
		      std::vector<particle> particles)
{
    // Parameter passing
    lp_grid lpg; unsigned int i,j,k;
    lpg.origin = origin;
    lpg.step = step;
    lpg.x = x; lpg.y = y; lpg.z = z;
    lpg.cell_count = x*y*z;
    lpg.particle_count = particles.size();

    // Make points in centers of cells and spatially sort them
    std::vector<Point> points;
    for (i=0; i<x; i++) {
	for (j=0; j<y; j++) {
	    for (k=0; k<z; k++) {
		points.push_back(Point(i, j, k));
	    }
	}
    }
    CGAL::hilbert_sort(points.begin(), points.end());

    // Allocate memory for the 3 arrays
    lpg.map = (unsigned int*) malloc(lpg.cell_count * sizeof(unsigned int));
    lpg.anchors = (unsigned int*) malloc((lpg.cell_count+1) * sizeof(unsigned int));
    lpg.particles = (particle**) malloc ((lpg.particle_count+1) * sizeof(particle*));

    // Populate MAP
    Point p;
    unsigned int cind;		// cind = cell index
    for (cind=0; cind<points.size(); cind++) {
	p = points[cind];
	lpg.map[linearize_address(lpg, p.x(), p.y(), p.z())] = cind;
    }

    // ANCHOR initialization
    for (i=0; i<x; i++) {
	for (j=0; j<y; j++) {
	    for (k=0; k<z; k++) {
		cind = lpg.map[linearize_address(lpg, i, j, k)];
		lpg.anchors[cind] = 0;
	    }
	}
    }
    lpg.anchors[lpg.cell_count] = lpg.particle_count;

    // Populate PARTICLES
    unsigned int pind;		// pind = particle index
    for (pind=0; pind<particles.size(); pind++) {
	// destination = lpg.anchors[target+1]-1 (last particle of target cell)
	unsigned int destination = lpg.anchors[lpg.map[particle_laddress(particles[pind], lpg)]+1]-1;
	// origin = the extra last slot in PARTICLES array
	lpg.particles[lpg.particle_count] = &particles[pind];
	move_particle_left(lpg, lpg.particle_count-1, destination);
    }
    return lpg;
}

int update_lp_grid (lp_grid lpg)
{
    for (unsigned int origin_cind=0; origin_cind<lpg.cell_count; origin_cind++) {
	for (unsigned int pind=lpg.anchors[origin_cind];
	     pind<lpg.anchors[origin_cind+1];
	     pind++) {
	    // destination_cind is the index of the cell to which the current particle
	    // actully belongs according to its position in space
	    unsigned int destination_cind = lpg.map[particle_laddress(*lpg.particles[pind], lpg)];
	    if (origin_cind != destination_cind) {
		if (origin_cind < destination_cind)
		    move_particle_right(lpg, pind, lpg.anchors[destination_cind]);
		else move_particle_left(lpg, pind, lpg.anchors[destination_cind+1]-1);
	    }
	}
    }
    return 0;
}

particle_range get_cell(lp_grid lpg, unsigned int i, unsigned int j, unsigned int k)
{
    particle_range pr;
    unsigned int cind = lpg.map[linearize_address(lpg, i, j, k)];
    pr.start = lpg.particles[lpg.anchors[cind]];
    pr.end = lpg.particles[lpg.anchors[cind+1]];
    return pr;
}
