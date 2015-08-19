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

int move_particle_right(lp_grid lpg, unsigned int origin, unsigned int destination)
{
    unsigned int i;
    particle* tmp_storage = lpg.particles[origin];
    for (i=origin+1; i<=destination; i++) {
	lpg.particles[i-1] = lpg.particles[i];
    }
    lpg.particles[destination] = tmp_storage;
    for (i=0; i<lpg.cell_count; i++) {
	if (lpg.anchors[i].cell_start>=(origin+1))
	    break;
    }
    while (lpg.anchors[i].cell_start<=destination) {
	lpg.anchors[i++].cell_start--;
    }
    return 0;
}

int move_particle_left(lp_grid lpg, unsigned int origin, unsigned int destination)
{
    unsigned int i;
    particle* tmp_storage = lpg.particles[origin];
    for (i=origin-1; i>=destination; i--) {
	lpg.particles[i+1] = lpg.particles[i];
    }
    lpg.particles[destination] = tmp_storage;
    for (i=0; i<lpg.cell_count; i++) {
	if (lpg.anchors[i].cell_start>=destination)
	    break;
    }
    while (lpg.anchors[i].cell_start<=(origin-1)) {
	lpg.anchors[i++].cell_start++;
    }
    return 0;
}


lp_grid make_lp_grid (btVector3 origin, float step,
		      unsigned int x, unsigned int y, unsigned int z,
		      std::vector<particle> particles)
{
    // Parameter passing
    unsigned int cell_count = x*y*z;
    lp_grid lpg; unsigned int i,j,k;
    lpg.particle_count = particles.size();
    lpg.x = x; lpg.y = y; lpg.z = z; lpg.cell_count = x*y*z;
    lpg.origin = origin; lpg.step = step;

    // Make points in centers of cells...
    std::vector<Point> points;
    for (i=0; i<x; i++) {
	for (j=0; j<y; j++) {
	    for (k=0; k<z; k++) {
		points.push_back(Point(i, j, k));
	    }
	}
    }
    // ...and let CGAL spatially sort them.
    CGAL::hilbert_sort(points.begin(), points.end());

    lpg.map = (unsigned int*) malloc(cell_count * sizeof(unsigned int));
    lpg.anchors = (anchor*) malloc((1+cell_count) * sizeof(anchor));
    lpg.particles = (particle**) malloc (lpg.particle_count * sizeof(particle*));

    // Populate MAP
    Point p; unsigned int aid;	// aid = anchor index
    for (aid=0; aid<points.size(); aid++) {
	p = points[aid];
	lpg.map[linearize_address(lpg, p.x(), p.y(), p.z())] = aid;
    }

    // ANCHOR initialization
    for (i=0; i<x; i++) {
	for (j=0; j<y; j++) {
	    for (k=0; k<z; k++) {
		aid = lpg.map[linearize_address(lpg, i, j, k)];
		lpg.anchors[aid].cell_start = 0;
		lpg.anchors[aid].cell_laddress = linearize_address(lpg, i, j, k);
	    }
	}
    }
    lpg.anchors[lpg.cell_count].cell_start = lpg.particle_count;
    lpg.anchors[lpg.cell_count].cell_laddress = lpg.cell_count+1;

    // Populate PARTICLES
    unsigned int pid;		// pid = particle index
    for (pid=0; pid<particles.size(); pid++) {
	// destination = last particle of target cell
	unsigned int destination = lpg.anchors[particle_laddress(particles[pid], lpg)+1].cell_start-1;
	// origin = last particle of last cell (always empty if grid not full)
	lpg.particles[lpg.particle_count-1] = &particles[pid];
	move_particle_left(lpg, lpg.particle_count-1, destination);
    }
    return lpg;
}

int update_lp_grid (lp_grid lpg)
{
    for (unsigned int cid=0; cid<lpg.cell_count; cid++) {
	for (unsigned int pid=lpg.anchors[cid].cell_start; pid<lpg.anchors[cid+1].cell_start; pid++) {
	    unsigned int origin_laddress = lpg.anchors[cid].cell_laddress;
	    unsigned int destination_laddress = particle_laddress(*lpg.particles[pid], lpg);
	    if (origin_laddress != destination_laddress) {
		if (lpg.map[origin_laddress] < lpg.map[destination_laddress])
		    move_particle_left(lpg, pid, lpg.anchors[cid].cell_start);
		else move_particle_right(lpg, pid, lpg.anchors[cid+1].cell_start-1);
	    }
	}
    }
    return 0;
}

particle_range get_cell(lp_grid lpg, unsigned int i, unsigned int j, unsigned int k)
{
    particle_range pr;
    unsigned int aid = lpg.map[linearize_address(lpg, i, j, k)];
    pr.start = lpg.particles[lpg.anchors[aid].cell_start];
    pr.end = lpg.particles[lpg.anchors[aid+1].cell_start-1];
    return pr;
}

    
