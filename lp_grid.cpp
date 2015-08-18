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

int linear_id (lp_grid lpg, int i, int j, int k)
{
    return (i*lpg.y*lpg.z + j*lpg.z + k);
}

index3 inv_linear_id (lp_grid lpg, int linear_id)
{
    index3 id;
    int yz = lpg.y * lpg.z;
    id.i = linear_id / yz;
    int yz_rem = linear_id % yz;
    id.j = yz_rem / lpg.z;
    id.k = yz_rem % lpg.z;
    return id;
}

int rotate(lp_grid lpg, unsigned int origin, unsigned int destination)
{
    if ((origin>=lpg.particle_count) or (destination>=lpg.particle_count))
	return 1;

    unsigned int i;
    particle* tmp_storage;

    if (origin<destination) {
	tmp_storage = lpg.linear_grid[origin];
	for (i=origin+1; i<=destination; i++) {
	    lpg.linear_grid[i-1] = lpg.linear_grid[i];
	}
	lpg.linear_grid[destination] = tmp_storage;
	for (i=0; i<lpg.cell_count; i++) {
	    if (lpg.anchors[i].linear_id>=(origin+1))
		break;
	}
	while (lpg.anchors[i].linear_id<=destination) {
	    lpg.anchors[i++].linear_id--;
	}
	return 0;
    }

    if (origin>destination) {
	tmp_storage = lpg.linear_grid[origin];
	for (i=origin-1; i>=destination; i--) {
	    lpg.linear_grid[i+1] = lpg.linear_grid[i];
	}
	lpg.linear_grid[destination] = tmp_storage;
	for (i=0; i<lpg.cell_count; i++) {
	    if (lpg.anchors[i].linear_id>=destination)
		break;
	}
	while (lpg.anchors[i].linear_id<=(origin-1)) {
	    lpg.anchors[i++].linear_id++;
	}
	return 0;

    }
}

lp_grid make_lp_grid (btVector3 origin, float step, int x, int y, int z, std::vector<particle> particles)
{
    unsigned int cell_count = x*y*z;
    lp_grid lpg; unsigned int i,j,k, cid, pid;
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

    lpg.linear_map = (unsigned int*) malloc(cell_count * sizeof(unsigned int));
    lpg.anchors = (anchor*) malloc(cell_count * sizeof(anchor));
    lpg.linear_grid = (particle**) malloc (lpg.particle_count * sizeof(particle*));

    // Populate LINEAR_MAP
    Point p;
    for (unsigned int pid=0; pid<points.size(); pid++) {
	p = points[pid];
	lpg.linear_map[linear_id(lpg, p.x(), p.y(), p.z())] = pid;
    }

    // Initialize ANCHOR.LINEAR_ID's to 0 and ANCHOR.GRID_ID's to respective indices
    for (i=0; i<x; i++) {
	for (j=0; j<y; j++) {
	    for (k=0; k<z; k++) {
		cid = linear_id(lpg, i, j, k);
		lpg.anchors[cid].linear_id = 0;
		lpg.anchors[cid].grid_id.i = i;
		lpg.anchors[cid].grid_id.j = j;
		lpg.anchors[cid].grid_id.k = k;
	    }
	}
    }

    // populate LINEAR_GRID
    btVector3 position;
    for (pid=0; pid<particles.size(); pid++) {
	position = particle_position(particles[pid]);
	i = floor(position.getX()/lpg.step);
	j = floor(position.getY()/lpg.step);
	k = floor(position.getZ()/lpg.step);
	lpg.linear_grid[lpg.particle_count-1] = &particles[pid];
	// origin = last particle of last cell (always empty if grid not full)
	// destination = last particle of target cell
	rotate(lpg, lpg.anchors[linear_id(lpg, i, j, k)+1].linear_id-1, lpg.particle_count-1);
    }
    return lpg;
}

int update_lp_grid (lp_grid lpg)
{
    return 0;
}
