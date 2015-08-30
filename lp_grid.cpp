#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/spatial_sort.h>

#include "lp_grid.h"

using CGAL::Cartesian;
using CGAL::Point_3;

typedef Cartesian<long> Kernel;
typedef Point_3<Kernel> Point;

long linearize_address (lp_grid lpg, long i, long j, long k)
{
    // if any index is out of bounds in any of the 3 dimensions of the grid the
    // linear address (LPG.CELL_COUNT) of the last cell containing the off-grid
    // particles is returned
    if (i<0 || i>=lpg.x || j<0 || j>=lpg.y || k<0 || k>=lpg.z)
	return lpg.cell_count;
    else return (i*lpg.y*lpg.z + j*lpg.z + k);
}

long particle_laddress (lp_grid lpg, particle* pp)
{
    btVector3 relative_position = particle_position(pp) - lpg.origin;
    long i = floor(relative_position.getX()/lpg.step);
    long j = floor(relative_position.getY()/lpg.step);
    long k = floor(relative_position.getZ()/lpg.step);
    return linearize_address(lpg, i, j, k);
}

// Letters in abbreviations "ici", "tpi" etc. stand for
// i=initial, t=terminal, c=cell, p=particle, i=index
int insert_particle(lp_grid lpg, particle* pp)
{
    // terminal cell and particle index
    long tci = lpg.map[particle_laddress(lpg, pp)];
    long tpi = lpg.anchors[tci+1];
    // right shift before insertion
    for (long pi=lpg.particle_count-1; pi>=tpi; pi--)
	lpg.particles[pi+1] = lpg.particles[pi];
    lpg.particles[tpi] = pp;
    // anchor adjustment
    for (long ci=tci+1; ci<=lpg.cell_count; ci++) 
	lpg.anchors[ci]++;
    return 0;
}

int allocate_lp_grid (lp_grid* lpg, aabb domain, fluid fluid)
{
    lpg->step = fluid.smoothing_radius;
    lpg->particle_count = fluid.particle_count;
    lpg->x = ceil((domain.max.getX() - domain.min.getX()) / lpg->step)+1;
    lpg->y = ceil((domain.max.getY() - domain.min.getY()) / lpg->step)+1;
    lpg->z = ceil((domain.max.getZ() - domain.min.getZ()) / lpg->step)+1;
    lpg->xss = ceil(sqrt(lpg->x));
    lpg->yss = ceil(sqrt(lpg->y));
    lpg->zss = ceil(sqrt(lpg->z));
    lpg->cell_count = lpg->x * lpg->y * lpg->z;
    lpg->origin = btVector3(domain.min.getX() - lpg->step/2,
			    domain.min.getY() - lpg->step/2,
			    domain.min.getZ() - lpg->step/2);

    // Allocate memory for the 3 arrays
    lpg->map	    = (long*) malloc((lpg->cell_count+1) * sizeof(long));
    lpg->anchors    = (long*) malloc((lpg->cell_count+1) * sizeof(long));
    lpg->particles  = (particle**) malloc ((lpg->particle_count+1) * sizeof(particle*));

    return 0;
}

int verify_lp_grid(lp_grid lpg)
{
    for (long ici=0; ici<lpg.cell_count; ici++) {
	for (long ipi=lpg.anchors[ici]; ipi<lpg.anchors[ici+1]; ipi++) {
	    long tci = lpg.map[particle_laddress(lpg, lpg.particles[ipi])];
	    if (ici != tci) 
		printf("ERROR pid=%lu ipi=%lu ici=%lu tci=%lu\n",
		       lpg.particles[ipi]->id, ipi, ici, tci);
	}
    }
    return 1;
}

lp_grid make_lp_grid (aabb domain, fluid fluid)
{
    // Grid parameter initialization/allocation
    lp_grid lpg; long i,j,k;
    allocate_lp_grid(&lpg, domain, fluid);
    printf("LP grid info:\n");
    printf("x=%lu, y=%lu, z=%lu\n", lpg.x, lpg.y, lpg.z);
    printf("xss=%lu, yss=%lu, zss=%lu\n", lpg.xss, lpg.yss, lpg.zss);
    printf("Cell count=%lu\n", lpg.cell_count);
    printf("Particle count=%lu\n", lpg.particle_count);
    printf("\n");

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
    long ci;
    for (ci=0; ci<points.size(); ci++) {
	lpg.map[linearize_address(lpg, points[ci].x(), points[ci].y(), points[ci].z())] = ci;
    }
    lpg.map[lpg.cell_count] = lpg.cell_count;

    // Initialize ANCHORS to 0
    for (ci=0; ci<=lpg.cell_count; ci++) lpg.anchors[ci] = 0;

    // Populate PARTICLES
    // for (long pi=0; pi<lpg.particle_count; pi++) insert_particle(lpg, fluid.particles + pi);
    for (particle* pp=fluid.particles; pp<fluid.particles+lpg.particle_count; pp++)
	insert_particle(lpg, pp);

    if (verify_lp_grid(lpg)) printf("New LP grid verified successfully!\n");
    else printf("New LP grid is wrong...\n");

    return lpg;
}


int move_particle(lp_grid lpg, long ipi, long ici, long tci)
{
    long ci, pi, tpi;
    particle* tmp_storage = lpg.particles[ipi];

    if (ici < tci) {
	long tpi = lpg.anchors[tci]-1;
	for (pi=ipi; pi<tpi; pi++) lpg.particles[pi] = lpg.particles[pi+1];
	for (ci=ici+1; ci<=tci; ci++) lpg.anchors[ci] -= 1;
    }
    else {
	tpi = lpg.anchors[tci+1];
	for (pi=ipi; pi>tpi; pi--) lpg.particles[pi] = lpg.particles[pi-1];
	for (ci=tci+1; ci<=ici; ci++) lpg.anchors[ci]++;
    }
    lpg.particles[tpi] = tmp_storage;
    printf("MOVE pid=%lu ipi=%lu tpi=%lu ici=%lu tci=%lu\n",
	   tmp_storage->id, ipi, tpi, ici, tci);
    return 0;
}

int update_lp_grid (lp_grid lpg)
{
    for (long ici=0; ici<lpg.cell_count; ici++) {
	for (long ipi=lpg.anchors[ici]; ipi<lpg.anchors[ici+1]; ipi++) {
	    long tci = lpg.map[particle_laddress(lpg, lpg.particles[ipi])];
	    if (ici != tci) move_particle(lpg, ipi, ici, tci);
	}
    }
    verify_lp_grid(lpg);
    return 0;
}

int delete_lp_grid(lp_grid lpg)
{
    free(lpg.map);
    free(lpg.anchors);
    free(lpg.particles);
    return 0;
}

cell get_cell(lp_grid lpg, long i, long j, long k)
{
    cell c;
    long ci = lpg.map[linearize_address(lpg, i, j, k)];
    if (ci == lpg.cell_count) 
	c.start = c.end = &lpg.particles[lpg.anchors[ci]];
    else {
	c.start = &lpg.particles[lpg.anchors[ci]];
	c.end = &lpg.particles[lpg.anchors[ci+1]];
    }
    return c;
}
