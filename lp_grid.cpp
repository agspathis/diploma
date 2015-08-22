#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/spatial_sort.h>

#include "lp_grid.h"
#include "utilities.h"

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

long particle_laddress (lp_grid lpg, particle p)
{
    btVector3 position = particle_position(p) - lpg.origin;
    long i = floor(position.getX()/lpg.step);
    long j = floor(position.getY()/lpg.step);
    long k = floor(position.getZ()/lpg.step);
    return linearize_address(lpg, i, j, k);
}

// Letters in abbreviations "ici", "tpi" etc. stand for
// i=initial, t=terminal, c=cell, p=particle, i=index
int insert_particle(lp_grid lpg, particle* pp)
{
    // terminal cell and particle index
    long tci = lpg.map[particle_laddress(lpg, *pp)];
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
    lpg->map	    = (long*) malloc(lpg->cell_count * sizeof(long));
    lpg->anchors    = (long*) malloc((lpg->cell_count+1) * sizeof(long));
    lpg->particles  = (particle**) malloc ((lpg->particle_count+1) * sizeof(particle*));

    return 0;
}

lp_grid make_lp_grid (aabb domain, float step, std::vector<particle> particles)
{
    // Grid parameter initialization/allocation
    lp_grid lpg; long i,j,k;
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
    long cind;			// cind = cell index
    for (cind=0; cind<points.size(); cind++) {
	lpg.map[linearize_address(lpg, points[cind].x(), points[cind].y(), points[cind].z())] = cind;
    }

    // Initialize ANCHORS to 0
    for (cind=0; cind<=lpg.cell_count; cind++) lpg.anchors[cind] = 0;

    // Populate PARTICLES
    for (long pind=0; pind<lpg.particle_count; pind++) // pind = particle index
	insert_particle(lpg, &particles[pind]);

    return lpg;
}


int move_particle(lp_grid lpg, long ipi, long ici, long tci)
{
    long ci, pi, tpi;
    particle* tmp_storage = lpg.particles[ipi];

    if (ici < tci) {
	long tpi = lpg.anchors[tci]-1;
	particle* tmp_storage = lpg.particles[ipi];
	for (pi=ipi+1; pi<=tpi; pi++)
	    lpg.particles[pi-1] = lpg.particles[pi];
	for (ci=ici+1; ci<=tci; ci++) lpg.anchors[ci]--;
    }
    else {
	tpi = lpg.anchors[tci+1];
	particle* tmp_storage = lpg.particles[ipi];
	for (pi=ipi-1; pi>=tpi; pi--)
	    lpg.particles[pi+1] = lpg.particles[pi];
	for (ci=tci+1; ci<=ici; ci++) lpg.anchors[ci]++;
    }
    lpg.particles[tpi] = tmp_storage;
    return 0;
}

int update_lp_grid (lp_grid lpg)
{
    for (long ici=0; ici<lpg.cell_count; ici++) {
	for (long ipi=lpg.anchors[ici]; ipi<lpg.anchors[ici+1]; ipi++) {
	    // tci is the index of the cell to which the current particle
	    // actully belongs according to its position in space
	    long tci = lpg.map[particle_laddress(lpg, *lpg.particles[ipi])];
	    move_particle(lpg, ipi, ici, tci);
	}
    }
    return 0;
}

particle_range get_cell(lp_grid lpg, long i, long j, long k)
{
    particle_range pr;
    long cind = lpg.map[linearize_address(lpg, i, j, k)];
    pr.start = lpg.particles[lpg.anchors[cind]];
    pr.end = lpg.particles[lpg.anchors[cind+1]];
    return pr;
}
