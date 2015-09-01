#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/spatial_sort.h>

#include "lp_grid.h"

using CGAL::Cartesian;
using CGAL::Point_3;

typedef Cartesian<long> Kernel;
typedef Point_3<Kernel> Point;

long linearize (lp_grid lpg, long i, long j, long k)
{
    // if any index is out of bounds in any of the 3 dimensions of the grid the
    // linear address (LPG.CELL_COUNT) of the last cell containing the off-grid
    // particles is returned
    if (i<0 || i>=lpg.x || j<0 || j>=lpg.y || k<0 || k>=lpg.z)
	return lpg.cell_count;
    else return (i*lpg.y*lpg.z + j*lpg.z + k);
}

anchor* particle_anchor (lp_grid lpg, particle* pp)
{
    btVector3 relative_position = particle_position(pp) - lpg.origin;
    long i = floor(relative_position.getX()/lpg.step);
    long j = floor(relative_position.getY()/lpg.step);
    long k = floor(relative_position.getZ()/lpg.step);
    return lpg.map[linearize(lpg, i, j, k)];
}

// Letters in abbreviations "ia", "tp" etc. throughout this file stand for
// i=initial, t=terminal/target, a=anchor pointer, p=particle pointer
void insert_particle(lp_grid lpg, particle* pp)
{
    anchor* ta = particle_anchor(lpg, pp);
    anchor tp = *(ta+1);
    for (anchor p=lpg.particles+lpg.particle_count; p>tp; p--)
	*p = *(p-1);
    *tp = pp;
    for (anchor* a=ta+1; a<=lpg.anchors+lpg.cell_count; a++) 
	(*a)++;
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
    lpg->map	    = (anchor**) malloc((lpg->cell_count+1) * sizeof(anchor*));
    lpg->anchors    = (anchor*) malloc((lpg->cell_count+1) * sizeof(anchor));
    lpg->particles  = (anchor) malloc ((lpg->particle_count+1) * sizeof(particle*));

    return 0;
}

void verify_lp_grid(lp_grid lpg)
{
    for (anchor* ia=lpg.anchors; ia<lpg.anchors+lpg.cell_count; ia++)
	for (anchor ip=*ia; ip<*(ia+1); ip++) {
	    anchor* ta = particle_anchor(lpg, *ip);
	    if (ia != ta) printf("ERROR pid=%lu ia=%p ip=%p ta=%p\n",
				 (*ip)->id, ia, ip, ta);
	}
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
    std::vector<Point> ps;
    for (i=0; i<lpg.x; i++) {
	for (j=0; j<lpg.y; j++) {
	    for (k=0; k<lpg.z; k++) {
		ps.push_back(Point(i, j, k));
	    }
	}
    }
    CGAL::hilbert_sort(ps.begin(), ps.end());

    // Initialize MAP
    for (long ci=0; ci<ps.size(); ci++) 
	lpg.map[linearize(lpg, ps[ci].x(), ps[ci].y(), ps[ci].z())] = lpg.anchors+ci;
    lpg.map[lpg.cell_count] = lpg.anchors+lpg.cell_count;
    // Initialize ANCHORS to LPG.PARTICLES (start of array)
    for (anchor* a=lpg.anchors; a<=lpg.anchors+lpg.cell_count; a++) *a = lpg.particles;
    // Populate PARTICLES
    for (particle* pp=fluid.particles; pp<fluid.particles+lpg.particle_count; pp++)
	insert_particle(lpg, pp);

    verify_lp_grid(lpg);
    return lpg;
}

void update_lp_grid (lp_grid lpg)
{
    anchor tp;
    particle* tmp_storage;
    for (anchor* ia=lpg.anchors; ia<lpg.anchors+lpg.cell_count; ia++)
	for (anchor ip=*ia; ip<*(ia+1); ip++) {
	    anchor* ta = particle_anchor(lpg, *ip);
	    if (ia != ta) {
		tmp_storage = *ip;
		if (ia < ta) {
		    tp = *ta-1;
		    for (anchor p=ip; p<tp; p++) *p = *(p+1);
		    for (anchor* a=ia+1; a<=ta; a++) (*a)--;
		    *tp = tmp_storage;
		    ip--;	// correction due to shift
		}
		else {
		    tp = *(ta+1);
		    for (anchor p=ip; p>tp; p--) *p = *(p-1);
		    for (anchor* a=ia; a>ta; a--) (*a)++;
		    *tp = tmp_storage;
		}
	    }
	}
    verify_lp_grid(lpg);
}

void delete_lp_grid(lp_grid lpg)
{
    free(lpg.map);
    free(lpg.anchors);
    free(lpg.particles);
}

cell get_cell(lp_grid lpg, long i, long j, long k)
{
    cell c;
    anchor* a = lpg.map[linearize(lpg, i, j, k)];
    if (a == lpg.map[lpg.cell_count]) {
	c.start = *a;
	c.end = &lpg.particles[lpg.particle_count];
    }
    else {
	c.start = *a;
	c.end = *(a+1);
    }
    return c;
}
