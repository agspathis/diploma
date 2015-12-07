#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/spatial_sort.h>

#include "lp_grid.h"

using CGAL::Cartesian;
using CGAL::Point_3;

typedef Cartesian<long> Kernel;
typedef Point_3<Kernel> Point;

// Function LINEARIZE implements standard x-major linear indexing (following
// vtk's choice of x-major indexing in structured_points dataset format). If any
// index is out of bounds in any of the 3 dimensions of the grid the linear
// address of the last cell (LPG.CELL_COUNT) containing the off-grid particles
// is returned. SDF denotes the subdivision factor of the addressed grid.
long linearize(lp_grid lpg, long i, long j, long k, int sdf)
{
    if (i<0 || i>=sdf*lpg.x || j<0 || j>=sdf*lpg.y || k<0 || k>=sdf*lpg.z)
	return sdf * lpg.cell_count;
    else return i + j*sdf*lpg.x + k*pow(sdf, 2)*lpg.x*lpg.y;
}

anchor* particle_anchor(lp_grid lpg, particle* pp)
{
    btVector3 relative_position = particle_position(pp) - lpg.origin;
    long i = floor(relative_position.getX()/lpg.step);
    long j = floor(relative_position.getY()/lpg.step);
    long k = floor(relative_position.getZ()/lpg.step);
    return lpg.map[linearize(lpg, i, j, k, 1)];
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
    lpg->cf_sdf = DEFAULT_CF_SDF;

    // allocate memory for the 3 arrays
    lpg->map = (anchor**) malloc((lpg->cell_count+1) * sizeof(anchor*));
    lpg->anchors = (anchor*) malloc((lpg->cell_count+1) * sizeof(anchor));
    lpg->particles = (anchor) malloc ((lpg->particle_count+1) * sizeof(particle*));
    lpg->color_field = (float*) malloc ((lpg->cell_count)
					* pow(lpg->cf_sdf, 3)
					* sizeof(float));

    return 0;
}

// Henceforth, letters 'i', 't' in variable names stand for 'initial' and
// 'terminal/target' respectively
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
    // grid parameter initialization/allocation
    lp_grid lpg; long i,j,k;
    allocate_lp_grid(&lpg, domain, fluid);
    printf("LP_GRID INFO:\n");
    printf("x=%lu, y=%lu, z=%lu\n", lpg.x, lpg.y, lpg.z);
    printf("xss=%lu, yss=%lu, zss=%lu\n", lpg.xss, lpg.yss, lpg.zss);
    printf("cell_count=%lu\n", lpg.cell_count);
    printf("\n");

    // make points in centers of cells and spatially sort them
    std::vector<Point> ps;
    for (i=0; i<lpg.x; i++) {
	for (j=0; j<lpg.y; j++) {
	    for (k=0; k<lpg.z; k++) {
		ps.push_back(Point(i, j, k));
	    }
	}
    }
    CGAL::hilbert_sort(ps.begin(), ps.end());

    // initialize MAP according to spatial sort
    for (long ci=0; ci<ps.size(); ci++)
	lpg.map[linearize(lpg, ps[ci].x(), ps[ci].y(), ps[ci].z(), 1)]
	    = lpg.anchors+ci;
    lpg.map[lpg.cell_count] = lpg.anchors+lpg.cell_count;
    // initialize array holding the particle count for each cell.
    ptrdiff_t cell_pc[lpg.cell_count+1] = {0};
    // scan particles and sum up the particle count for each cell in CELL_PC
    for (particle* pp=fluid.particles; pp<fluid.particles+lpg.particle_count; pp++)
	cell_pc[particle_anchor(lpg, pp) - lpg.anchors]++;
    // initialize anchors
    size_t anchor_offset=0;
    for (size_t ci = 0; ci<lpg.cell_count+1; ci++) {
	lpg.anchors[ci] = lpg.particles + anchor_offset;
	anchor_offset += cell_pc[ci];
    }
    // Populate particles by filling pointers to particles in the cells where
    // anchors point, and incrementing the respective anchor each time. At the
    // end, each anchor points to the first particle of the next cell (the last
    // one to unallocated memory)
    for (particle* pp=fluid.particles; pp<fluid.particles+lpg.particle_count; pp++) {
	anchor* ta = particle_anchor(lpg, pp);
	**ta = pp;
	(*ta)++;
    }
    // set each anchor to point where the previous one does
    for (size_t ai = lpg.cell_count; ai>0; ai--)
	lpg.anchors[ai] = lpg.anchors[ai-1];
    lpg.anchors[0] = lpg.particles;

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
		    for (anchor* a=ia+1; a<=ta; a++) (*a)--;
		    for (anchor p=ip; p<tp; p++) *p = *(p+1);
		    *tp = tmp_storage;
		    ip--;	// correction due to shift
		}
		else {
		    tp = *(ta+1);
		    for (anchor* a=ia; a>ta; a--) (*a)++;
		    for (anchor p=ip; p>tp; p--) *p = *(p-1);
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
    free(lpg.color_field);
}

cell get_cell(lp_grid lpg, long i, long j, long k)
{
    cell c;
    anchor* a = lpg.map[linearize(lpg, i, j, k, 1)];
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
