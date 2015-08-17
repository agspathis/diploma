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

lp_grid make_lp_grid (btVector3 origin, float step, int x, int y, int z, std::vector<particle> particles)
{
    lp_grid lpg; int i,j,k;
    lpg.x = x; lpg.y = y; lpg.z = z;
    lpg.origin = origin;
    lpg.step = step;

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

    lpg.linear_map = (int*) malloc(x*y*z*sizeof(unsigned int));
    lpg.anchors = (int*) malloc(x*y*z*sizeof(unsigned int));
    lpg.linear_grid = (particle**) malloc (particles.size() * sizeof(particle*));
    Point p;
    for (int pid=0; pid<points.size(); pid++) {
	p = points[pid]; i = p.x(); j=p.y(); k=p.z();
    }
    
    return lpg;
}

// std::vector<particle*> grid_ref (grid g, int i, int j, int k)
// {
//     if ((i<0) or (j<0) or (k<0) or (i>g.x) or (j>g.y) or (k>g.z))
// 	return std::vector<particle*>();
//     return *g.linear_map[linear_id(g, i, j, k)];
// }

int populate_grid (lp_grid &lpg, std::vector<particle> particles)
{
    int i, j, k;
    btVector3 position;

    for (int id=0; id<particles.size(); id++) {
	position = particle_position(particles[id]);
	i = floor(position.getX()/lpg.step);
	j = floor(position.getY()/lpg.step);
	k = floor(position.getZ()/lpg.step);
    }
    return 0;
}

// int update_grid (grid &g)
// {
//     for (int cd=0; g.linear_grid.size(); cd++) {
//     }
//     return 0;
// }
