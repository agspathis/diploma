#include <iostream>
#include <cfloat>

#include "terrain.h"

// Read .obj file
aabb read_obj(const char* filename, std::vector<btVector3>& vertices, std::vector<btVector3>& faces)
{
    aabb aabb;
    float x_min = +FLT_MAX;
    float y_min = +FLT_MAX;
    float z_min = +FLT_MAX;
    float x_max = -FLT_MAX;
    float y_max = -FLT_MAX;
    float z_max = -FLT_MAX;
    float x, y, z;
    int i, j, k;
    char line[128];
    FILE* objfile;
    if (!(objfile = fopen(filename, "rt"))) {
	std::cout << "file not found\n";
	aabb.min = btVector3(0, 0, 0);
	aabb.max = btVector3(0, 0, 0);
	return aabb;
    }
    while (fgets(line, 128, objfile)) {
	switch (line[0]) {
	case 'v':
	    sscanf(&line[1],"%f %f %f", &x, &y, &z);
	    vertices.push_back(btVector3(x, y, z));
	    if (x < x_min) x_min = x;
	    if (y < y_min) y_min = y;
	    if (z < z_min) z_min = z;
	    if (x > x_max) x_max = x;
	    if (y > y_max) y_max = y;
	    if (z > z_max) z_max = z;
	    break;
	case 'f':
	    sscanf(&line[1],"%d %d %d", &i, &j, &k);
	    faces.push_back(btVector3(--i, --j, --k));
	    break;
	default:
	    continue;
	};
    }
    fclose(objfile);

    aabb.min = btVector3(x_min, y_min, z_min);
    aabb.max = btVector3(x_max, y_max, z_max);

    return aabb;
}

// Convert geometry from .obj to a triangular mesh
btTriangleMesh* import_obj(const char* filename, aabb& aabb)
{
    std::vector<btVector3> vertices;
    std::vector<btVector3> faces;
    aabb = read_obj(filename, vertices, faces);
    btTriangleMesh* triangleMesh = new btTriangleMesh();
    btVector3 f;
    std::cout << faces.size() << "\n";
    for(int i=0; i<faces.size(); i++) {
	f = faces[i];
	triangleMesh->addTriangle(vertices[f.getX()], vertices[f.getY()], vertices[f.getZ()]);
    }
    return triangleMesh;
}
