#ifndef TERRAIN_H
#define TERRAIN_H

#include <vector>
#include <btBulletDynamicsCommon.h>

struct aabb {
    btVector3 min;
    btVector3 max;
};

struct index3 {
    int i, j, k;
};

// Read .obj file
aabb read_obj(const char* filename, std::vector<btVector3>& vertices, std::vector<btVector3>& faces);

// Convert geometry from .obj to a triangular mesh
btTriangleMesh* import_obj(const char* filename, aabb& aabb);

#endif
