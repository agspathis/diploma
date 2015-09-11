#ifndef TERRAIN_H
#define TERRAIN_H

#include <btBulletDynamicsCommon.h>

typedef btVector3 vertex;

struct aabb {
    btVector3 min;
    btVector3 max;
};

struct face {
    long v0i;
    long v1i;
    long v2i;
};

struct model {
    long vertex_count;
    long face_count;
    vertex* vertices;
    face* faces;
    aabb maabb;
};

struct terrain {
    long vertex_count;
    long face_count;
    aabb taabb;
    btTriangleMesh* triangle_mesh;
    btCollisionShape* shape;
    btRigidBody* rigid_body;
};

float aabb_volume(aabb aabb);

terrain make_terrain_obj(const char* filename);

void delete_terrain(terrain t);

#endif
