#ifndef TERRAIN_H
#define TERRAIN_H

#include <stdint.h>
#include <btBulletDynamicsCommon.h>

// typedef btVector3 vertex;

struct aabb {
    btVector3 min;
    btVector3 max;
};

struct vertex {
    float x;
    float y;
    float z;
};

struct face {
    long v0i;
    long v1i;
    long v2i;
};

struct terrain {
    aabb taabb;
    long vertex_count;
    long face_count;
    vertex* vertices;
    face* faces;
    float* forces;
    btTriangleMesh* triangle_mesh;
    btCollisionShape* shape;
    btRigidBody* rigid_body;
};

float aabb_volume(aabb aabb);

terrain make_terrain_obj(const char* filename, float scale_factor);

void delete_terrain(terrain t);

void collect_terrain_forces(btDynamicsWorld* dynamics_world, terrain t);

#endif
