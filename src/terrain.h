#ifndef TERRAIN_H
#define TERRAIN_H

#include <vector>
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
    btTriangleMesh* triangle_mesh;
    btRigidBody* rigid_body;
};

struct terrain_impulse {
    btVector3 position;
    float impulse;
};

float aabb_volume(aabb aabb);

terrain make_terrain_obj(const char* filename, float scale_factor);

terrain terrain_boundary(terrain t);

void delete_terrain(terrain t);

void collect_terrain_impulses(btDynamicsWorld* dynamics_world,
			      std::vector<terrain_impulse>& terrain_impulses);

#endif
