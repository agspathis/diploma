#ifndef TERRAIN_H
#define TERRAIN_H

#include <vector>
#include <btBulletDynamicsCommon.h>

struct aabb {
    btVector3 min;
    btVector3 max;
};

struct terrain {
    long vertex_count;
    long face_count;
    aabb terrain_aabb;
    btTriangleMesh* triangle_mesh;
    btCollisionShape* shape;
    btDefaultMotionState* motion_state;
    btRigidBody* rigid_body;
};

terrain import_obj(char* filename);

int delete_terrain(terrain t);

#endif
