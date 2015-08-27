#include "utilities.h"
#include "terrain.h"

aabb read_obj(char* filename, std::vector<btVector3> &vertices, std::vector<btVector3> &faces)
{
    aabb aabb;
    float x_min = +FLT_MAX;
    float y_min = +FLT_MAX;
    float z_min = +FLT_MAX;
    float x_max = -FLT_MAX;
    float y_max = -FLT_MAX;
    float z_max = -FLT_MAX;
    float x, y, z;
    long i, j, k;
    char line[128];
    FILE* objfile;
    if (!(objfile = fopen(filename, "rt"))) {
	printf("file not found\n");
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

// Import model from .obj file to a TERRAIN structure
terrain import_obj(char* filename)
{
    terrain t;
    std::vector<btVector3> vertices;
    std::vector<btVector3> faces;
    t.terrain_aabb = read_obj(filename, vertices, faces);
    t.vertex_count = vertices.size();
    t.face_count = faces.size();
    t.triangle_mesh = new btTriangleMesh();
    for(long fi=0; fi<faces.size(); fi++)
	t.triangle_mesh->addTriangle(vertices[faces[fi].getX()],
				     vertices[faces[fi].getY()],
				     vertices[faces[fi].getZ()]);

    // construct rigid body for simulation
    t.shape = new btBvhTriangleMeshShape(t.triangle_mesh,true);
    t.motion_state = new btDefaultMotionState(btTransform(btQuaternion(0, 0, 0, 1),
							  btVector3(0, 0, 0)));
    btRigidBody::btRigidBodyConstructionInfo t_ci(0, t.motion_state, t.shape, btVector3(0, 0, 0));
    t_ci.m_restitution = 1.0;
    t.rigid_body = new btRigidBody(t_ci);

    // print terrain aabb
    printf("Terrain AABB info:\n");
    print_aabb(t.terrain_aabb);
    printf("\n");

    return t;
}

int delete_terrain(terrain t)
{
    delete t.shape;
    delete t.triangle_mesh;
    delete t.motion_state;
    delete t.rigid_body;
    return 0;
}
