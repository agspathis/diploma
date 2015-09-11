#include "utilities.h"
#include "terrain.h"

#define TERRAIN_SCALING_FACTOR 0.04

model read_obj(const char* filename)
{
    model m;
    float x_min = +FLT_MAX, y_min = +FLT_MAX, z_min = +FLT_MAX;
    float x_max = -FLT_MAX, y_max = -FLT_MAX, z_max = -FLT_MAX;
    float x, y, z;		// vertex coordinates
    long i, j, k;		// face vertex indices
    long vertex_count=0, face_count=0;
    char line[128];
    FILE* obj_file;

    if (!(obj_file = fopen(filename, "rt"))) {
	printf("file not found\n");
	m.maabb.min = btVector3(0, 0, 0);
	m.maabb.max = btVector3(0, 0, 0);
	return m;
    }
    // count vertices, faces
    while (fgets(line, 128, obj_file)) {
	switch (line[0]) {
	case 'v': vertex_count++; break;
	case 'f': face_count++; break;
	default: continue;
	};
    }
    m.vertices = (vertex*) malloc(vertex_count * sizeof(vertex));
    m.faces = (face*) malloc(face_count * sizeof(face));
    m.vertex_count = vertex_count;
    m.face_count = face_count;

    vertex_count = face_count = 0;
    rewind(obj_file);
    while (fgets(line, 128, obj_file)) {
	switch (line[0]) {
	case 'v':
	    sscanf(&line[1],"%f %f %f", &x, &y, &z);
	    m.vertices[vertex_count++] = btVector3(x, y, z);
	    if (x < x_min) x_min = x;
	    if (y < y_min) y_min = y;
	    if (z < z_min) z_min = z;
	    if (x > x_max) x_max = x;
	    if (y > y_max) y_max = y;
	    if (z > z_max) z_max = z;
	    break;
	case 'f':
	    sscanf(&line[1],"%d %d %d", &i, &j, &k);
	    m.faces[face_count].v0i = i-1;
	    m.faces[face_count].v1i = j-1;
	    m.faces[face_count].v2i = k-1;
	    face_count++;
	    break;
	default:
	    continue;
	};
    }
    fclose(obj_file);

    m.maabb.min = btVector3(x_min, y_min, z_min);
    m.maabb.max = btVector3(x_max, y_max, z_max);

    return m;
}

// reposition the model such that T.TAABB.MIN = (0, 0, 0)
void dock_scale_model(model* m, float factor)
{
    vertex min = m->maabb.min;
    for (vertex* vp=m->vertices; vp<m->vertices + m->vertex_count; vp++) {
	*vp -= min;
	*vp *= factor;
    }
    m->maabb.min -= min;
    m->maabb.min *= factor;
    m->maabb.max -= min;
    m->maabb.max *= factor;
}

void delete_model(model m)
{
    free(m.vertices);
    free(m.faces);
}

terrain make_terrain_obj(const char* filename)
{
    terrain t;
    model m = read_obj(filename);
    dock_scale_model(&m, TERRAIN_SCALING_FACTOR);
    t.taabb = m.maabb;
    t.triangle_mesh = new btTriangleMesh();
    for(long fi=0; fi<m.face_count; fi++) {
    	face f = m.faces[fi];
    	t.triangle_mesh->addTriangle(m.vertices[f.v0i],
    				     m.vertices[f.v1i],
    				     m.vertices[f.v2i]);
    }

    // construct rigid body for simulation
    t.shape = new btBvhTriangleMeshShape(t.triangle_mesh, true);
    btDefaultMotionState* t_motion_state = new btDefaultMotionState
	(btTransform(btQuaternion(0, 0, 0, 1), btVector3(0, 0, 0)));
    btRigidBody::btRigidBodyConstructionInfo t_ci
	(0, t_motion_state, t.shape, btVector3(0, 0, 0));
    // t_ci.m_restitution = 0.8;	// fluid adhesion to terrain
    t.rigid_body = new btRigidBody(t_ci);

    // print terrain aabb
    printf("TERRAIN INFO:\n");
    print_aabb(t.taabb);
    printf("vertex_count = %lu\n", m.vertex_count);
    printf("face_count = %lu\n", m.face_count);
    printf("\n");

    delete_model(m);
    return t;
}

void delete_terrain(terrain t)
{
    delete t.shape;
    delete t.triangle_mesh;
    delete t.rigid_body->getMotionState();
    delete t.rigid_body;
}

float aabb_volume(aabb aabb)
{
    aabb.max -= aabb.min;
    return (aabb.max.getX() * aabb.max.getY() * aabb.max.getZ());
}
