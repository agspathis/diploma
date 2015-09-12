#include "utilities.h"
#include "terrain.h"

terrain read_obj(const char* filename)
{
    terrain t;
    float x_min = +FLT_MAX, y_min = +FLT_MAX, z_min = +FLT_MAX;
    float x_max = -FLT_MAX, y_max = -FLT_MAX, z_max = -FLT_MAX;
    float vx, vy, vz;		// vertex coordinates
    long fi, fj, fk;		// face vertex indices
    long vertex_count=0, face_count=0;
    char line[128];
    FILE* obj;

    if (!(obj = fopen(filename, "rt"))) {
	printf("file not found\n");
	t.taabb.min = btVector3(0, 0, 0);
	t.taabb.max = btVector3(0, 0, 0);
	return t;
    }
    // count vertices, faces
    while (fgets(line, 128, obj)) {
	switch (line[0]) {
	case 'v': vertex_count++; break;
	case 'f': face_count++; break;
	default: continue;
	};
    }
    t.vertices = (vertex*) malloc(vertex_count * sizeof(vertex));
    t.faces = (face*) malloc(face_count * sizeof(face));
    t.vertex_count = vertex_count;
    t.face_count = face_count;

    rewind(obj);
    vertex_count = face_count = 0;
    while (fgets(line, 128, obj)) {
	switch (line[0]) {
	case 'v':
	    sscanf(&line[1],"%f %f %f", &vx, &vy, &vz);
	    if (vx < x_min) x_min = vx;
	    if (vy < y_min) y_min = vy;
	    if (vz < z_min) z_min = vz;
	    if (vx > x_max) x_max = vx;
	    if (vy > y_max) y_max = vy;
	    if (vz > z_max) z_max = vz;
	    t.vertices[vertex_count].x = vx;
	    t.vertices[vertex_count].y = vy;
	    t.vertices[vertex_count].z = vz;
	    vertex_count++;
	    break;
	case 'f':
	    sscanf(&line[1],"%ld %ld %ld", &fi, &fj, &fk);
	    t.faces[face_count].v0i = fi-1;
	    t.faces[face_count].v1i = fj-1;
	    t.faces[face_count].v2i = fk-1;
	    face_count++;
	    break;
	default:
	    continue;
	};
    }
    fclose(obj);

    t.taabb.min = btVector3(x_min, y_min, z_min);
    t.taabb.max = btVector3(x_max, y_max, z_max);

    return t;
}

// reposition terrain such that T.TAABB.MIN = (0, 0, 0) and scale by FACTOR
void dock_scale_terrain(terrain* tp, float factor)
{
    btVector3 min = tp->taabb.min;
    float x_min = min.getX();
    float y_min = min.getY();
    float z_min = min.getZ();
    for (vertex* vp=tp->vertices; vp<tp->vertices + tp->vertex_count; vp++) {
	vp->x = (vp->x - x_min) * factor;
	vp->y = (vp->y - y_min) * factor;
	vp->z = (vp->z - z_min) * factor;
    }

    tp->taabb.min = btVector3(0, 0, 0);
    tp->taabb.max -= min;
    tp->taabb.max *= factor;
}

terrain make_terrain_obj(const char* filename, float scale_factor)
{
    terrain t = read_obj(filename);
    dock_scale_terrain(&t, scale_factor);
    t.triangle_mesh = new btTriangleMesh();
    for(face* fp=t.faces; fp<t.faces + t.face_count; fp++)
    	t.triangle_mesh->addTriangle(btVector3(t.vertices[fp->v0i].x,
					       t.vertices[fp->v0i].y,
					       t.vertices[fp->v0i].z),
				     btVector3(t.vertices[fp->v1i].x,
					       t.vertices[fp->v1i].y,
					       t.vertices[fp->v1i].z),
				     btVector3(t.vertices[fp->v2i].x,
					       t.vertices[fp->v2i].y,
					       t.vertices[fp->v2i].z));

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
    printf("vertex_count = %lu\n", t.vertex_count);
    printf("face_count = %lu\n", t.face_count);
    printf("\n");

    return t;
}

void delete_terrain(terrain t)
{
    free(t.vertices);
    free(t.faces);
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
