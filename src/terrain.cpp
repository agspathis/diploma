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
    btCollisionShape* t_shape = new btBvhTriangleMeshShape(t.triangle_mesh, true);
    btDefaultMotionState* t_motion_state = new btDefaultMotionState
	(btTransform(btQuaternion(0, 0, 0, 1), btVector3(0, 0, 0)));
    btRigidBody::btRigidBodyConstructionInfo t_ci
	(0, t_motion_state, t_shape, btVector3(0, 0, 0));
    t_ci.m_restitution = 0.0;
    t_ci.m_friction = 0.0;
    t_ci.m_rollingFriction = 0.0;
    t.rigid_body = new btRigidBody(t_ci);

    // print terrain aabb
    printf("TERRAIN INFO:\n");
    print_aabb(t.taabb);
    printf("vertex_count = %lu\n", t.vertex_count);
    printf("face_count = %lu\n", t.face_count);
    printf("\n");

    return t;
}

// creates the aabb boundary of terrain T as a separate terrain object
terrain terrain_boundary(terrain t)
{
    terrain tb;
    tb.taabb = t.taabb;
    float x_max = t.taabb.max.getX();
    float y_max = t.taabb.max.getY();
    float z_max = t.taabb.max.getZ();
    tb.vertex_count = 8;
    tb.face_count = 12;
    tb.vertices = (vertex*) malloc(tb.vertex_count * sizeof(vertex));
    tb.faces = (face*) malloc(tb.face_count * sizeof(face));

    tb.vertices[0].x = 0;     tb.vertices[0].y = 0;     tb.vertices[0].z = z_max;
    tb.vertices[1].x = 0;     tb.vertices[1].y = 0;     tb.vertices[1].z = 0;
    tb.vertices[2].x = 0;     tb.vertices[2].y = y_max; tb.vertices[2].z = 0;
    tb.vertices[3].x = 0;     tb.vertices[3].y = y_max; tb.vertices[3].z = z_max;
    tb.vertices[4].x = x_max; tb.vertices[4].y = 0;     tb.vertices[4].z = z_max;
    tb.vertices[5].x = x_max; tb.vertices[5].y = 0;     tb.vertices[5].z = 0;
    tb.vertices[6].x = x_max; tb.vertices[6].y = y_max; tb.vertices[6].z = 0;
    tb.vertices[7].x = x_max; tb.vertices[7].y = y_max; tb.vertices[7].z = z_max;

    // hardcoded topology
    tb.faces[ 0].v0i = 3; tb.faces[ 0].v1i = 2; tb.faces[ 0].v2i = 1;
    tb.faces[ 1].v0i = 1; tb.faces[ 1].v1i = 0; tb.faces[ 1].v2i = 3;
    tb.faces[ 2].v0i = 1; tb.faces[ 2].v1i = 5; tb.faces[ 2].v2i = 4;
    tb.faces[ 3].v0i = 4; tb.faces[ 3].v1i = 0; tb.faces[ 3].v2i = 1;
    tb.faces[ 4].v0i = 2; tb.faces[ 4].v1i = 6; tb.faces[ 4].v2i = 5;
    tb.faces[ 5].v0i = 5; tb.faces[ 5].v1i = 1; tb.faces[ 5].v2i = 2;
    tb.faces[ 6].v0i = 7; tb.faces[ 6].v1i = 6; tb.faces[ 6].v2i = 2;
    tb.faces[ 7].v0i = 2; tb.faces[ 7].v1i = 3; tb.faces[ 7].v2i = 7;
    tb.faces[ 8].v0i = 4; tb.faces[ 8].v1i = 7; tb.faces[ 8].v2i = 3;
    tb.faces[ 9].v0i = 3; tb.faces[ 9].v1i = 0; tb.faces[ 9].v2i = 4;
    tb.faces[10].v0i = 5; tb.faces[10].v1i = 6; tb.faces[10].v2i = 7;
    tb.faces[11].v0i = 7; tb.faces[11].v1i = 4; tb.faces[11].v2i = 5;

    tb.triangle_mesh = new btTriangleMesh();
    for(face* fp=tb.faces; fp<tb.faces + tb.face_count; fp++)
	tb.triangle_mesh->addTriangle(btVector3(tb.vertices[fp->v0i].x,
						tb.vertices[fp->v0i].y,
						tb.vertices[fp->v0i].z),
				      btVector3(tb.vertices[fp->v1i].x,
						tb.vertices[fp->v1i].y,
						tb.vertices[fp->v1i].z),
				      btVector3(tb.vertices[fp->v2i].x,
						tb.vertices[fp->v2i].y,
						tb.vertices[fp->v2i].z));

    // construct rigid body for simulation
    btCollisionShape* tb_shape = new btBvhTriangleMeshShape(tb.triangle_mesh, true);
    btDefaultMotionState* tb_motion_state = new btDefaultMotionState
	(btTransform(btQuaternion(0, 0, 0, 1), btVector3(0, 0, 0)));
    btRigidBody::btRigidBodyConstructionInfo tb_ci
	(0, tb_motion_state, tb_shape, btVector3(0, 0, 0));
    tb_ci.m_restitution = 0.0;
    tb_ci.m_friction = 0.0;
    tb_ci.m_rollingFriction = 0.0;
    tb.rigid_body = new btRigidBody(tb_ci);

    return tb;
}

void delete_terrain(terrain t)
{
    free(t.vertices);
    free(t.faces);
    delete t.triangle_mesh;
    delete t.rigid_body->getCollisionShape();
    delete t.rigid_body->getMotionState();
    delete t.rigid_body;
}

float aabb_volume(aabb aabb)
{
    aabb.max -= aabb.min;
    return (aabb.max.getX() * aabb.max.getY() * aabb.max.getZ());
}

void collect_terrain_impulses(btDynamicsWorld* dynamics_world,
			      std::vector<terrain_impulse>& terrain_impulses)
{
    long manifold_count = dynamics_world->getDispatcher()->getNumManifolds();
    for (long mi=0; mi<manifold_count; mi++) {
	btPersistentManifold* cm = dynamics_world->getDispatcher()->getManifoldByIndexInternal(mi);
	const btCollisionObject* co0 = static_cast<const btCollisionObject*>(cm->getBody0());
	const btCollisionObject* co1 = static_cast<const btCollisionObject*>(cm->getBody1());

	// filter fluid/terrain collisions
	if ((co0->getCollisionShape()->getShapeType() == TRIANGLE_MESH_SHAPE_PROXYTYPE and
	     co1->getCollisionShape()->getShapeType() == SPHERE_SHAPE_PROXYTYPE)
	    or
	    (co0->getCollisionShape()->getShapeType() == SPHERE_SHAPE_PROXYTYPE and
	     co1->getCollisionShape()->getShapeType() == TRIANGLE_MESH_SHAPE_PROXYTYPE)) {

	    for (int cci=0; cci<cm->getNumContacts(); cci++)
	    {
		btManifoldPoint& cpt = cm->getContactPoint(cci);
		const btVector3& normal = cpt.m_normalWorldOnB;
		float impulse = cpt.getAppliedImpulse();
		if (cpt.getAppliedImpulse() != 0) {
		    terrain_impulse ti;
		    ti.position = cpt.getPositionWorldOnB();
		    ti.impulse = impulse;
		    terrain_impulses.push_back(ti);
		}
	    }
	}
    }
}
