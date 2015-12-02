/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		classes defining primitives in the scene

***********************************************************/

#include "util.h"
#include "mesh.h"



// All primitives should provide a intersection function.
// To create more primitives, inherit from SceneObject.
// Namely, you can create, Sphere, Cylinder, etc... classes
// here.
class SceneObject {
public:
	virtual int get_type() = 0;
	// Returns true if an intersection occured, false otherwise.
	virtual bool intersect( Ray3D&, const Matrix4x4&, const Matrix4x4& ) = 0;
	virtual Texture* get_texture() = 0;
};

// Example primitive you can create, this is a unit square on
// the xy-plane.
class UnitSquare : public SceneObject {
public:
	UnitSquare( Texture* texture ) : _texture(texture) {}
	int get_type() { return SQUARE; }
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
	Texture* get_texture() { return _texture; }

private:
	Texture* _texture;
};

class UnitSphere : public SceneObject {
public:
	UnitSphere( Texture* texture ) : _texture(texture) {}
	int get_type() { return SPHERE; }
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
	Texture* get_texture() { return _texture; }

private:
	Texture* _texture;
};

class UnitCylinder : public SceneObject {
public:
	UnitCylinder( Texture* texture ) : _texture(texture) {}
	int get_type() { return CYLINDER; }
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
	bool intersectCylinderCap( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, double t_value );
	Texture* get_texture() { return _texture; }

private:
	Texture* _texture;
};

class Triangle : public SceneObject {
public:
	Triangle( Point3D v0, Point3D v1, Point3D v2, Vector3D normal,
		Texture* texture ) : _v0(v0), _v1(v1), _v2(v2), _normal(normal), _texture(texture) {}
	int get_type() { return TRIANGLE; }
	bool intersect(Ray3D& r, Ray3D& ray, const Matrix4x4& worldToModel, const Matrix4x4& modelToWorld);
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
	Point3D get_vertex1() { return _v0; }
	Point3D get_vertex2() { return _v1; }
	Point3D get_vertex3() { return _v2; }
	Texture* get_texture() { return _texture; }

private:
	Point3D _v0;
	Point3D _v1;
	Point3D _v2;
	Vector3D _normal;
	Texture* _texture;
};

class TriangleSimple : public SceneObject {
public:
	TriangleSimple() : _triangle(*(new Point3D(0,1,0)),*(new Point3D(1,0,0)),*(new Point3D(-1,0,0)),*(new Vector3D(0,0,1)), NULL) {}
	int get_type() { return _triangle.get_type(); }
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel, const Matrix4x4& modelToWorld ) { return _triangle.intersect(ray,worldToModel,modelToWorld); }
	Texture* get_texture() { return _triangle.get_texture(); }

private:
	Triangle _triangle;
};

class MeshObject : public SceneObject {
public:
	MeshObject(Mesh* mesh, Texture* texture) : _mesh(mesh), _add_to_scene(true), _texture(texture),
		_bounding_vol(mesh->get_bounding_vol()) {}
	int get_type() { return MESH; }
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
	Mesh * get_mesh() { return _mesh; }
	bool add_to_scene() { return _add_to_scene; }
	void set_add_to_scene(bool add) { _add_to_scene = add; }
	unsigned char* get_rbuffer() { return NULL; }
	unsigned char* get_gbuffer() { return NULL; }
	unsigned char* get_bbuffer() { return NULL; }
	void set_rbuffer(unsigned char* rbuffer) { return; }
	void set_gbuffer(unsigned char* gbuffer) { return; }
	void set_bbuffer(unsigned char* bbuffer) { return; }
	Texture* get_texture() { return _texture; }
private:
	Mesh* _mesh;
	bool _add_to_scene;
	BoundingVolume* _bounding_vol;
	Texture* _texture;
};
