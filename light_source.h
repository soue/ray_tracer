/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		   light source classes

***********************************************************/

#include "util.h"
// #include "scene_object.h"

// Constants
//
// Shading mode
enum { PHONG, AMBIENT_DIFFUSE, SIGNATURE };
const int SHADE_MODE = PHONG;
const bool SHADOWS_ENABLED = false;
const int SHADOW_RAYS = 1;



// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray
// differently.
class LightSource {
public:
	virtual void shade( Ray3D& ) = 0;
	virtual Point3D get_position() = 0;
	virtual Point3D get_sample_position() = 0; // sample different points from the area light
	virtual Colour get_amb() const = 0;

};


// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col),
		_col_diffuse(col), _col_specular(col) {}
	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) :
		_pos(pos), _col_ambient(ambient), _col_diffuse(diffuse),
		_col_specular(specular) {}
	void shade( Ray3D& ray );
	Colour get_amb() const {return _col_ambient;}

	Point3D get_position() { return _pos; }
	Point3D get_sample_position() { return _pos; }

private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse;
	Colour _col_specular;
};

// A rectangular area (extended) light is defined by the top-left, top-right,
// bottom-left and bottom-right vertices, the normal and colour.
class AreaLight : public LightSource {
public:
	AreaLight( Point3D top_left, Point3D top_right,
	Point3D bottom_left, Point3D bottom_right, Colour col ) :
		_top_left(top_left), _top_right(top_right),
		_bottom_left(bottom_left), _bottom_right(bottom_right),
		_normal(get_light_normal()),
		_col_ambient(col), _col_diffuse(col), _col_specular(col) {}
	AreaLight( Point3D top_left, Point3D top_right,
	Point3D bottom_left, Point3D bottom_right,
	Colour ambient, Colour diffuse, Colour specular ) :
		_top_left(top_left), _top_right(top_right),
		_bottom_left(bottom_left), _bottom_right(bottom_right),
		_normal(get_light_normal()),
		_col_ambient(ambient), _col_diffuse(diffuse), _col_specular(specular) {}
	void shade( Ray3D& ray );
	Colour get_amb() const {return _col_ambient;}

	Point3D get_position() { return get_light_center(); }
	Point3D get_sample_position();
	Vector3D get_light_normal();
	Point3D get_light_center();

private:
	Colour _col_ambient;
	Colour _col_diffuse;
	Colour _col_specular;
	Vector3D _normal;
	Point3D _top_left;
	Point3D _top_right;
	Point3D _bottom_left;
	Point3D _bottom_right;
};



// A rectangular area (extended) light is defined by the top-left, top-right,
// bottom-left and bottom-right vertices, the normal and colour.
class SphereLight : public LightSource {
public:
	SphereLight( Point3D center, double radius, Colour col ) :
		_center(center), _radius(radius),
		_col_ambient(col), _col_diffuse(col), _col_specular(col) {}
	void shade( Ray3D& ray );
	Colour get_amb() const {return _col_ambient;}
	Point3D get_position() { return _center; }
	Point3D get_sample_position();
	Point3D get_light_center();

private:
	Colour _col_ambient;
	Colour _col_diffuse;
	Colour _col_specular;
	Point3D _center;
	double _radius;
};
