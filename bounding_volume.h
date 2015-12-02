/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		classes defining primitives in the scene

***********************************************************/

#include "util.h"


class BoundingVolume {
public:
	BoundingVolume() {}
	~BoundingVolume() {}
	virtual bool intersect( Ray3D&, const Matrix4x4&, const Matrix4x4& ) = 0;
	virtual Point3D get_min() = 0;
	virtual Point3D get_max() = 0;
};

class BoundingBox : public BoundingVolume {
public:
	BoundingBox( Point3D p_min, Point3D p_max ) : _min(p_min), _max(p_max) {}
	~BoundingBox() { return; }
	bool intersect( Ray3D&, const Matrix4x4&, const Matrix4x4& );
	Point3D get_min() {return _min;}
	Point3D get_max() {return _max;}
private:
	Point3D _min;
	Point3D _max;
};