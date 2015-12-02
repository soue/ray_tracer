#include <cmath>
#include <iostream>
#include "bounding_volume.h"


// BoundingBox::BoundingBox( Point3D p_min, Point3D p_max ) {
// 	_min = p_min;
// 	_max = p_max;
// }

bool BoundingBox::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// Vector3D n = Vector3D(0, 0, 1);
	
	std::cout << "bounding box intersection\n";
	// std::cout << "min point=" << _min << ", max point=" << _max << "\n";
	
	// double xmin = _min[0];
// 	double ymin = _min[1];
// 	double zmin = _min[2];
// 	double xmax = _max[0];
// 	double ymax = _max[1];
// 	double zmax = _max[2];
	
	double t_value;
	
	// Transform ray into object space
	Point3D objOrig = worldToModel * ray.origin;
	Vector3D objDir = worldToModel * ray.dir;
	
	// Find t_value for intersection
	//
	// Equation of plane that square lies on:
	//   a(x-x0) + b(y-y0) + c(z-z0) = 0
	// Equation of the line of intersection:
	//   line(t) = ray_origin + t * ray_dir
	// At intersection, z=0
	//   t = (z - ray.origin_z) / ray.dir_z
	//   t = -ray.origin_z / ray.dir_z
	if (objDir[2] != 0) {
		t_value = -1 * objOrig[2] / objDir[2];
	} else {
		// Square lies on the same plane
		return false;
	}

	// Calculate intersection point
	// double intersectionX = objOrig[0] + t_value * objDir[0];
// 	double intersectionY = objOrig[1] + t_value * objDir[1];
// 	double intersectionZ = objOrig[2] + t_value * objDir[2];
	Point3D intersection = objOrig + t_value * objDir;
	
	// Check if intersection lies within the bounds of the square
	// if (intersection[0] >= xmin && intersection[0] <= xmax && 
// 		intersection[1] >= ymin && intersection[1] <= ymax &&
// 		intersection[2] >= zmin && intersection[2] <= zmax) {
// 		return true;
// 	}
	
	std::cout << "..ending intersect bbox\n";
	
	return false;
}