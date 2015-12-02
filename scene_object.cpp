/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"
#include "bmp_io.h"


bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0),
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point,
	// intersection.normal, intersection.none, intersection.t_value.
	//
	// HINT: Remember to first transform the ray into object space
	// to simplify the intersection test.
	// Assume unit square defined on xy-plane
	Vector3D n = Vector3D(0, 0, 1);
	double xmin = -0.5;
	double ymin = -0.5;
	double xmax = 0.5;
	double ymax = 0.5;
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
	if (t_value > M_EPSILON) {
		// Calculate intersection point
		double intersectionX = objOrig[0] + t_value * objDir[0];
		double intersectionY = objOrig[1] + t_value * objDir[1];

		// Check if intersection lies within the bounds of the square
		if (intersectionX >= xmin && intersectionX <= xmax &&
			intersectionY >= ymin && intersectionY <= ymax) {
			// Check for nearest intersection
			if (ray.intersection.none == false && t_value >= ray.intersection.t_value) {
				return false;
			}
			// Intersection found (convert back to world coordinates)
			ray.intersection.none = false;
			ray.intersection.point = modelToWorld * Point3D(intersectionX, intersectionY, 0);
			ray.intersection.normal = transNorm(worldToModel, n);
			ray.intersection.normal.normalize();
			ray.intersection.t_value = t_value;
			ray.intersection.obj = this;
			return true;
		}
	}
	return false;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
	const Matrix4x4& modelToWorld ) {

	Point3D objOrig = worldToModel*ray.origin;
	Vector3D objDir = worldToModel*ray.dir;

	double A = objDir.dot(objDir);

	Vector3D objOrigAsVector = Vector3D(objOrig[0],objOrig[1],objOrig[2]);
	double B = objOrigAsVector.dot(objDir);
	double C = objOrigAsVector.dot(objOrigAsVector)-1;

	double D = B*B -A*C;

	if(D < 0) {
		return false;
	}

	double t1 = -B/A + sqrt(D)/A;
	double t2 = -B/A - sqrt(D)/A;

	if(t1 < 0 && t2 < 0){
			return false;
	}

	if(t1 > 0 && t2 < 0){
			Point3D intersection = objOrig + t1*objDir;
			if( (ray.intersection.none || ray.intersection.t_value > t1) && t1>0.00000001){
					ray.intersection.t_value = t1;
					ray.intersection.point = modelToWorld*intersection;
					ray.intersection.normal = transNorm(worldToModel, intersection - Point3D(0,0,0));
					ray.intersection.normal.normalize();
					ray.intersection.none = false;
					ray.intersection.obj = this;
					return true;
			}
	}

	if(t1 > 0 && t2 > 0){
			Point3D intersection = objOrig + t2*objDir;
			if( (ray.intersection.none || ray.intersection.t_value > t2) && t2>0.00000001){
					ray.intersection.t_value = t2;
					ray.intersection.point = modelToWorld*intersection;
					ray.intersection.normal = transNorm(worldToModel, intersection - Point3D(0,0,0));
					ray.intersection.normal.normalize();
					ray.intersection.none = false;
					ray.intersection.obj = this;
					return true;
			}
	}

	return false;

}

// Helper function to determine if there is an intesection with the caps of the cylinder
bool UnitCylinder::intersectCylinderCap( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, double t_value ) {

	Point3D cylinderCenter(0, 0, 0);
	// Vector3D objOrig = worldToModel * ray.origin - cylinderCenter;
	Point3D objOrig = worldToModel * ray.origin;
	Vector3D origAsVector = Vector3D(objOrig[0],objOrig[1],objOrig[2]);
	Vector3D objDir = worldToModel * ray.dir;

	// Check for earlier intersection
	if ((ray.intersection.none || t_value < ray.intersection.t_value) &&
			t_value > M_EPSILON) {

		// Calculate intersection point
		// Point3D capIntersect = worldToModel * ray.origin + t_value * objDir;
		Point3D capIntersect = objOrig + t_value * objDir;
	
		// Calculate x^2 + y^2 = r^2
		Vector3D capIntersectVect = Vector3D(capIntersect[0],capIntersect[1],capIntersect[2]);
		double rSquared = capIntersectVect.dot(capIntersectVect) - capIntersectVect[2]*capIntersectVect[2];
	
		// Check intersection with cap
		if (rSquared < 1) {
			// Found intersection, now check if intersection is with
			// z=1 or z=-1 plane and set the normal accordinly
			Vector3D n;
			if (capIntersect[2] < 0.5){ // Ray intersects z=0 plane, otherwise z=+1
				n = Vector3D(0, 0, -1);
			} else {
				n = Vector3D(0, 0, 1);
			}
	
			ray.intersection.none = false;
			ray.intersection.t_value = t_value;
			ray.intersection.point = modelToWorld * capIntersect;
			ray.intersection.normal = transNorm(worldToModel, n);
			ray.intersection.obj = this;
			return true;
		}
	}

	// No intersection with cap
	return false;
}

bool UnitCylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {

	// Perform intersection of unit cylinder centered at the origin:
	// 1. solve the quadratic equation (use the discriminent to see if there is an intersection)
	// 2. if discriminent is greater than zero:
	// 	i. check if the intersection is with the infinite cylinder
	//  i. if not, check if the intersection is with the top or bottom caps
	bool has_intersection = false;
	// Transform to object coordinates
	Point3D cylinderCenter(0, 0, 0);
	// Vector3D objOrig = worldToModel * ray.origin - cylinderCenter;
	Point3D objOrig = worldToModel * ray.origin;
	Vector3D objDir = worldToModel * ray.dir;

	// Ignore z coordinate to transform problem into 2D
	Vector3D origAsVector = Vector3D(objOrig[0],objOrig[1],objOrig[2]);
	Vector3D xyOrig = origAsVector;
	Vector3D xyDir = objDir;
	xyOrig[2] = 0.0;
	xyDir[2] = 0.0;

	// Check intersection to infinite cylinder by solving the quadratic equation
	// (as done in UnitSphere::intersect)
	double A = xyDir.dot(xyDir);
	double B = 2 * xyOrig.dot(xyDir);
	double C = xyOrig.dot(xyOrig) - 1;
	double D = (B*B) - (4*A*C);

	// Discriminent greater than zero means we have an intersection
	if (D < 0) {
		return false;
	}

	double t1 = ( -B + sqrt(D) ) / (2*A);
	double t2 = ( -B - sqrt(D) ) / (2*A);

	// Check if ray facing away from cylinder
	if (t1 < 0) {
		return false;
	}
	// Set t_value to the near face
	double t_value = t1;
	if (t2 >= 0) {
		t_value = t2;
	}
	
	// Check if this is the nearest intersection with infinite cylinder
	if ((ray.intersection.none || t_value < ray.intersection.t_value) &&
		t_value > M_EPSILON) {
		
		Point3D point = objOrig + t_value * objDir;
		Vector3D n = point - cylinderCenter;

		// Check if ray intersects cylinder side between z=-1 and z=+1
		if (point[2] > 0 && point[2] < 1){
			ray.intersection.none = false;
			ray.intersection.t_value = t_value;
			ray.intersection.point = modelToWorld * point;
			ray.intersection.normal = transNorm(worldToModel, n);
			ray.intersection.obj = this;
			// return true;
			has_intersection = true;
		}
	}

	// Check for intersection with top or bottom caps of cylinder (2D unit circles)
	// First see if ray direction is parallel to cylinder caps
	if (objDir[2] == 0 && !has_intersection){
		return false;
	}
	// Calculate t_values for z=1 and z=0 using the equation of the line and plane
	// (as done in UnitSquare::intersect)
	t1 = (1 - objOrig[2]) / objDir[2];
	// t2 = (-1 - objOrig[2]) /objDir[2];
	t2 = (objOrig[2]) /objDir[2];

	// Check if ray lies outside of the area between the two caps
	if(t1 > 0 && t2 > 0){
		double t_value = fmin(t1,t2);
		return (intersectCylinderCap(ray, worldToModel, modelToWorld, t_value) || has_intersection);
	// Check if ray lies inside the area between the two caps
	} else if ((t1 > 0 && t2 < 0) || (t1 < 0 && t2 > 0)) {
		double t_value = fmax(t1, t2);
		return (intersectCylinderCap(ray, worldToModel, modelToWorld, t_value) || has_intersection);
	} else {
		// Case: t1 < 0 && t2<0
		// This means ray lies outside the area between the caps and facing away from them
		return has_intersection;
	}
	return has_intersection;
}


//http://geomalgorithms.com/a06-_intersect-2.html
bool Triangle::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {

	// Transform ray from world to model coordinates
	Point3D rayOrig = worldToModel * ray.origin;
	Vector3D rayDir = worldToModel * ray.dir;


	Vector3D n = (_v1 - _v0).cross(_v2 - _v0);


	Vector3D u = (_v1 - _v0);

	Vector3D v = (_v2 - _v0);


	Vector3D w0 = (rayOrig -  _v0);

	double a = - (n.dot(w0));

	double b = n.dot(rayDir);

	if (fabs(b) < 0.00001) {     // ray is  parallel to triangle plane
		return false;
	}

	// get intersect point of ray with triangle plane
	double r = a / b;
	if (r < 0.0){                    // ray goes away from triangle
			return false;                   // => no intersect
	}
	// for a segment, also test if (r > 1.0) => no intersect

	Point3D point =  rayOrig + r * rayDir;

	float    uu, uv, vv, wu, wv, D;
    uu = u.dot(u);
    uv = u.dot(v);
    vv = v.dot(v);
		Vector3D w = point - _v0;
    wu = w.dot(u);
    wv = w.dot(v);
    D = uv * uv - uu * vv;


			// get and test parametric coords
	float s, t;
	s = (uv * wv - vv * wu) / D;
	if (s < 0.0 || s > 1.0){        // I is outside T
			return false;
	}
	t = (uv * wu - uu * wv) / D;
	if (t < 0.0 || (s + t) > 1.0){  // I is outside T
			return false;
	}



	if (ray.intersection.none || r <= ray.intersection.t_value) {
		ray.intersection.none = false;
		ray.intersection.point = modelToWorld * point;
		ray.intersection.t_value = r;
		ray.intersection.normal = transNorm(worldToModel, n);
		ray.intersection.normal.normalize();
		ray.intersection.obj = this;
		return true;
	}


	return false;
}


//
// bool Triangle::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
// 		const Matrix4x4& modelToWorld ) {
//
// 	// Transform ray from world to model coordinates
// 	Point3D rayOrig = worldToModel * ray.origin;
// 	Vector3D rayDir = worldToModel * ray.dir;
//
// 	Vector3D normal = (_v2 - _v0).cross(_v1 - _v0);
// 	normal.normalize();
//
// 	// Check if ray is parallel to triangle
// 	// face might no be non-planar.. average two triangle normals
// 	double determinant = rayDir.dot(normal);
//
// 	if (fabs(determinant) < M_EPSILON) {
// 		return false;
// 	}
//  	// Compute t_value of the equation of the line of intersection and
//  	// find the point of intersection
// 	double t_value = ((_v0 - rayOrig).dot(normal)) / determinant;
// 	Point3D point = rayOrig + t_value * rayDir;
//
//
// 	Vector3D sideOneA = (_v1 - _v0);
// 	Vector3D sideOneB = (_v0 - _v1);
//
// 	Vector3D sideTwoA = (_v2 - _v1);
// 	Vector3D sideTwoB = (_v1 - _v2);
//
// 	Vector3D sideThreeA = (_v0 - _v2);
// 	Vector3D sideThreeB = (_v2 - _v0);
//
// 	Vector3D normalOne = *(new Vector3D(sideOneA[1],sideOneB[0],0.0));
// 	double one = (point - _v0).dot(normalOne);
//
// 	Vector3D normalTwo = *(new Vector3D(sideTwoA[1],sideTwoB[0],0.0));
// 	double two = (point - _v1).dot(normalTwo);
//
// 	Vector3D normalThree = *(new Vector3D(sideThreeA[1],sideThreeB[0],0.0));
// 	double three = (point - _v2).dot(normalThree);
//
// 		// Triangle intersection test (3 half planes)
// 		if (one > 0 & two > 0 & three > 0) {
// 			// Found intersection, set ray attributes accordingly
// 			if ((ray.intersection.none || t_value <= ray.intersection.t_value)) {
// 				ray.intersection.none = false;
// 				ray.intersection.point = modelToWorld * point;
// 				ray.intersection.t_value = t_value;
// 				ray.intersection.normal = transNorm(worldToModel, normal);
// 				ray.intersection.normal.normalize();
// 				ray.intersection.obj = this;
// 				return true;
// 			}
// 		}
// 		return false;
// }

bool planeIntersection(Point3D origin, Vector3D dir, int plane, double min, double max,
	double * Tnear, double * Tfar) {

	if (dir[plane] == 0) {
		if (origin[plane] < min || origin[plane] > max) {
			return false;
		} else {
			*Tnear = fmax(*Tnear,min);
			*Tfar = fmin(*Tfar,max);
		}
	} else {
		double T1 = (min-origin[plane])/dir[plane];
		double T2 = (max-origin[plane])/dir[plane];
		if (T1>T2) {
			double tmp = T1;
			T1 = T2;
			T2 = tmp;
		}
		*Tnear = fmax(*Tnear,T1);
		*Tfar = fmin(*Tfar,T2);
	}

	if (*Tnear>*Tfar || *Tfar < 0) {
		return false;
	}

	return true;
}

bool MeshObject::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld ) {

	BoundingVolume* bv = _mesh->get_bounding_vol();
	Point3D _min = bv->min;
	Point3D _max = bv->max;

	double xmin = _min[0];
	double ymin = _min[1];
	double zmin = _min[2];
	double xmax = _max[0];
	double ymax = _max[1];
	double zmax = _max[2];
	double t_value;

	// Transform ray into object space
	Point3D origin = worldToModel * ray.origin;
	Vector3D dir = worldToModel * ray.dir;

	double Tnear = -INFINITY;
	double Tfar = INFINITY;

	if (planeIntersection(origin, dir, 0, xmin, xmax, &Tnear, &Tfar)) {
		if (planeIntersection(origin, dir, 1, ymin, ymax, &Tnear, &Tfar)) {
			if (planeIntersection(origin, dir, 2, zmin, zmax, &Tnear, &Tfar)) {
				return true;
			}
		}
	}

	return false;
}
