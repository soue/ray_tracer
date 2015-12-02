/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include <cstdlib>
#include "light_source.h"
#include "scene_object.h"



void applyTexture( Ray3D& ray ) {
	if (ray.intersection.obj->get_texture() != NULL) {
		Texture* texture = ray.intersection.obj->get_texture();
		int objType = ray.intersection.obj->get_type();
		if (objType == SQUARE) {
			Point3D p = ray.intersection.point;
			double u = p[0]/p[2];
			double v = p[1]/p[2];
			std::cout << "p=" << p << " ; u,v=" << u << ", " << v << "\n";
			ray.col = ray.col + texture->getUVMapping(u,v);
		} else if (objType == SPHERE) {
			//sphere origin=(x_c,y_c,z_c)
			//x=x_c+rcos(phi)sin(theta)
			//y=y_c_rsin(phi)cos(theta)
			//z=z_c+rcos(theta)
			//


			//theta=arccos((z-z_c)/r)
			//phi=arctan((y-y_c)/(x-x_c))
			//
			//u=phi/2pi
			//v=(pi-theta)/pi
			Vector3D N = ray.intersection.normal;
			double u = 0.5 + atan2(N[1],N[0])/(2.0*M_PI);
			double v = 0.5 - asin(N[2])/M_PI;
			ray.col = ray.col + texture->getUVMapping(u,v);
		} else if (objType == CYLINDER) {
			//h=z/height
			//theta=arctan(y/x)
			//
			//u=theta/2pi
			//v=h
			double radius = 1.0;
			Point3D p = ray.intersection.point;
			double height = p[2]/4.0*radius;
			double theta = atan2(p[1],p[0]);
			double v = theta/(2.0*M_PI);
			double u = height;
			ray.col = ray.col + texture->getUVMapping(u,v);
		} else if (objType == TRIANGLE) {
			//compute barycentric coordinates
			Triangle* triangle = (Triangle*)ray.intersection.obj;
			Point3D p1 = triangle->get_vertex1();
			Point3D p2 = triangle->get_vertex2();
			Point3D p3 = triangle->get_vertex3();
			Point3D p = ray.intersection.point;

			//vectors from each vertex to intersection point
			Vector3D A1 = p1 - p;
			Vector3D A2 = p2 - p;
			Vector3D A3 = p3 - p;

			//areas of triangles created by vectors to intersection point
			// double a = (((p1-p2).cross(p1-p3)).length())/2.0;
			double a1 = ((A2.cross(A3)).length())/(2.0);
			double a2 = ((A3.cross(A1)).length())/(2.0);
			double a3 = ((A1.cross(A2)).length())/(2.0);

			double u = p[0]*a1/(a1+a2+a3);
			double v = p[1]*a2/(a1+a2+a3);
			ray.col = ray.col + texture->getUVMapping(u,v);
		}
	}
}



void PointLight::shade( Ray3D& ray ) {
	ray.col = Colour(0.0, 0.0, 0.0);

	// TODO: implement this function to fill in values for ray.col
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray
	// is available.  So be sure that traverseScene() is called on the ray
	// before this function.

	// Specular alpha constant
	double alpha = ray.intersection.mat->specular_exp;
	// Normal
	Vector3D n = ray.intersection.normal;
	n.normalize();
	// Light direction
	Vector3D s = this->get_position() - ray.intersection.point;
	s.normalize();
	// Ray direction
	Vector3D d = ray.dir;
	d.normalize();
	// Reflection direction
	Vector3D m = s - 2 * n.dot(s) * n;
	m.normalize();

	// Calculate ambient, diffuse and specular terms for phong shading
	Colour ambientTerm = this->_col_ambient * ray.intersection.mat->ambient;
	Colour diffuseTerm = fmax(0, n.dot(s)) * this->_col_diffuse * ray.intersection.mat->diffuse;
	Colour specularTerm = pow(fmax(0, m.dot(d)), alpha) * this->_col_specular * ray.intersection.mat->specular;

	applyTexture(ray);

	if (SHADOWS_ENABLED) {
		if (ray.intersection.lightRayHits == 0) {
			// Only ambient term for point in full shadow
			ray.col = ray.col + ambientTerm;
		} else {
			if (SHADE_MODE == AMBIENT_DIFFUSE) {
				ray.col = ray.col + ambientTerm;
			} else if (SHADE_MODE == PHONG) {
				ray.col = ray.col + ambientTerm + specularTerm;
			}
			ray.col = ray.col + diffuseTerm;
		}
	} else {
		if (SHADE_MODE == AMBIENT_DIFFUSE) {
			ray.col = ray.col + ambientTerm;
		} else if (SHADE_MODE == PHONG) {
			ray.col = ray.col + ambientTerm + specularTerm;
		}
		ray.col = ray.col + diffuseTerm;
	}

	if (ray.intersection.mat->isTexture){
		ray.col = ray.col + ray.intersection.mat->texture->getColour(ray.intersection.point[0], ray.intersection.point[1]);
	}

 	ray.col.clamp();
}

void AreaLight::shade( Ray3D& ray ) {
	ray.col = Colour(0.0, 0.0, 0.0);


	// Specular alpha constant
	double alpha = ray.intersection.mat->specular_exp;
	// Normal
	Vector3D n = ray.intersection.normal;
	n.normalize();
	// Light direction
	Vector3D s = this->get_position() - ray.intersection.point;
	s.normalize();
	// Ray direction
	Vector3D d = ray.dir;
	d.normalize();
	// Reflection direction
	Vector3D m = s - 2 * n.dot(s) * n;
	m.normalize();

	// Calculate ambient, diffuse and specular terms for phong shading
	Colour ambientTerm = this->_col_ambient * ray.intersection.mat->ambient;
	Colour diffuseTerm = fmax(0, n.dot(s)) * this->_col_diffuse * ray.intersection.mat->diffuse;
	Colour specularTerm = pow(fmax(0, m.dot(d)), alpha) * this->_col_specular * ray.intersection.mat->specular;

	//Check this!
	applyTexture(ray);
	if (SHADOWS_ENABLED) {
		// Include ambient term and only some diffuse and specular depending on how much the
		// point is in shadow
		// ray.col = ambientTerm +
		// 		  ((ray.intersection.lightRayHits/SHADOW_RAYS) * (diffuseTerm + specularTerm));
		double one = ray.intersection.lightRayHits;
		double two = SHADOW_RAYS;
		double example = one/two;
		ray.col = ambientTerm+(example)*(diffuseTerm + specularTerm);
	} else {
		ray.col = ambientTerm + diffuseTerm + specularTerm;
	}

	if (ray.intersection.mat->isTexture){
		ray.col = ray.col + ray.intersection.mat->texture->getColour(ray.intersection.point[0], ray.intersection.point[1]) ;
	}

	ray.col.clamp();
}


// Point3D AreaLight::get_sample_position() {
// 	// Generate random sample point on from area light plane to cast multiple shadow rays.
// 	// This is for implementing soft shadowing.
//
// 	double x = ((rand() % (201))-100.0)/100;
//
// 	double y = ((rand() % (201))-100.0)/100;
//
// //	printf("%d\n",(rand() % (30001)));
//
// 	double z = _top_left[2];
// 	Point3D sample_point2 = Point3D(x, y, z);
// 	printf("(%f,%f,%f)\n",sample_point2[0],sample_point2[1],sample_point2[2]);
// 	return sample_point2;
// }


Point3D AreaLight::get_sample_position() {
	// Generate random sample point on from area light plane to cast multiple shadow rays.
	// This is for implementing soft shadowing.
	double x1 = _top_right[0];
	double y1 = _top_left[1];
	int x2 = _top_right[0] - _top_left[0];
	int y2 = _top_right[1] - _bottom_right[1];
	// printf("(%d,%d)\n",x2,y2);
	double sample1 =(rand() % ((x2)*100)+1);
	double sample2 =(rand() % ((y2)*100)+1);
	// printf("(%f,%f)\n",sample1,sample2);



	double x = (rand() % ((x2)*100+1) + _top_left[0]*100)/100;
	double y = (rand() % ((y2)*100+1) + _bottom_right[1]*100)/100;
	double z;

	if (abs(_normal[1]) < M_EPSILON)
		z = _top_right[2];
	else
		z = (-_normal[0] * (x - _top_left[0])  - _normal[1] * (y - _top_left[1])) /
				_normal[1] + _top_left[2];

	Point3D sample_point2 = Point3D(x, y, z);
	// printf("(%f,%f,%f)\n",sample_point2[0],sample_point2[1],sample_point2[2]);

	return sample_point2;
}

// Point3D AreaLight::get_sample_position() {
// 	// Generate random sample point on from area light plane to cast multiple shadow rays.
// 	// This is for implementing soft shadowing.
// 	double x1 = _top_right[0];
// 	double y1 = _top_left[1];
// 	int x2 = _top_right[0] - _top_left[0];
// 	int y2 = _bottom_right[1] - _top_right[1];
// 	printf("(%d,%d)\n",x2,y2);
//
// 	double x = rand() % x2 + x1;
// 	double y = rand() % y2 + y1;
// 	double z;
//
// 	if (abs(_normal[1]) < M_EPSILON)
// 		z = _top_right[2];
// 	else
// 		z = (-_normal[0] * (x - _top_left[0])  - _normal[1] * (y - _top_left[1])) /
// 				_normal[1] + _top_left[2];
//
// 	Point3D sample_point2 = Point3D(x, y, z);
// 	printf("(%f,%f,%f)\n",sample_point2[0],sample_point2[1],sample_point2[2]);
//
// 	return sample_point2;
// }


Vector3D AreaLight::get_light_normal() {
	Vector3D l1 = this->_top_left - this->_top_right;
	Vector3D l2 = this->_top_left - this->_bottom_left;
	Vector3D n = l1.cross(l2);
	n.normalize();
	return n;
}

Point3D AreaLight::get_light_center() {
	// Get the center point of the area light
	int p1 = (_top_right[0] + _top_left[0])/2;
	int p2 = (_bottom_right[1] + _top_right[1])/2;
	return Point3D(p1,p2,_top_left[2]);
}

void SphereLight::shade( Ray3D& ray ) {

	ray.col = Colour(0.0, 0.0, 0.0);
	// Specular alpha constant
	double alpha = ray.intersection.mat->specular_exp;
	// Normal
	Vector3D n = ray.intersection.normal;
	n.normalize();
	// Light direction
	Vector3D s = this->get_position() - ray.intersection.point;
	s.normalize();
	// Ray direction
	Vector3D d = ray.dir;
	d.normalize();
	// Reflection direction
	Vector3D m = s - 2 * n.dot(s) * n;
	m.normalize();

	// Calculate ambient, diffuse and specular terms for phong shading
	Colour ambientTerm = this->_col_ambient * ray.intersection.mat->ambient;
	Colour diffuseTerm = fmax(0, n.dot(s)) * this->_col_diffuse * ray.intersection.mat->diffuse;
	Colour specularTerm = pow(fmax(0, m.dot(d)), alpha) * this->_col_specular * ray.intersection.mat->specular;

	//Check this!

	// // applyTexture(ray);
	// if (SHADOWS_ENABLED) {
	// 	// Include ambient term and only some diffuse and specular depending on how much the
	// 	// point is in shadow
	// 	// ray.col = ambientTerm +
	// 	// 		  ((ray.intersection.lightRayHits/SHADOW_RAYS) * (diffuseTerm + specularTerm));
	// 	double one = ray.intersection.lightRayHits;
	// 	double two = SHADOW_RAYS;
	// 	double example = one/two;
	// 	ray.col = (example)*(ambientTerm+diffuseTerm+specularTerm);
	// } else {
	// 	ray.col = ambientTerm + diffuseTerm + specularTerm;
	// }


	if (SHADOWS_ENABLED) {
		if (ray.intersection.lightRayHits == 0) {
			// Only ambient term for point in full shadow
			ray.col = ray.col + ambientTerm;
		} else {
				double one = ray.intersection.lightRayHits;
				double two = SHADOW_RAYS;
				double example = one/two;


			if (SHADE_MODE == AMBIENT_DIFFUSE) {
				ray.col = ray.col + ambientTerm;
			} else if (SHADE_MODE == PHONG) {
				ray.col = ray.col + ambientTerm + specularTerm;
			}
			ray.col = ray.col + diffuseTerm;
			ray.col = (example)*(ray.col);


		}
	} else {
		if (SHADE_MODE == AMBIENT_DIFFUSE) {
			ray.col = ray.col + ambientTerm;
		} else if (SHADE_MODE == PHONG) {
			ray.col = ray.col + ambientTerm + specularTerm;
		}
		ray.col = ray.col + diffuseTerm;
	}

	if (ray.intersection.mat->isTexture){
		ray.col = ray.col + ray.intersection.mat->texture->getColour(ray.intersection.point[0], ray.intersection.point[1]);
	}

	if (ray.intersection.mat->isTexture){
	//	printf("Need to add textures to arealight \n");
	}

	// ray.col = (1 - ray.intersection.mat->opacity) * ray.col; // transparency effect on shadows
	ray.col.clamp();
}




Point3D SphereLight::get_sample_position() {
	double x;
	double y;
	double z;

	double u = rand();
	double v = rand();

	x = _radius*sin(u)*cos(v)+_center[0];
	y = _radius*cos(u)*cos(v)+_center[1];
	z = _radius*sin(v)+_center[2];


	Point3D sample_point = Point3D(x, y, z);
	// printf("(%f,%f,%f)\n",sample_point[0],sample_point[1],sample_point[2]);
	// printf("(%f)\n",total);

	return sample_point;
}


Point3D SphereLight::get_light_center() {

	return _center;
}
