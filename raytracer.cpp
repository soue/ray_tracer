/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h,
		and the main function which specifies the
		scene to be rendered.

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>


int intersection_count;
int skips;

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

Vector3D computeNormal(Point3D v0, Point3D v1, Point3D v2) {
	Vector3D n = (v1-v0).cross(v2-v0);
	n.normalize();
	return n;
}

Vector3D computeNormal(Point3D v0, Point3D v1, Point3D v2, Point3D v3) {
	Vector3D n1 = (v1-v0).cross(v2-v1);
	n1.normalize();
	Vector3D n2 = (v2-v0).cross(v3-v2);
	n2.normalize();
	return (0.5*(n1+n2));
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent,
		SceneObject* obj, Material* mat ) {

	bool add_mesh = false;

	if (obj->get_type() == MESH) {
		MeshObject * meshObj = (MeshObject*)obj;
		if (meshObj->add_to_scene()) {
			meshObj->set_add_to_scene(false);
			vector< vector<int> > face_vertices = meshObj->get_mesh()->get_face_vertices();
 			vector<Point3D> vertices = meshObj->get_mesh()->get_vertices();
			SceneDagNode* mesh = addObject( meshObj, mat );

			for (int f = 0; f < face_vertices.size(); f++) {
				if (face_vertices[f].size() < 4) {
					Vector3D n = computeNormal(vertices[face_vertices[f][0]-1],
												vertices[face_vertices[f][1]-1],
												vertices[face_vertices[f][2]-1]);
					Triangle * triangle = new Triangle(vertices[face_vertices[f][0]-1],
													   vertices[face_vertices[f][1]-1],
													   vertices[face_vertices[f][2]-1],
													   n,
													   meshObj->get_texture());
					addObject(mesh, triangle, mat);
				} else {
					Vector3D n = computeNormal(vertices[face_vertices[f][0]-1],
												vertices[face_vertices[f][1]-1],
												vertices[face_vertices[f][2]-1],
												vertices[face_vertices[f][3]-1]);
					Triangle * triangle1 = new Triangle(vertices[face_vertices[f][0]-1],
													   vertices[face_vertices[f][1]-1],
													   vertices[face_vertices[f][2]-1],
													   n,
													   meshObj->get_texture());
					Triangle * triangle2 = new Triangle(vertices[face_vertices[f][0]-1],
													   vertices[face_vertices[f][2]-1],
													   vertices[face_vertices[f][3]-1],
													   n,
													   meshObj->get_texture());
					addObject(mesh, triangle1, mat);
					addObject(mesh, triangle2, mat);
				}
			}
			return mesh;
		} else {
			add_mesh = true;
		}
	}
	if (obj->get_type() != MESH || add_mesh) {
		SceneDagNode* node = new SceneDagNode( obj, mat );
		node->parent = parent;
		node->next = NULL;
		node->child = NULL;

		// Add the object to the parent's child list, this means
		// whatever transformation applied to the parent will also
		// be applied to the child.
		if (parent->child == NULL) {
			parent->child = node;
		}
		else {
			parent = parent->child;
			while (parent->next != NULL) {
				parent = parent->next;
			}
			parent->next = node;
		}
		return node;
	}
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;

	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation;
			angle = -angle;
		}
		else {
			node->invtrans = rotation*node->invtrans;
		}
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;

	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation;
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans;
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;

	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale;
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans;
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view,
		Vector3D up ) {
	Matrix4x4 mat;
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat;
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	SceneDagNode *childPtr;

	intersection_count++;

	bool skipMesh = false;
	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel;
	if (node->obj) {
		// Perform intersection.
		if (node->obj->get_type() == MESH) {
			if (!node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
				skipMesh = true;
			}
		} else if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		if (!skipMesh || childPtr->obj->get_type() != TRIANGLE) {

			skipMesh = false;
			traverseScene(childPtr, ray);
		} else {
			skips++;
		}
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;

		//If shadows are enabled we want to know if the intersection point has direct line of sight with a line source.
		if (SHADOWS_ENABLED) {

			int numShadowRays = SHADOW_RAYS;

			// Make sure shadows for this ray is set to zero
			ray.intersection.lightRayHits = 0;
			for (int sray = 0; sray < numShadowRays; sray++) {
				Point3D lightPosition = curLight->light->get_sample_position();
				Vector3D s =  ray.intersection.point - lightPosition;
				Ray3D shadowRay = Ray3D(lightPosition, s);
				traverseScene(_root, shadowRay);
				if (!shadowRay.intersection.none && (ray.intersection.point - shadowRay.intersection.point).length() < 0.0001) {
					ray.intersection.lightRayHits++;
				}else{
				}
			}
			if (ray.intersection.lightRayHits > 0) {
				curLight->light->shade(ray);
			}

		} else if (!SHADOWS_ENABLED) {
			curLight->light->shade(ray);
		}
		curLight = curLight->next;
	}
}


void Raytracer::initPixelBuffer(bool ignoreAAFlag) {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	if (AA_ENABLED && !ignoreAAFlag){
		numbytes = numbytes * AA_SAMPLES * AA_SAMPLES;
	}
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

void Raytracer::initTextureBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	if (AA_ENABLED){
		numbytes = numbytes * AA_SAMPLES * AA_SAMPLES;
	}
	_rtbuffer = new unsigned char[numbytes];
	_gtbuffer = new unsigned char[numbytes];
	_btbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rtbuffer[i*_scrWidth+j] = 0;
			_gtbuffer[i*_scrWidth+j] = 0;
			_btbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushTextureBuffer() {
	delete _rtbuffer;
	delete _gtbuffer;
	delete _btbuffer;
}

int* Raytracer::getAvePixelCol(int row, int col, unsigned char* rbuffer,
							   unsigned char* gbuffer, unsigned char* bbuffer){
	int row1 = row * AA_SAMPLES;
	int row2 = row1 + AA_SAMPLES - 1;
	int col1 = col * AA_SAMPLES;
	int col2 = col1 + AA_SAMPLES - 1;
	int pix = (row2 - row1) * (col2 - col1);
	int ind;
	int c[] = {0, 0, 0};

	for (int i = row1; i < row2; i++){
		for (int j = col1; j < col2; j++){
			ind = i * _scrWidth * AA_SAMPLES + j;
			c[0] += rbuffer[ind];
			c[1] += gbuffer[ind];
			c[2] += bbuffer[ind];
		}
	}

	for (int k = 0; k < 3; k++) {
		c[k] = c[k]/pix;
	}
	return c;
}

void Raytracer::applyAntiAliasing() {
	unsigned char* tmp_rbuffer = _rbuffer;
	unsigned char* tmp_gbuffer = _gbuffer;
	unsigned char* tmp_bbuffer = _bbuffer;
	// Reset pixel buffers to original image size
	initPixelBuffer(true);
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			int * c = getAvePixelCol(i,j, tmp_rbuffer, tmp_gbuffer, tmp_bbuffer);
			_rbuffer[i*_scrWidth+j] = c[0];
			_gbuffer[i*_scrWidth+j] = c[1];
			_bbuffer[i*_scrWidth+j] = c[2];
		}
	}
}

Point3D Raytracer::applyDOFJitter(Point3D point, int size) {
	Point3D jitteredPoint = point;
	double jitter;
	for (int j = 0; j < size; j++) {
		jitter = (rand() % 101) / 500.0;
		jitteredPoint[j] = jitteredPoint[j] + jitter;
	}
	return jitteredPoint;
}

Vector3D Raytracer::applyGlossyJitter(Ray3D ray, Vector3D dir) {
	// double jitter1 = (rand() % 101) / 100.0;
	// double jitter2 = (rand() % 101) / 100.0;
	// double theta1 = 2 * M_PI * jitter2;
	// double theta2 = acos(pow((1-jitter1), ray.intersection.mat->reflectance));
	// double x = sin(theta1) * cos(theta2);
	// double y = sin(theta1) * sin(theta2);
	// Vector3D u = dir.cross(ray.intersection.normal);
	// Vector3D v = dir.cross(u);
	// Vector3D reflectedDir = (sin(theta2) * cos(theta1) * u) + (sin(theta2) * sin(theta1) *v) + dir;
	// reflectedDir.normalize();

	Vector3D reflectedDir = dir;
	double jitter;
	for (int j = 0; j < 3; j++) {
		jitter = (rand() % 101) / 1000.0;
		reflectedDir[j] = reflectedDir[j] + jitter;
	}
	reflectedDir.normalize();
	return reflectedDir;
}

// Modified to take in a depth and refraction index parameters
// depth keeps track of the number of light bounces, maxes out at MAX_RAY_DEPTH and stops
// refract_index keeps track of the refraction index of the current ray
Colour Raytracer::shadeRay( Ray3D& ray, int depth, double refract_index ) {
	Colour col(0.0, 0.0, 0.0);


	if (depth < MAX_RAY_DEPTH) {
		traverseScene(_root, ray);

		// Don't bother shading if the ray didn't hit anything.
		if (!ray.intersection.none) {
			computeShading(ray);
			col = ray.col;

			// Recursive ray tracing for reflection
			if (REFLECTION_ENABLED && ray.intersection.mat->reflectance > 0) {
				ray.dir.normalize();
				Vector3D b = ray.origin - ray.intersection.point;
				b.normalize();
				ray.intersection.normal.normalize();
				Vector3D reflectedDir = -b + 2*ray.intersection.normal.dot(b)*ray.intersection.normal;
				// Vector3D reflectedDir =  -1 * ray.dir + 2 * ray.intersection.normal.dot(ray.dir) * ray.intersection.normal;
				Ray3D reflectedRay;
				reflectedRay.dir = reflectedDir;
				reflectedRay.origin = ray.intersection.point;

				Colour reflectedCol(0.0, 0.0, 0.0);
				for (int g = 0; g < GLOSSY_RAYS; g++) {
					if (GLOSSY_RAYS > 1) {
						// Apply jitter to reflected direction for glossy effect
						reflectedRay.dir = applyGlossyJitter(ray, reflectedDir);
					}
					reflectedCol = reflectedCol +
						(ray.intersection.mat->specular) * shadeRay(reflectedRay, depth+1, refract_index);
				}
				col = (1-ray.intersection.mat->reflectance)*col + ((1.0/GLOSSY_RAYS) * (1/(depth+1)) * reflectedCol);
			}

			// Refraction ray tracing for translucent objects
			if (REFRACTION_ENABLED && ray.intersection.mat->opacity > 0) {

				// Snell's law: sin(theta_1)/sin(theta_2) = c1/c2
 				// Assume light travels between air and an object (or the reverse)

				Vector3D rayDir = ray.dir;
				rayDir.normalize();
				Vector3D n = ray.intersection.normal;
				n.normalize();

				// incident ray cos theta
				double cosI = n.dot(rayDir);
				double r_index, sinthetat, nDotNegD;
				Vector3D refractedDir;

				bool internal = false;
				if (cosI < 0.0) { // air to object
					r_index = 1.0 / ray.intersection.mat->light_speed;
					nDotNegD = n.dot(-rayDir);
					double thetai = acos(nDotNegD);
					sinthetat = r_index * sin(thetai);
					if(fabs(sinthetat) > 1.0/r_index){
						// printf("YES INTERAL REFRACTION ALSO ADD IT FOR entering!\n");
						internal = true;

					}
					sinthetat = sinthetat*sinthetat;
				} else { // object to air
					r_index = ray.intersection.mat->light_speed / 1.0;
					n = -n;
					n.normalize();
					nDotNegD = n.dot(-rayDir);
					sinthetat = r_index*r_index*(1-(nDotNegD*nDotNegD));
					if(sinthetat > 1.0/r_index){
						// printf("YES INTERAL REFRACTION ALSO ADD IT FOR LEAVING!\n");
						internal = true;
					}
				}

				// if (internal) {
// 
// 					ray.dir.normalize();
// 					Vector3D b = ray.origin - ray.intersection.point;
// 					b.normalize();
// 					ray.intersection.normal.normalize();
// 					Vector3D reflectedDir = -b + 2*ray.intersection.normal.dot(b)*ray.intersection.normal;
// 					// Vector3D reflectedDir =  -1 * ray.dir + 2 * ray.intersection.normal.dot(ray.dir) * ray.intersection.normal;
// 					Ray3D reflectedRay;
// 					reflectedRay.dir = reflectedDir;
// 					reflectedRay.origin = ray.intersection.point;
// 					col = (1-ray.intersection.mat->reflectance)*col + 
// 							(ray.intersection.mat->opacity) * shadeRay(reflectedRay, depth+1, refract_index);
// 				} else {
					double cosT2 = 1.0 - pow(r_index, 2) * (1 - pow(nDotNegD, 2));
					if (cosT2 > 0.0) {
						// Shade refracted ray
						refractedDir = ((r_index * nDotNegD - sqrt(cosT2)) * n) - (r_index * (-rayDir));
						refractedDir.normalize();

						Ray3D refractedRay = Ray3D(ray.intersection.point + M_EPSILON * refractedDir, refractedDir);
						Colour refractiveCol = shadeRay(refractedRay, depth + 1, r_index);
						col = (1-ray.intersection.mat->opacity)*col+ray.intersection.mat->opacity * refractiveCol;
					}

				// }
			}
		} else if (ENV_MAPPING != NULL) {
				ray.dir.normalize();
				col = ENV_MAPPING->getEnvironmentColour(ray.dir);
		}
	}
	col.clamp();
	return col;
}

void Raytracer::applyMotion(SceneDagNode* node){
	if(node->isBlur){
		translate(node, 0.25 * (node->blurTranslation));
	}
	SceneDagNode* childPtr = node->child;
	while(childPtr != NULL){
		applyMotion(childPtr);
		childPtr = childPtr->next;
	}
}

void Raytracer::render( int width, int height, Point3D eye, Vector3D view,
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_refract_index = REFRACT_AIR;
	_scrWidth = width;
	_scrHeight = height;
	int aaWidth = width;
	int aaHeight = height;
	if (AA_ENABLED){
		aaWidth = aaWidth * AA_SAMPLES;
		aaHeight = aaHeight * AA_SAMPLES;
	}
	double factor = (double(aaHeight)/2)/tan(fov*M_PI/360.0);
	initPixelBuffer(false);
	viewToWorld = initInvViewMatrix(eye, view, up);

	//sample over time for motion blur
	//to keep runtime low we sample 4 times. Could increase this to make it look better
	int total = MOTION_BLUR_CYCLES;
	int number = ((MOTION_BLUR_CYCLES)*(MOTION_BLUR_CYCLES+1))/2;

	//This is a weight average.
	int sum_number = (1 - pow(2, total+1))/(1-(2))-1;

	for(int blur = 0; blur < MOTION_BLUR_CYCLES; blur++){

		// Construct a ray for each pixel.
		for (int i = 0; i < aaHeight; i++) {
			for (int j = 0; j < aaWidth; j++) {
				// Sets up ray origin and direction in view space,
				// image plane is at z = -1.
				Point3D origin(0, 0, 0);
				Point3D imagePlane;
				imagePlane[0] = (-double(aaWidth)/2 + 0.5 + j)/factor;
				imagePlane[1] = (-double(aaHeight)/2 + 0.5 + i)/factor;
				imagePlane[2] = -1;

				// TODO: Convert ray to world space and call
				// shadeRay(ray) to generate pixel colour.
				Ray3D ray;
				Colour col(0.0, 0.0, 0.0);

				// Ray properties if DOF is not set
				int dofRays = DOF_RAYS;
				Vector3D rayDir = imagePlane - origin;
				Point3D point = viewToWorld * imagePlane;

				// Ray properties if DOF is set
				if (DOF_RAYS > 1) {
					// DOF requires multiple rays casted and calculating the points in focus
					double Dtol = (Point3D(0,0,-1) - origin).length();
					double DtoP = rayDir.length();
					rayDir.normalize();
					point = eye + (DtoP/(Dtol/(Dtol+FOCAL_LEN)))*rayDir;
        }



				// Generate pixel colour by summing up the colour returned
				// by casting dofRays (1 if DOF not enabled)
				Point3D jitteredOrigin = origin;
				for (int r = 0; r < dofRays; r++){
					if (DOF_RAYS > 1){
						// Generate random jitters for eye position between 0-1
						jitteredOrigin = applyDOFJitter(origin, 2);
					}
					ray.origin =  viewToWorld * jitteredOrigin;
					ray.dir = point - ray.origin;
					ray.dir.normalize();
					col = col + shadeRay(ray, 0, _refract_index);
				}

				// Divide the summed colour by the total number of rays casted
				double d = dofRays;
				col = (1.0/d) * col;
				col.clamp();
				//for motion blur
				//We use a weight average.
				if(MOTION_BLUR_CYCLES > 0){
					int cycle_number = blur+1;
					double num = pow(2, cycle_number);
					double dom = sum_number;
					double percentage = (num/sum_number);
					col = (percentage) * col;
					// printf("percentage: %f %f %f \n",percentage,num,dom);
				}

				col.clamp();

				_rbuffer[i*aaWidth+j] += int(col[0]*255);
				_rbuffer[i*aaWidth+j] = _rbuffer[i*aaWidth+j] <= 255 ? _rbuffer[i*aaWidth+j] : 255;
				_gbuffer[i*aaWidth+j] += int(col[1]*255);
				_gbuffer[i*aaWidth+j] = _gbuffer[i*aaWidth+j] <= 255 ? _gbuffer[i*aaWidth+j] : 255;
				_bbuffer[i*aaWidth+j] += int(col[2]*255);
				_bbuffer[i*aaWidth+j] = _bbuffer[i*aaWidth+j] <= 255 ? _bbuffer[i*aaWidth+j] : 255;

			}
		}

		std::cout << "A" << std::endl;
		applyMotion(_root);
		std::cout << "B" << std::endl;
	}

	// Antialiasing super sampling - get the average value form the 3x3 grid of samples and
	// assign it to the appropriate pixel
	if (AA_ENABLED) {
		applyAntiAliasing();
	}
	flushPixelBuffer(fileName);
}

int main(int argc, char* argv[]){

	intersection_count = 0;
	skips = 0;

// Build your scene and setup your camera here, by calling
// functions from Raytracer.  The code here sets up an example
// scene and renders it from two different view points, DO NOT
// change this if you're just implementing part one of the
// assignment.
	Raytracer raytracer;
	int width = 320;
	int height = 240;
	// //
	// int width = 120;
	// int height = 80;

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}




// Materials for shading

	Texture brick("textures/brick.bmp");
	Material brickRed(Colour(0.2295, 0.08825, 0.0275),
	 				 Colour(0.5508, 0.2118, 0.066),
	 				 Colour(0.580594, 0.223257, 0.0695701),
						1000.2,		//specularity
						0.0, 		//reflection coefficient
						0.0, 		//opacity
						0.0);		//refraction
	brickRed.isTexture = true;
	brickRed.texture = &brick;

	Material glass( Colour(0.05, 0.05, 0.5),
				Colour(0.05, 0.05, 0.5),
				Colour(0.5, 0.5, 0.5),
				55,
				0.1,
				0.75,
				1.00);
	Material glass2( Colour(0.05, 0.05, 0.5),
				Colour(0.05, 0.05, 0.5),
				Colour(0.5, 0.5, 0.5),
				55,
				0.1,
				0.90,
				REFRACT_WATER);
	Material glass3( Colour(0.05, 0.05, 0.5),
				Colour(0.05, 0.05, 0.5),
				Colour(0.5, 0.5, 0.5),
				55,
				0.2,
				0.90,
				1.01);

Material glass4( Colour(0.05, 0.05, 0.5),
			Colour(0.05, 0.05, 0.5),
			Colour(0.5, 0.5, 0.5),
			55,
			0.0,
			0.70,
			1.01);

	Material mirror( Colour(0.05, 0.05, 0.5),
				Colour(0.05, 0.05, 0.5),
				Colour(0.5, 0.5, 0.5),
				102,
				1.0,
				0.0,
				1.00);

	Material mirror2(Colour(0.0, 0.0, 0.0),
	Colour(0.0, 0.0, 0.0),
	Colour(0.7, 0.7, 0.7),
	102,
	1.0,
	0.0,
	1.00);

	Material gold( Colour(0.3, 0.3, 0.3),
	 			   Colour(0.75164, 0.60648, 0.22648),
	 			   Colour(0.628281, 0.555802, 0.366065),
	 			   51.2, 			//specular exponent
	 			   0.0, 			//reflection coefficient
	 			   0.0, 			//opacity
	 			   1.0); 			//refraction index

	Material jade( Colour(0.2, 0.2, 0.2),
				   Colour(0.54, 0.89, 0.63),
				   Colour(0.316228, 0.316228, 0.316228),
				   100.0,
	 			   0.0,
	 			   0.0,
	 			   REFRACT_AIR);

	Material iceblue( Colour(0, 0, 0),
		Colour(100.0 / 255, 149.0 / 255, 237.0 / 255),
		Colour(0.316228, 0.316228, 0.316228),
				   100.0,
	 			   0.5,
	 			   0.0,
	 			   REFRACT_AIR);


	Material copper( Colour(0.2295, 0.08825, 0.0275),
	 				 Colour(0.5508, 0.2118, 0.066),
	 				 Colour(0.580594, 0.223257, 0.0695701),
	 				 20.2,
	 				 0.0,
	 				 0.0,
	 				 REFRACT_WATER);
// Defines a point light source.


	raytracer.addLightSource(new PointLight(Point3D(0,2, 2),Colour(0.9, 0.9, 0.9)));

// Apply some transformations to the unit square.
	double factor1[3] = { 2.0, 2.0, 2.0 };
	double factor2[3] = { 26.0, 26.0, 26.0 };
	double factor3[3] = { 1.0, 1.0, 2.0 };
	double factor4[3] = { 1.2, 1.2, 1.2 };
	double factor5[3] = { 0.8, 0.8, 0.8 };
	double factor6[3] = { 1000.0, 1000.0, 1000.0 };
	double factor7[3] = { 1.0, 1.0, 3.0 };
	//
	// SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(NULL), &mirror2 );
// 
// 	raytracer.translate(sphere2,Vector3D(3, 0, -4));
// 
// 
// 	Mesh* m0 = new Mesh("obj/bunny.obj");
// 
// 	SceneDagNode* sphere5 = raytracer.addObject(new MeshObject(m0, NULL), &glass4 );
// 	// sphere->isBlur = true;
// 	// sphere->blurTranslation = Vector3D(0.1,0,0);
// 	raytracer.translate(sphere5,Vector3D(-3, -1, -4));
// 	raytracer.rotate(sphere5, 'y', -180);
// 	raytracer.scale(sphere5, Point3D(0, 0, 0), factor1);
	//
	//
	// Mesh* m = new Mesh("obj/bunny.obj");
// 	SceneDagNode* bunny1 = raytracer.addObject( new MeshObject(m, NULL), &brickRed );
// 	raytracer.translate(bunny1, Vector3D(0.7, -1, -4.2));
// 
// 	Mesh* m2 = new Mesh("obj/cube_quad.obj");
// 	SceneDagNode* bunny2 = raytracer.addObject( new MeshObject(m2, NULL), &jade );
// 	raytracer.translate(bunny2, Vector3D(-0.7, -1, -4.2));
// 	raytracer.rotate(bunny2, 'y', -180);

	SceneDagNode* c = raytracer.addObject( new UnitCylinder(NULL), &glass );
	raytracer.translate(c, Vector3D(0, 0, -6));
	// raytracer.rotate(c, 'x', -90);
	raytracer.scale(c, Point3D(0, 0, 0), factor7);

	// SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(NULL), &copper );
	// raytracer.translate(sphere2, Vector3D(0, 0, -1.8));

	// SceneDagNode* plane = raytracer.addObject( new UnitSquare(NULL), &iceblue );
// 	raytracer.translate(plane, Vector3D(0, -1, 0));
// 	raytracer.rotate(plane, 'x', -90);
// 	raytracer.scale(plane, Point3D(0, 0, 0), factor6);


	// SceneDagNode* globe = raytracer.addObject( new UnitSphere(NULL), &glass3 );
	// raytracer.translate(globe, Vector3D(0, 0, -2));
	// raytracer.scale(globe, Point3D(0, 0, 0), factor1);





	std::clock_t    start;

	Vector3D up(0, 1, 0);
	double fov = 60;


//
// // Camera parameters.
// 	Point3D eye(0, 0, 1);
// 	Vector3D view(0, 0, -1);
//
  start = std::clock();
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");

	// Render it from a different point of view.
	Point3D eye3(4, 2, 1);
	Vector3D view3(-4, -2, -6);
	raytracer.render(width, height, eye3, view3, up, fov, "view2.bmp");


	Point3D eye4(-4, 2, 1);
	Vector3D view4(4, -2, -6);
	raytracer.render(width, height, eye4, view4, up, fov, "view3.bmp");

	std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	std::cout << "intersection count=" << intersection_count << "\n";
	std::cout << "skip count=" << skips << "\n";

	return 0;
}
