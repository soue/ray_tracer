/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		This file contains the interface and
		datastructures of the raytracer.
		Simple traversal and addition code to
		the datastructures are given to you.

***********************************************************/

#include "util.h"
#include "scene_object.h"
#include "light_source.h"

// Scene redering controls
const bool REFLECTION_ENABLED = false;
const bool REFRACTION_ENABLED = true;
const bool AA_ENABLED = false;

const int MOTION_BLUR_CYCLES = 1;	//number of cycles; 1 for no motion blue.
const int AA_SAMPLES = 1;			//anti-aliasing samples per pixel
const int MAX_RAY_DEPTH = 5;		//recursive depth for calling shadeRay
const int DOF_RAYS = 1;				//1 for no DOF, > 1 for DOF
const int GLOSSY_RAYS = 1;			//1 for perfect reflection, >1 for 'glossy' reflection

const double REFRACT_AIR = 1.0;
const double REFRACT_WATER = 1.33;
const double REFRACT_GLASS = 1.5;
const double REFRACT_GLASS2 = 1.01;
const double REFRACT_LOW = 0.667;

const double FOCAL_LEN = 4.0;

 //static Texture * ENV_MAPPING = new Texture("environments/villa.bmp");
 static Texture * ENV_MAPPING = new Texture("environments/skybox.bmp");
// static Texture * ENV_MAPPING = NULL;

// Linked list containing light sources in the scene.
struct LightListNode {
	LightListNode() : light(NULL), next(NULL) {}
	LightListNode( LightSource* light, LightListNode* next = NULL ) :
		light(light), next(next) {}
	~LightListNode() {
		if (!light) delete light;
	}
	LightSource* light;
	LightListNode* next;
};

// The scene graph, containing objects in the scene.
struct SceneDagNode {
	SceneDagNode() :
		obj(NULL), mat(NULL),
		next(NULL), parent(NULL), child(NULL), isBlur(false) {
	}

	SceneDagNode( SceneObject* obj, Material* mat ) :
		obj(obj), mat(mat), next(NULL), parent(NULL), child(NULL) {
		}

	~SceneDagNode() {
		if (!obj) delete obj;
		if (!mat) delete mat;
	}

	// Pointer to geometry primitive, used for intersection.
	SceneObject* obj;
	// Pointer to material of the object, used in shading.
	Material* mat;
	// Each node maintains a transformation matrix, which maps the
	// geometry from object space to world space and the inverse.
	Matrix4x4 trans;
	Matrix4x4 invtrans;

	// Internal structure of the tree, you shouldn't have to worry
	// about them.
	SceneDagNode* next;
	SceneDagNode* parent;
	SceneDagNode* child;


	//for motion blur
	bool isBlur;
	Vector3D blurTranslation;
};

class Raytracer {
public:
	Raytracer();
	~Raytracer();

	// Renders an image fileName with width and height and a camera
	// positioned at eye, with view vector view, up vector up, and
	// field of view fov.
	void render( int width, int height, Point3D eye, Vector3D view,
			Vector3D up, double fov, char* fileName );

	// Add an object into the scene, with material mat.  The function
	// returns a handle to the object node you just added, use the
	// handle to apply transformations to the object.
	SceneDagNode* addObject( SceneObject* obj, Material* mat ) {
		return addObject(_root, obj, mat);
	}

	// Add an object into the scene with a specific parent node,
	// don't worry about this unless you want to do hierarchical
	// modeling.  You could create nodes with NULL obj and mat,
	// in which case they just represent transformations.
	SceneDagNode* addObject( SceneDagNode* parent, SceneObject* obj,
			Material* mat );

	// Add a light source.
	LightListNode* addLightSource( LightSource* light );

	// Transformation functions are implemented by right-multiplying
	// the transformation matrix to the node's transformation matrix.

	// Apply rotation about axis 'x', 'y', 'z' angle degrees to node.
	void rotate( SceneDagNode* node, char axis, double angle );

	// Apply translation in the direction of trans to node.
	void translate( SceneDagNode* node, Vector3D trans );

	// Apply scaling about a fixed point origin.
	void scale( SceneDagNode* node, Point3D origin, double factor[3] );

private:
	// Allocates and initializes the pixel buffer for rendering, you
	// could add an interesting background to your scene by modifying
	// this function (ignoreAAFlag will ignore the antialiasing flag
	// if set to true).
	void initPixelBuffer(bool ignoreAAFlag);

	// Saves the pixel buffer to a file and deletes the buffer.
	void flushPixelBuffer(char *file_name);

	// Texture buffers (like pixel buffers)
	void initTextureBuffer();
	void flushTextureBuffer();

	// Return the colour of the ray after intersection and shading, call
	// this function recursively for reflection and refraction.
	// The depth parameter keeps track of the recursive calls, rendering stops at MAX_RAY_DEPTH.
	// The refract_index keeps track of the refraction index for the current ray casted
	Colour shadeRay( Ray3D& ray, int depth, double refract_index );

	//for motion blur
	void applyMotion(SceneDagNode*);


	// Constructs a view to world transformation matrix based on the
	// camera parameters.
	Matrix4x4 initInvViewMatrix( Point3D eye, Vector3D view, Vector3D up );

	// Traversal code for the scene graph, the ray is transformed into
	// the object space of each node where intersection is performed.
	void traverseScene( SceneDagNode* node, Ray3D& ray );

	// After intersection, calculate the colour of the ray by shading it
	// with all light sources in the scene.
	void computeShading( Ray3D& ray );

	// Given pixel(i,j) calculate the average pixel colour among the samples
	int *getAvePixelCol( int i, int j, unsigned char* rbuffer,
			unsigned char* gbuffer, unsigned char* bbuffer );

	// Apply anti-aliasing method (supersampling of size AA_SAMPLES*AA_SAMPLES)
	void applyAntiAliasing();

	// Apply jitter for anti-aliasing depth of field; jitter the origin of the rays casted
	Point3D applyDOFJitter( Point3D point, int size );
	// Apply jitter for anti-aliasing glossy reflection; jitter the direction of the multiple
	// reflected rays
	Vector3D applyGlossyJitter( Ray3D ray, Vector3D point );

	// Width and height of the viewport.
	int _scrWidth;
	int _scrHeight;

	// Light list and scene graph.
	LightListNode *_lightSource;
	SceneDagNode *_root;

	// Pixel buffer.
	unsigned char* _rbuffer;
	unsigned char* _gbuffer;
	unsigned char* _bbuffer;

	//textures
	unsigned char* _rtbuffer;
	unsigned char* _gtbuffer;
	unsigned char* _btbuffer;
	long int _theight;
	unsigned long int _twidth;

	// Refraction index
	double _refract_index;

	// Maintain global transformation matrices similar to OpenGL's matrix
	// stack.  These are used during scene traversal.
	Matrix4x4 _modelToWorld;
	Matrix4x4 _worldToModel;
};
