#include "util.h"
// #include "bounding_volume.h"

#include <vector>
#include <string>
#include <stdio.h>
#include <sstream>
#include <cstring>

using std::vector;
using std::string;
using std::stringstream;

const int BUFFER = 256;

struct BoundingVolume {
	BoundingVolume(Point3D pmin, Point3D pmax) : min(pmin), max(pmax) {}
	Point3D min;
	Point3D max;
};

class Mesh {
public:
	Mesh(const char * file_name);
	~Mesh() {};
	vector< vector<int> > get_face_vertices() { return _face_vertices; }
	vector< vector<int> > get_face_textures() { return _face_textures; }
	vector< vector<int> > get_face_normals() { return _face_normals; }
	vector<Point3D> get_vertices() { return _vertices; }
	vector<Point3D> get_texture_coords() { return _texture_coords; }
	vector<Vector3D> get_normals() { return _normals; }
	BoundingVolume* get_bounding_vol() { return _bounding_vol; }
	vector<string> split(const char * s, char delim);

private:
	vector<Point3D> _vertices;
	vector<Vector3D> _normals;
	vector<Point3D> _texture_coords;
	vector< vector<int> > _face_vertices;
	vector< vector<int> > _face_textures;
	vector< vector<int> > _face_normals;
	BoundingVolume* _bounding_vol;
};
