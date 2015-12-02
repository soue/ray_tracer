#include "mesh.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstdio>


vector<string> Mesh::split(const char * s, const char delim) {
	vector<string> result;
	stringstream sstream(s);
	string token;	
	string newline("\n");
	while (getline(sstream, token, delim)) {
		if (!token.empty() && token.compare(newline) != 0) {
			result.push_back(token);

		}
	}
	return result;
}

Mesh::Mesh(const char * file_name) {
	FILE * f;
	char buffer[BUFFER];
	Point3D p;
	double xmin = INFINITY;
	double ymin = INFINITY;
	double zmin = INFINITY;
	double xmax = -INFINITY;
	double ymax = -INFINITY;
	double zmax = -INFINITY;
	if ((f = fopen(file_name, "r")) == NULL) {
		std::cout << "Can't open file " << file_name << "\n";
	} else {
		while (fgets(buffer, sizeof(buffer), f) != NULL) {
			if (strlen(buffer) == 0 || buffer[0] == '#') {
				continue;
			}
			vector<string> line = split(buffer, ' ');

			if (line.size() > 0) {
				if (line[0] == "v") {
					p = Point3D(atof(line[1].c_str()),
								atof(line[2].c_str()),
								atof(line[3].c_str()));
					xmin = fmin(xmin,p[0]);
					ymin = fmin(ymin,p[1]);
					zmin = fmin(zmin,p[2]);
					xmax = fmax(xmax,p[0]);
					ymax = fmax(ymax,p[1]);
					zmax = fmax(zmax,p[2]);
				
					_vertices.push_back(p);
				} else if (line[0] == "f") {
					vector<int> vertices;
					vector<int> textures;
					vector<int> normals;
					for (int f = 1; f < line.size(); f++) {
						int index = line[f].find("/");
						if (index > 0) {
							vertices.push_back(atoi(line[f].substr(0,index).c_str()));
							string next = line[f].substr(index+1,line[f].length());
							int index2 = next.find("/");
							if (index2 < 0) {
								textures.push_back(atoi(next.c_str()));
								normals.push_back(0);
							} else if (index2 > 0) {
								textures.push_back(atoi(next.substr(0,index2).c_str()));
								normals.push_back(atoi(next.substr(index2+1,next.length()).c_str()));
							} else {
								normals.push_back(atoi(next.substr(index2+1,next.length()).c_str()));
								textures.push_back(0);
							}
							
						} else {
							vertices.push_back(atoi(line[f].c_str()));
						}
					}
					_face_vertices.push_back(vertices);
					_face_textures.push_back(textures);
					_face_normals.push_back(normals);
					
				} else if (line[0] == "vn") {
					Vector3D normal = Vector3D(atof(line[1].c_str()),
												atof(line[2].c_str()),
												atof(line[3].c_str()));
					normal.normalize();
					_normals.push_back(normal);
				} else if (line[0] == "vt") {
					_texture_coords.push_back(Point3D(atof(line[1].c_str()),
													  atof(line[2].c_str()),
													  0));
				}	
			}
		}
		_bounding_vol = new BoundingVolume(Point3D(xmin,ymin,zmin),
									   	   Point3D(xmax,ymax,zmax));
	}
	fclose(f);
}

