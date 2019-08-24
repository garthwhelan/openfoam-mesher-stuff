#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<boost/format.hpp>
#include<algorithm>
#include<cmath>

class point {
public:
  float p[3];
  bool operator==(const point& rhs);
  point(float x, float y, float z);
};

class face {
public:
  std::vector<int> point_inds;
  int owner;
  int neighbour;
  //==, < are not be intuitive
  bool operator==(const face& rhs);
  face();
  face(int owner, int neighbour, std::vector<int> point_inds);
  void combine_faces(face& f);
};

enum PATCH_TYPE {PATCH,EMPTY,SYMMETRYPLANE,WALL,WEDGE,CYCLIC,PROCESSOR};
struct patch {
  std::vector<int> faces;
  std::string name;
  PATCH_TYPE pt;
};

class Mesh {
public:
  std::vector<point> points;
  std::vector<face> faces;
  std::vector<patch> patches;
  int ncells;
  void remove_duplicate_points();
  void remove_duplicate_faces();
  void write_mesh();
  void order_mesh();

};

Mesh make_cartesian_mesh(int xdim, int ydim, int zdim);






