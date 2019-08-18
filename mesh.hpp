#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<boost/format.hpp>
#include<algorithm>
#include<cmath>

struct point {
  float p[3];
};

struct face {
  std::vector<int> point_inds;
  int owner;
  int neighbour;
};

enum PATCH_TYPE {PATCH,EMPTY,SYMMETRYPLANE,WALL,WEDGE,CYCLIC,PROCESSOR};
struct patch {
  std::vector<int> faces;
  std::string name;
  PATCH_TYPE pt;
};

struct Mesh {
  std::vector<point> points;
  std::vector<face> faces;
  std::vector<patch> patches;
  int ncells;
};

Mesh make_cartesian_mesh(int xdim, int ydim, int zdim);

void write_mesh(Mesh M);


