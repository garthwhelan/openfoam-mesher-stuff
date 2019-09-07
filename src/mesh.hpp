#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<boost/format.hpp>
#include<algorithm>
#include<cmath>

/*
  mesh has points, faces(w/owners and neighbors), and patches (of faces + additional info)
  some rules: each face only in at most one patch
  would like to be able to index faces by cells, possibly by both owners and neighbors
  all faces are unique

  map from faces to owners and neighbours?

  for each point
  for each face (and points associated)
  for each cell (and faces and points associated)
 */

struct point {
  float x,y,z;
  bool operator==(const point& rhs);
  point(float xi, float yi, float zi) : x{xi},y{yi},z{zi} {}
};

class face {
public:
  std::vector<int> point_inds;
  int owner;
  int neighbour;
  //==, < are not be intuitive
  bool operator==(const face& rhs);
  face() : owner{-2},neighbour{-2},point_inds{} {}
  face(int o, int n, std::vector<int> pinds) : owner{o},neighbour{n},point_inds{pinds} {}
  void combine_faces(face& f);
};

enum class PT {PATCH,EMPTY,SYMMETRYPLANE,WALL,WEDGE,CYCLIC,PROCESSOR};
struct patch {
  std::vector<int> faces;
  std::string name;
  PT pt;//patch type
};

class Mesh {
public:
  std::vector<point> points;
  std::vector<face> faces;
  std::vector<patch> patches;
  Mesh() : points{},faces{},patches{} {}
  void remove_duplicate_points();
  void remove_duplicate_faces();
  void write_mesh();
  void order_mesh();

  void remove_intersecting_blocks(bool (*fn) (point));
  //map to boundary for cartesian_mesh
  
  int ncells() const;
  static Mesh combine_meshes(const Mesh& M1, const Mesh& M2);
  static Mesh make_3D_cartesian_mesh(int xdim, int ydim, int zdim);
  static Mesh make_2D_cartesian_mesh(int xdim, int ydim);
};








