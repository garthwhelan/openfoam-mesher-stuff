#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<boost/format.hpp>
#include<algorithm>
#include<cmath>

#define MAX_POINTS_IN_FACE

/*
https://cfd.direct/openfoam/user-guide/v7-mesh-description/
 */
  
struct point {
  float x,y,z;
  bool operator==(const point& rhs);
  point(float x, float y, float z) : x{x},y{y},z{z} {}
  void print();
};

class face {
public:
  std::vector<int> point_inds;
  //this owner/neighbour only to be used to simplify loops over faces
  int owner,neighbour;
  //==, < are not intuitive
  bool operator==(const face& rhs);
  face() : point_inds{},owner{-2},neighbour{-1} {}
  face(std::vector<int> point_inds) : owner{-2},neighbour{-1},point_inds{point_inds} {}
  face(int owner, int neighbour, std::vector<int> point_inds) : owner{owner},neighbour{neighbour},point_inds{point_inds} {}
  void combine_faces(face& f);
  float area();
};

enum class PT {PATCH,EMPTY,SYMMETRYPLANE,WALL,WEDGE,CYCLIC,PROCESSOR};

class patch {
private:
  static int patch_count;//to give patches unique names
public:

  std::vector<int> face_inds;
  std::string patch_name;
  PT patch_type;
  patch() : face_inds{},patch_name{"patch_"+std::to_string(patch_count)},patch_type{PT::PATCH} {patch_count++;}
  patch(std::vector<int> face_inds, std::string patch_name, PT patch_type) : face_inds{face_inds},patch_name{patch_name},patch_type{patch_type} {};
};

class cell {
public:
  //faces each have exactly one owner and at most one neighbour
  std::vector<int> owns;
  std::vector<int> neighbours;
  cell() : owns{},neighbours{} {}
};

class Mesh {
public:
  std::vector<point> points;
  std::vector<face> faces;
  std::vector<cell> cells;
  std::vector<patch> patches;
  Mesh() : points{},faces{},cells{},patches{} {}
  int ncells() const;
  //also removes unused points
  void remove_duplicate_points();
  void remove_duplicate_faces();
  void remove_cells(std::vector<int>,patch&);
  void write_mesh();
  bool faceind_in_patch(int);
  bool pointind_in_mesh_faces(int);
  void cleanup();

  // void order_mesh();

  // void remove_intersecting_blocks(bool (*fn) (point));
  // //map to boundary for cartesian_mesh
  
  static Mesh combine_meshes(const Mesh& M1, const Mesh& M2);
  static Mesh make_3D_cartesian_mesh(int xdim, int ydim, int zdim);
  static Mesh make_2D_cartesian_mesh(int xdim, int ydim);
};
