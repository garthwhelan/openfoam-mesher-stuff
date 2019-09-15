#pragma once

#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<boost/format.hpp>
#include<algorithm>
#include<cmath>
#include "point.hpp"

/*
https://cfd.direct/openfoam/user-guide/v7-mesh-description/
 */



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

struct cell {
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

  bool faceind_in_patch(int);
  bool pointind_in_mesh_faces(int);

  void remove_invalid_faces();
  void remove_invalid_cells();
  void correctly_orient_faces();
  void cleanup_boundary();
  void cleanup();

  point face_CoM(int);
  point cell_CoM(int);
  point face_normal(int);
  
  void remove_duplicate_points();
  void remove_duplicate_faces();

  void write_mesh();

  static Mesh make_3D_cartesian_mesh(int xdim, int ydim, int zdim);
  static Mesh make_2D_cartesian_mesh(int xdim, int ydim);
  static Mesh combine_meshes(const Mesh&, const Mesh&);
  void remove_cells(std::vector<int>,patch&);

};
