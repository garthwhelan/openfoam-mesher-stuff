#include "point.hpp"
#include "mesh.hpp"

bool Mesh::is_pointind_used(int i) {
  for(face f : this->faces) {
    for(int p_ind : f.point_inds) {
      if(i==p_ind) return true;
    }
  }
  return false;
}

bool Mesh::pointind_in_boundary(int ind) {
  for(patch p : this->patches) {
    if(p.patch_name=="frontAndBackPlanes") break;
    for(int face_ind : p.face_inds) {
      for(int point_ind : this->faces[face_ind].point_inds) {
        if(ind == point_ind) return true;
      }
    }
  }
  return false;
}

point Mesh::face_normal(int ind) {
  face f = this->faces[ind];
  point p = point::cross(this->points[f.point_inds[1]]-this->points[f.point_inds[0]],
                         this->points[f.point_inds[2]]-this->points[f.point_inds[1]]);
  return p/point::len(p);
}

point Mesh::face_CoM(int ind) {
  point p {0,0,0};
  int count = 0;
  face f=this->faces[ind];
  for(int point_ind : f.point_inds) {
    p+=this->points[point_ind];
    count++;
  }
  p/=count;
  return p;
}

point Mesh::cell_CoM(int ind) {
  point p {0,0,0};
  int count = 0;
  for(int face_ind : this->cells[ind].owns) {
    p+=face_CoM(face_ind);
    count++;
  }
  for(int face_ind : this->cells[ind].neighbours) {
    p+=face_CoM(face_ind);
    count++;
  }
  p/=count;
  return p;
}

std::vector<int> Mesh::get_faces_by_point(int ind) {
  std::vector<int> ret_val {};
  for(int i = 0; i < this->faces.size(); i++) {
    for(int point_ind : this->faces[i].point_inds) {
      if(point_ind==ind) {
        ret_val.push_back(i);
        break;
      }
    }
  }
  return ret_val;
}

