#include "point.hpp"
#include "mesh.hpp"

/*
  all mesh qualities to be calculated by point index, 
  and to be focused on things in terms of faces
  lower values are higher quality
 */

double Mesh::flow_alignment(int point_ind, std::vector<int> face_inds, point (*v_dir)(point)) {
  //12 faces for hexahedral mesh, expect 4 to be aligned, 8 to be perpendicular
  double total_flow_alignment = 0;
  for(int face_ind : face_inds) {
    total_flow_alignment+=this->face_flow_alignment(face_ind,v_dir);
  }
  return total_flow_alignment;
}

double Mesh::face_flow_alignment(int face_ind,point (*v_dir)(point)) {
  point f_CoM = this->face_CoM(face_ind);
  double A = 1.0-point::dot(this->face_normal(face_ind),v_dir(f_CoM));
  double B = 1.0-point::len(point::cross(this->face_normal(face_ind),v_dir(f_CoM)));
  if (A<B) { return A; } else { return B; }

}

double Mesh::skewness(int point_ind, std::vector<int> face_inds) {
  double total_skewness = 0;//lower is better
  for(int face_ind : face_inds) {
    total_skewness+=this->face_skewness(face_ind);
  }
  return total_skewness;
}

double Mesh::face_skewness(int face_ind) {
  double max_cos_theta = 0;
  face f = this->faces[face_ind];
  for(int i = 0; i < f.point_inds.size(); i++) {
    point vec1=this->points[f.point_inds[(0+i)%f.point_inds.size()]]-this->points[f.point_inds[(1+i)%f.point_inds.size()]];
    point vec2=this->points[f.point_inds[(1+i)%f.point_inds.size()]]-this->points[f.point_inds[(2+i)%f.point_inds.size()]];
    double cos_theta=point::dot(vec1,vec2)/(point::len(vec1)*point::len(vec2));
    if(cos_theta>max_cos_theta) max_cos_theta=cos_theta;
  }
  return max_cos_theta;
}

// //also compare face normal direction to this
// double Mesh::nonorthogonality_to_neighbour_cells(int ind) {
//   cell c = this->cells[ind];
//   double max_nonortho = 0.0;
//   for(int own : c.owns) {
//     point vec1 = this->cell_CoM(ind)-this->face_CoM(own);
//     point vec2 = this->cell_CoM(ind)-this->cell_CoM(this->faces[own].neighbour);
//     double res = 1-point::dot(vec1,vec2)/pow(point::dot(vec1,vec1)*point::dot(vec2,vec2),0.5);
//     if(res > max_nonortho) max_nonortho=res;
//   }
//   for(int neigh : c.neighbours) {
//     point vec1 = this->cell_CoM(ind)-this->face_CoM(neigh);
//     point vec2 = this->cell_CoM(ind)-this->cell_CoM(this->faces[neigh].owner);
//     double res = 1-point::dot(vec1,vec2)/pow(point::dot(vec1,vec1)*point::dot(vec2,vec2),0.5);
//     if(res > max_nonortho) max_nonortho=res;
//   }
//   return max_nonortho;
// }


// double Mesh::face_aspect_ratio(int ind) {//also should be improved
//   double res1=point::len(this->points[this->faces[ind].point_inds[0]]-this->points[this->faces[ind].point_inds[1]]);
//   double res2=point::len(this->points[this->faces[ind].point_inds[0]]-this->points[this->faces[ind].point_inds[2]]);
//   if(res1>res2) {
//     return res1;
//   } else {
//     return res2;
//   }
// }
