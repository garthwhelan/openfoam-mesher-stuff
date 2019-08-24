#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<boost/format.hpp>
#include<algorithm>
#include "mesh.hpp"
#include<map>
#include<utility> //for pair
//TODO: other patch types
//fix order_mesh
const float small_number = 1e-4;

int gen_ind(int i, int j, int k, int nx,int ny, int nz) {
  return i*(ny*nz)+j*(nz)+k;
}

bool near(float a, float b) {
  return (std::abs(a-b)<small_number);
}

bool point::operator==(const point& rhs) {
  return near(this->p[0],rhs.p[0])&&near(this->p[1],rhs.p[1])&&near(this->p[2],rhs.p[2]);
}

bool face::operator==(const face& rhs) {
  return std::is_permutation(this->point_inds.begin(),this->point_inds.end(),rhs.point_inds.begin());
}

point::point(float x, float y, float z) {
  p[0]=x;
  p[1]=y;
  p[2]=z;
}

face::face() {
}

face::face(int owner, int neighbour, std::vector<int> point_inds) {
  this->owner=owner;
  this->neighbour=neighbour;
  this->point_inds=point_inds;
}

bool face_order(face f1, face f2) {
  return f1.owner>f2.owner;
}

//TODO: refactor/make so checkMesh doesn't complain about
//upper triangularity-should probably just reindex cells
//workaround: use renumberMesh shell utility from openFoam
void Mesh::order_mesh() {
  std::vector<face> internal_faces;
  for(face f : this->faces) {
    if(f.neighbour!=-1) internal_faces.push_back(f);
  }
  std::sort(internal_faces.begin(),internal_faces.end(),face_order);
          
  std::vector<std::vector<face>> patch_face_lists;
  for(patch p : this->patches) {
    std::vector<face> face_list;
    for(int face_ind : p.faces) {
      face_list.push_back(this->faces[face_ind]);
    }
    std::sort(face_list.begin(),face_list.end(),face_order);
    patch_face_lists.push_back(face_list);
  }
  
  this->faces.clear();
  this->faces.insert(this->faces.begin(),internal_faces.begin(),internal_faces.end());

  for(int i = 0; i < this->patches.size(); i++) {
    std::vector<face> face_list = patch_face_lists[i];
     int initial_size = this->faces.size();
    this->faces.insert(this->faces.end(),face_list.begin(),face_list.end());
    int end_size = this->faces.size();
    this->patches[i].faces.clear();
    for(int j = initial_size; j < end_size; j++) this->patches[i].faces.push_back(j);
  }
}

//for combining meshes or patches
void Mesh::remove_duplicate_points() {  
  std::vector<point> point_array;
  for(int i = 0; i < this->points.size(); i++) {
    int p_loc = std::distance(point_array.begin(),std::find(point_array.begin(),point_array.end(),this->points[i]));
    if(p_loc != point_array.size()) {
      for(face &f : this->faces) {
        std::replace(f.point_inds.begin(),f.point_inds.end(),i,p_loc);
      }
    } else {
      point_array.push_back(this->points[i]);
      for(face &f : this->faces) {
        std::replace(f.point_inds.begin(),f.point_inds.end(),i,(int)point_array.size()-1);
      }
    }
  }
  this->points.clear();
  this->points=point_array;
}

void face::combine_faces(face& f) {
  if(this->owner>f.owner) {
    this->neighbour=f.owner;
  } else {
    this->neighbour = this->owner;
    this->owner = f.owner;
    std::reverse(this->point_inds.begin(),this->point_inds.end());
  }
  f.owner=-1;
}

//only works when duplicate faces are in patches
void Mesh::remove_duplicate_faces() {
  //search through the patches for duplicate faces
  for(int i = 0; i < this->faces.size(); i++) {
    if(this->faces[i].neighbour==-1) {
      int f_ind = std::distance(this->faces.begin(),std::find(this->faces.begin()+i+1,this->faces.end(),this->faces[i]));
      if(f_ind!=this->faces.size()) { //pair of faces that are the same
        for(patch &p : this->patches) {//combine so remove both from patches
          p.faces.erase(std::remove_if(p.faces.begin(),p.faces.end(),[i,f_ind](int ind){return ((ind==i)or(ind==f_ind));}),p.faces.end());
        }
        //get orientation right and mark second one for deletion
        this->faces[i].combine_faces(this->faces[f_ind]);
      }
    }
  }
  //build up array of unique faces and reindex patches
  std::vector<face> face_array;
  for(int i = 0; i < this->faces.size(); i++) {
    if(this->faces[i].owner!=-1) {
      face_array.push_back(this->faces[i]);
      for(patch &p : this->patches) {
        for(int &find : p.faces) {
          if(find==i) find=face_array.size()-1;
        }
      }
    }
  }
  this->faces.clear();
  this->faces=face_array;
}

  

// 
  // for(int i = 0; i < this->faces.size(); i++) {

  //   bool in_f_arr = (std::find(face_array.begin(),face_array.end(),this->faces[i])!=face_array.end());
  //   if(in_f_arr) continue;
    
  //   face_array.push_back(faces[i]);
    

    
Mesh make_cartesian_mesh(int xdim, int ydim, int zdim) {
  
  Mesh M;

  M.ncells=xdim*ydim*zdim;
  
  for(int i = 0; i < xdim+1; i++) {
    for(int j = 0; j < ydim+1; j++) {
      for(int k = 0; k < zdim+1; k++) {
        (M.points).push_back(point(i,j,k));
      }
    }
  }
  
  std::vector<face> xfaces((xdim+1)*ydim*zdim);
  for(int i = 0; i < xdim+1; i++) {
    for(int j = 0; j < ydim; j++) {
      for(int k = 0; k < zdim; k++) {
        face A;
        (A.point_inds).push_back(gen_ind(i,j,k,xdim+1,ydim+1,zdim+1));
        (A.point_inds).push_back(gen_ind(i,j+1,k,xdim+1,ydim+1,zdim+1));
        (A.point_inds).push_back(gen_ind(i,j+1,k+1,xdim+1,ydim+1,zdim+1));
        (A.point_inds).push_back(gen_ind(i,j,k+1,xdim+1,ydim+1,zdim+1));
        if(i==0) {
          A.owner=gen_ind(i,j,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else if(i==xdim){
          std::reverse(A.point_inds.begin(),A.point_inds.end());
          A.owner=gen_ind(i-1,j,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else {
          std::reverse(A.point_inds.begin(),A.point_inds.end());
          A.owner=gen_ind(i-1,j,k,xdim,ydim,zdim);
          A.neighbour=gen_ind(i,j,k,xdim,ydim,zdim);
        }
        xfaces[gen_ind(i,j,k,xdim+1,ydim,zdim)]=A;
      }}}
  (M.faces).insert(M.faces.end(),xfaces.begin(),xfaces.end());
    
  std::vector<face> yfaces(xdim*(ydim+1)*zdim);
  for(int i = 0; i < xdim; i++) {
    for(int j = 0; j < ydim+1; j++) {
      for(int k = 0; k < zdim; k++) {
        face A;
        (A.point_inds).push_back(gen_ind(i,j,k,xdim+1,ydim+1,zdim+1));
        (A.point_inds).push_back(gen_ind(i,j,k+1,xdim+1,ydim+1,zdim+1));
        (A.point_inds).push_back(gen_ind(i+1,j,k+1,xdim+1,ydim+1,zdim+1));
        (A.point_inds).push_back(gen_ind(i+1,j,k,xdim+1,ydim+1,zdim+1));
        if(j==0) {          
          A.owner=gen_ind(i,j,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else if(j==ydim){
          std::reverse(A.point_inds.begin(),A.point_inds.end());
          A.owner=gen_ind(i,j-1,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else {
          std::reverse(A.point_inds.begin(),A.point_inds.end());
          A.owner=gen_ind(i,j-1,k,xdim,ydim,zdim);
          A.neighbour=gen_ind(i,j,k,xdim,ydim,zdim);
        }
        yfaces[gen_ind(i,j,k,xdim,ydim+1,zdim)]=A;
      }}}
  (M.faces).insert(M.faces.end(),yfaces.begin(),yfaces.end());
  
  std::vector<face> zfaces(xdim*ydim*(zdim+1));
  for(int i = 0; i < xdim; i++) {
    for(int j = 0; j < ydim; j++) {
      for(int k = 0; k < zdim+1; k++) {
        face A;
        (A.point_inds).push_back(gen_ind(i,j,k,xdim+1,ydim+1,zdim+1));
        (A.point_inds).push_back(gen_ind(i+1,j,k,xdim+1,ydim+1,zdim+1));
        (A.point_inds).push_back(gen_ind(i+1,j+1,k,xdim+1,ydim+1,zdim+1));
        (A.point_inds).push_back(gen_ind(i,j+1,k,xdim+1,ydim+1,zdim+1));
        if(k==0) {
          A.owner=gen_ind(i,j,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else if(k==zdim){
          std::reverse(A.point_inds.begin(),A.point_inds.end());
          A.owner=gen_ind(i,j,k-1,xdim,ydim,zdim);
          A.neighbour=-1;
        } else {
          std::reverse(A.point_inds.begin(),A.point_inds.end());
          A.owner=gen_ind(i,j,k-1,xdim,ydim,zdim);
          A.neighbour=gen_ind(i,j,k,xdim,ydim,zdim);
        }
        zfaces[gen_ind(i,j,k,xdim,ydim,zdim+1)]=A;
      }}}
  (M.faces).insert(M.faces.end(),zfaces.begin(),zfaces.end());

  //face order is all (nx+1)*ydim*zdim xfaces then nx*(ydim+1)*zdim yfaces, then nx*ydim*(zdim+1) zfaces
  patch xu_patch;  
  patch xd_patch;
  xu_patch.name="xu";
  xd_patch.name="xd";
  xu_patch.pt=PATCH;
  xd_patch.pt=PATCH;
  for(int j = 0; j < ydim; j++) {
    for(int k = 0; k < zdim; k++) {
      xu_patch.faces.push_back(gen_ind(0,j,k,xdim+1,ydim,zdim));
      xd_patch.faces.push_back(gen_ind(xdim,j,k,xdim+1,ydim,zdim));
    }}
  M.patches.push_back(xu_patch);
  M.patches.push_back(xd_patch);

  patch yu_patch;
  patch yd_patch;
  yu_patch.name="yu";
  yd_patch.name="yd";
  yu_patch.pt=PATCH;
  yd_patch.pt=PATCH;
  for(int i = 0; i < xdim; i++) {
    for(int k = 0; k < zdim; k++) {
      yu_patch.faces.push_back((xdim+1)*ydim*zdim+gen_ind(i,0,k,xdim,ydim+1,zdim));
      yd_patch.faces.push_back((xdim+1)*ydim*zdim+gen_ind(i,ydim,k,xdim,ydim+1,zdim));
    }}
  M.patches.push_back(yu_patch);
  M.patches.push_back(yd_patch);

  patch zu_patch;
  patch zd_patch;
  zu_patch.name="zu";
  zd_patch.name="zd";
  zu_patch.pt=PATCH;
  zd_patch.pt=PATCH;
  for(int i = 0; i < xdim; i++) {
    for(int j = 0; j < ydim; j++) {
      zu_patch.faces.push_back((xdim+1)*ydim*zdim+xdim*(ydim+1)*zdim+gen_ind(i,j,0,xdim,ydim,zdim+1));
      zd_patch.faces.push_back((xdim+1)*ydim*zdim+xdim*(ydim+1)*zdim+gen_ind(i,j,zdim,xdim,ydim,zdim+1));
    }}
  M.patches.push_back(zu_patch);
  M.patches.push_back(zd_patch);
  
  //order_mesh(M);
  
  return M;
}

void Mesh::write_mesh() {
  std::string line;
  //points
  std::ifstream points_template("outputfile_templates/points_template");
  std::ofstream points_out("polyMesh/points");
  while(getline(points_template,line)) {
    points_out << line << '\n';
  }
  points_out << boost::format("%d\n(\n") % (int)this->points.size();
  for(int ind = 0; ind < this->points.size(); ind++) {
    points_out << boost::format("(%f %f %f)\n") % this->points[ind].p[0] % this->points[ind].p[1] % this->points[ind].p[2];
  }
  points_out << ")\n";
  points_template.close();
  points_out.close();

  std::ifstream faces_template("outputfile_templates/faces_template");
  std::ifstream owner_template("outputfile_templates/owner_template");
  std::ifstream neighbour_template("outputfile_templates/neighbour_template");
  std::ifstream boundary_template("outputfile_templates/boundary_template");

  std::ofstream faces_out("polyMesh/faces");
  std::ofstream owner_out("polyMesh/owner");
  std::ofstream neighbour_out("polyMesh/neighbour");
  std::ofstream boundary_out("polyMesh/boundary");

  while(getline(faces_template,line)) {
    faces_out << line << '\n';
  }

  int internalfaces = this->faces.size();
  for(patch A : this->patches) {
    internalfaces -= A.faces.size();
  }
  
  for(int i = 0; i < 17; i++) {
    getline(owner_template,line);
    if(i==12) {
      owner_out << boost::format("\t\tnote        \"nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d\";\n") % (int)this->points.size() % this->ncells % (int)this->faces.size() % internalfaces;
    } else {
      owner_out << line << "\n";
    }
    getline(neighbour_template,line);
    if(i==12) {
      neighbour_out << boost::format("\t\tnote        \"nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d\";\n") % (int)this->points.size() % this->ncells % (int)this->faces.size() % internalfaces;
    } else {
      neighbour_out << line << "\n";
    }
  }
  
  faces_out << boost::format("%d\n(\n") % (int)this->faces.size();
  owner_out << boost::format("%d\n(\n") % (int)this->faces.size();
  neighbour_out << boost::format("%d\n(\n") % internalfaces;
    
  while(getline(boundary_template,line)) {
    boundary_out << line << '\n';
  }
    
  int start_face = 0;
    
  for(int ind = 0; ind < this->faces.size(); ind++) {
    face A = this->faces[ind];
    if(A.neighbour!=-1) {
      faces_out << boost::format("4(%d %d %d %d)\n") % A.point_inds[0] % A.point_inds[1] % A.point_inds[2] % A.point_inds[3];
      owner_out << boost::format("%d\n") % A.owner;
      neighbour_out << boost::format("%d\n") % A.neighbour;
      start_face++;
    }
  }

  boundary_out << boost::format("%d\n(\n") % this->patches.size();
  int patch_count = 1;
  for(patch P : this->patches) {

    switch(P.pt) {
    case PATCH:
      boundary_out << boost::format("\t%s\n\t{\n\t\ttype patch;\n\t\tnFaces %d;\n\t\tstartFace %d;\n\t}\n") %P.name % P.faces.size() % start_face;
      break;
    case EMPTY:
      boundary_out << boost::format("\t%s\n\t{\n\t\ttype empty;\n\t\tinGroups 1(empty);\n\t\tnFaces %d;\n\t\tstartFace %d;\n\t}\n") % P.name % P.faces.size() % start_face;
      break;
    case WALL:
      boundary_out << boost::format("\t%s\n\t{\n\t\ttype wall;\n\t\tinGroups 1(wall);\n\t\tnFaces %d;\n\t\tstartFace %d;\n\t}\n") % P.name % P.faces.size() % start_face;
      break;
    default:
      printf("Other patch types not implemented\n");
      break;
    }

    start_face+=P.faces.size();
    patch_count++;
    for(int ind : P.faces) {
      faces_out << boost::format("4(%d %d %d %d)\n") % this->faces[ind].point_inds[0] % this->faces[ind].point_inds[1] % this->faces[ind].point_inds[2] % this->faces[ind].point_inds[3];
      owner_out << boost::format("%d\n") % this->faces[ind].owner;
    }
  }
  boundary_out << ")";
  faces_out << ")";
  owner_out << ")";
  neighbour_out << ")";
    
  faces_template.close();
  owner_template.close();
  neighbour_template.close();
  boundary_template.close();

  faces_out.close();
  owner_out.close();
  neighbour_out.close();
  boundary_out.close();  
}

//break patch into 2 subpatches
/*void partition_patch(Mesh M) {
  
  }*/

//combine 2 patches with different faces but the same points
//void merge_patches(Mesh M, patch1, patch2)
