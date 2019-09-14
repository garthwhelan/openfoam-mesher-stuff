#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<boost/format.hpp>
#include<algorithm>
#include "mesh.hpp"
//TODO: other patch types
//fix order_mesh
const float small_number = 1e-6;
int patch::patch_count = 0;


void point::print() {
  std::cout << "[" << this->x << "," << this->y << ","<< this->z << "]";
}

int gen_ind(int i, int j, int k, int nx,int ny, int nz) {
  return i*(ny*nz)+j*(nz)+k;
}

bool near(float a, float b) {
  return (std::abs(a-b)<small_number);
}

bool point::operator==(const point& rhs) {
  return near(this->x,rhs.x)&&near(this->y,rhs.y)&&near(this->z,rhs.z);
}

bool face::operator==(const face& rhs) {
  return std::is_permutation(this->point_inds.begin(),this->point_inds.end(),rhs.point_inds.begin());
}

template<class T>
void concat_arrs(T& out, const T& in1, const T& in2) {
  out.reserve(in1.size()+in2.size());
  out.insert(out.end(),in1.begin(),in1.end());
  out.insert(out.end(),in2.begin(),in2.end());
}

int Mesh::ncells() const {
  return this->cells.size();
}

Mesh Mesh::combine_meshes(const Mesh& M1, const Mesh& M2) {
  Mesh M;
  concat_arrs<std::vector<point>>(M.points,M1.points,M2.points);
  concat_arrs<std::vector<face>>(M.faces,M1.faces,M2.faces);
  int M1points = M1.points.size();
  int M1ncells = M1.ncells();
  std::transform(M.faces.begin()+M1.faces.size(),M.faces.end(),M.faces.begin()+M1.faces.size(),
                 [M1points,M1ncells](face f){
                   for(int &pind : f.point_inds) {
                     pind+=M1points;
                   }
                   f.owner+=M1ncells;
                   if(f.neighbour!=-1) f.neighbour+=M1ncells;
                   return f;
                 });
  M.remove_duplicate_points();

  for(int i = 0; i < M1.ncells()+M2.ncells(); i++) M.cells.push_back(cell());
    
  for(int i = 0; i < M.faces.size(); i++) {
    M.cells[M.faces[i].owner].owns.push_back(i);
    if(M.faces[i].neighbour!=-1)
      M.cells[M.faces[i].neighbour].neighbours.push_back(i);
  }
  
  int M1faces = M1.faces.size();
  concat_arrs<std::vector<patch>>(M.patches,M1.patches,M2.patches);
  std::transform(M.patches.begin(),M.patches.begin()+M1.patches.size(),M.patches.begin(),
                  [](patch p){
                    p.patch_name+="1";//so patch names are unique
                    return p;});
  std::transform(M.patches.begin()+M1.patches.size(),M.patches.end(),M.patches.begin()+M1.patches.size(),
                 [M1faces](patch p){
                   for(int &find : p.face_inds) {
                     find+=M1faces;
                   }
                   p.patch_name+="2";
                   return p;});
  M.remove_duplicate_faces();
  M.cleanup();
  return M;
}

bool Mesh::pointind_in_mesh_faces(int i) {
  for(face f : this->faces) {
    for(int p_ind : f.point_inds) {
      if(i==p_ind) return true;
    }
  }
  return false;
}

//for combining meshes or patches
void Mesh::remove_duplicate_points() {  
  std::vector<point> point_array;
  for(int i = 0; i < this->points.size(); i++) {
    if(not pointind_in_mesh_faces(i)) continue;
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

//only those duplicate faces in patches
void Mesh::remove_duplicate_faces() {
  //search through the patches for duplicate faces
  for(int i = 0; i < this->faces.size(); i++) {
    if(this->faces[i].neighbour==-1) {
      int f_ind = std::distance(this->faces.begin(),std::find(this->faces.begin()+i+1,this->faces.end(),this->faces[i]));
      if(f_ind!=this->faces.size()) { //pair of faces that are the same
        for(patch &p : this->patches) {//combine so remove both from patches
          p.face_inds.erase(std::remove_if(p.face_inds.begin(),p.face_inds.end(),[i,f_ind](int ind){return ((ind==i)or(ind==f_ind));}),p.face_inds.end());
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
        for(int &find : p.face_inds) {
          if(find==i) find=face_array.size()-1;
        }
      }
    }
  }
  this->faces.clear();
  this->faces=face_array;
  this->patches.erase(std::remove_if(this->patches.begin(),this->patches.end(),[](patch p){return p.face_inds.size()==0;}),this->patches.end());
  //redo the inds for cells
  for(cell &c : this->cells) {
    c.owns.clear();
    c.neighbours.clear();
  }
  for(int i = 0; i < this->faces.size(); i++) {
    int owner = this->faces[i].owner;
    int neighbour = this->faces[i].neighbour;
    this->cells[owner].owns.push_back(i);
    if(neighbour!=-1) this->cells[neighbour].neighbours.push_back(i);
  }
}

bool Mesh::faceind_in_patch(int f_ind) {
  for(patch p : this->patches) {
    if(std::find(p.face_inds.begin(),p.face_inds.end(),f_ind)!=p.face_inds.end())
      return true;
  }
  return false;
}

void Mesh::cleanup() {
  //remove faces which are unaccessable
  //remove cells with no faces
  //make sure all indices are contiguous
  //reindex patches
  //put faces outside of patches into a new patch
  ////////////
  std::vector<face> valid_faces;  
  for(int i = 0; i < this->faces.size(); i++) {

    if((this->faces[i].owner==-1)and(this->faces[i].neighbour!=-1)) {
      this->faces[i].owner=this->faces[i].neighbour;
      this->faces[i].neighbour=-1;
      std::reverse(this->faces[i].point_inds.begin(),this->faces[i].point_inds.end());
    }
    
    if(this->faces[i].owner==-1) {
      for(patch &p : this->patches) {
        p.face_inds.erase(std::remove_if(p.face_inds.begin(),p.face_inds.end(),[i](int ind){return ind==i;}),p.face_inds.end());
      }
      for(cell &c : this->cells) {
        c.owns.erase(std::remove_if(c.owns.begin(),c.owns.end(),[i](int ind){return ind==i;}),c.owns.end());
        c.neighbours.erase(std::remove_if(c.neighbours.begin(),c.neighbours.end(),[i](int ind){return ind==i;}),c.neighbours.end());
      }
    } else {//!=-1    
      valid_faces.push_back(this->faces[i]);
      for(patch &p : this->patches) {
        for(int &f_ind : p.face_inds) {
          if(f_ind==i) f_ind=valid_faces.size()-1;
        }}
      for(cell &c : this->cells) {
        for(int &own : c.owns) {
          if(own==i) own=valid_faces.size()-1;
        }
        for(int &neigh : c.neighbours) {
          if(neigh==i) neigh=valid_faces.size()-1;
        }
      }      
    }
  }
  this->faces.clear();
  this->faces=valid_faces;

  //remove empty patches
  this->patches.erase(std::remove_if(this->patches.begin(),this->patches.end(),
      [](patch p){return p.face_inds.size()==0;}),this->patches.end());
  
  std::vector<cell> valid_cells;
  for(int i = 0; i < this->cells.size(); i++) {
    if(this->cells[i].owns.size()+this->cells[i].neighbours.size()==0) {
      //no face will have this as owner or neighbour
    } else {//valid cell
      valid_cells.push_back(this->cells[i]);
      for(face &f : this->faces) {
        if(f.owner==i) f.owner=valid_cells.size()-1;
        if(f.neighbour==i) f.neighbour=valid_cells.size()-1;
      }
    }
  }
  this->cells.clear();
  this->cells=valid_cells;

  //put faces outside of patches into them
  bool need_new_patch = false;
  for(int i = 0; i < this->faces.size(); i++) {
    if((this->faces[i].neighbour==-1)and(not(this->faceind_in_patch(i)))) {
      need_new_patch=true;
      break;
    }
  }
  if(need_new_patch) {
    this->patches.push_back(patch());
    for(int i = 0; i < this->faces.size(); i++) {
      if((this->faces[i].neighbour==-1)and(not(this->faceind_in_patch(i)))) {
        this->patches[this->patches.size()-1].face_inds.push_back(i);
      }
    }    
  }
}
//for face attached to neighbour, set owner/neighbour to -1
//remove if both are negative remove that face
void Mesh::remove_cells(std::vector<int> cell_inds, patch& p) {
  for(int cell_ind : cell_inds) {
    for(int own : this->cells[cell_ind].owns) {
      this->faces[own].owner=-1;
    }
    for(int neigh : this->cells[cell_ind].neighbours) {
      this->faces[neigh].neighbour=-1;
    }
    this->cells[cell_ind].owns.clear();
    this->cells[cell_ind].neighbours.clear();
  }
  this->cleanup();
  for(int i = 0; i < this->faces.size(); i++) {
    if((faces[i].neighbour==-1)and(not this->faceind_in_patch(i)))
      p.face_inds.push_back(i);
  }
}

Mesh Mesh::make_3D_cartesian_mesh(int xdim, int ydim, int zdim) {
  
  Mesh M;

  M.points.reserve((xdim+1)*(ydim+1)*(zdim+1));
  for(int i = 0; i < xdim+1; i++) {
    for(int j = 0; j < ydim+1; j++) {
      for(int k = 0; k < zdim+1; k++) {
        M.points.push_back(point(i,j,k));
      }}}

  M.cells.reserve(xdim*ydim*zdim);
  for(int i = 0; i < xdim; i++) {
    for(int j = 0; j < ydim; j++) {
      for(int k = 0; k < zdim; k++) {
        M.cells.push_back(cell());
      }}}

  for(int i = 0; i < xdim+1; i++) {
    for(int j = 0; j < ydim; j++) {
      for(int k = 0; k < zdim; k++) {        
        int owner,neighbour,ind;
        std::vector<int> point_inds {gen_ind(i,j,k,xdim+1,ydim+1,zdim+1),gen_ind(i,j+1,k,xdim+1,ydim+1,zdim+1),gen_ind(i,j+1,k+1,xdim+1,ydim+1,zdim+1),gen_ind(i,j,k+1,xdim+1,ydim+1,zdim+1)};
        if(i==0) {
          owner=gen_ind(i,j,k,xdim,ydim,zdim);
          neighbour=-1;
        } else if(i==xdim){
          std::reverse(point_inds.begin(),point_inds.end());
          owner=gen_ind(i-1,j,k,xdim,ydim,zdim);
          neighbour=-1;
        } else {
          std::reverse(point_inds.begin(),point_inds.end());
          owner=gen_ind(i-1,j,k,xdim,ydim,zdim);
          neighbour=gen_ind(i,j,k,xdim,ydim,zdim);
        }
        ind = gen_ind(i,j,k,xdim+1,ydim,zdim);
        M.faces.push_back(face(owner,neighbour,point_inds));
        M.cells[owner].owns.push_back(ind);
        if(neighbour!=-1) M.cells[neighbour].neighbours.push_back(ind);
      }}}

  for(int i = 0; i < xdim; i++) {
    for(int j = 0; j < ydim+1; j++) {
      for(int k = 0; k < zdim; k++) {
        int owner,neighbour,ind;
        std::vector<int> point_inds {gen_ind(i,j,k,xdim+1,ydim+1,zdim+1),gen_ind(i,j,k+1,xdim+1,ydim+1,zdim+1),gen_ind(i+1,j,k+1,xdim+1,ydim+1,zdim+1),gen_ind(i+1,j,k,xdim+1,ydim+1,zdim+1)};
        if(j==0) {          
          owner=gen_ind(i,j,k,xdim,ydim,zdim);
          neighbour=-1;
        } else if(j==ydim){
          std::reverse(point_inds.begin(),point_inds.end());
          owner=gen_ind(i,j-1,k,xdim,ydim,zdim);
          neighbour=-1;
        } else {
          std::reverse(point_inds.begin(),point_inds.end());
          owner=gen_ind(i,j-1,k,xdim,ydim,zdim);
          neighbour=gen_ind(i,j,k,xdim,ydim,zdim);
        }
        ind = (xdim+1)*ydim*zdim + gen_ind(i,j,k,xdim,ydim+1,zdim);
        M.faces.push_back(face(owner,neighbour,point_inds));
        M.cells[owner].owns.push_back(ind);
        if(neighbour!=-1) M.cells[neighbour].neighbours.push_back(ind);
      }}}
  
  for(int i = 0; i < xdim; i++) {
    for(int j = 0; j < ydim; j++) {
      for(int k = 0; k < zdim+1; k++) {
        int owner,neighbour,ind;
        std::vector<int> point_inds {gen_ind(i,j,k,xdim+1,ydim+1,zdim+1),gen_ind(i+1,j,k,xdim+1,ydim+1,zdim+1),gen_ind(i+1,j+1,k,xdim+1,ydim+1,zdim+1),gen_ind(i,j+1,k,xdim+1,ydim+1,zdim+1)};
        if(k==0) {
          owner=gen_ind(i,j,k,xdim,ydim,zdim);
          neighbour=-1;
        } else if(k==zdim){
          std::reverse(point_inds.begin(),point_inds.end());
          owner=gen_ind(i,j,k-1,xdim,ydim,zdim);
          neighbour=-1;
        } else {
          std::reverse(point_inds.begin(),point_inds.end());
          owner=gen_ind(i,j,k-1,xdim,ydim,zdim);
          neighbour=gen_ind(i,j,k,xdim,ydim,zdim);
        }
        ind = (xdim+1)*ydim*zdim + xdim*(ydim+1)*zdim+gen_ind(i,j,k,xdim,ydim,zdim+1);
        M.faces.push_back(face(owner,neighbour,point_inds));
        M.cells[owner].owns.push_back(ind);
        if(neighbour!=-1) M.cells[neighbour].neighbours.push_back(ind);
      }}}

  //face order is all (nx+1)*ydim*zdim xfaces then nx*(ydim+1)*zdim yfaces, then nx*ydim*(zdim+1) zfaces
  std::vector<int> xu_face_inds,xd_face_inds;
  for(int j = 0; j < ydim; j++) {
    for(int k = 0; k < zdim; k++) {
      xu_face_inds.push_back(gen_ind(0,j,k,xdim+1,ydim,zdim));
      xd_face_inds.push_back(gen_ind(xdim,j,k,xdim+1,ydim,zdim));
    }}
  M.patches.push_back(patch(xu_face_inds,"xu",PT::PATCH));
  M.patches.push_back(patch(xd_face_inds,"xd",PT::PATCH));
  
  std::vector<int> yu_face_inds,yd_face_inds;
  for(int i = 0; i < xdim; i++) {
    for(int k = 0; k < zdim; k++) {
      yu_face_inds.push_back((xdim+1)*ydim*zdim+gen_ind(i,0,k,xdim,ydim+1,zdim));
      yd_face_inds.push_back((xdim+1)*ydim*zdim+gen_ind(i,ydim,k,xdim,ydim+1,zdim));
    }}
  M.patches.push_back(patch(yu_face_inds,"yu",PT::PATCH));
  M.patches.push_back(patch(yd_face_inds,"yd",PT::PATCH));

  std::vector<int> zu_face_inds,zd_face_inds;
  for(int i = 0; i < xdim; i++) {
    for(int j = 0; j < ydim; j++) {
      zu_face_inds.push_back((xdim+1)*ydim*zdim+xdim*(ydim+1)*zdim+gen_ind(i,j,0,xdim,ydim,zdim+1));
      zd_face_inds.push_back((xdim+1)*ydim*zdim+xdim*(ydim+1)*zdim+gen_ind(i,j,zdim,xdim,ydim,zdim+1));
    }}
  M.patches.push_back(patch(zu_face_inds,"zu",PT::PATCH));
  M.patches.push_back(patch(zd_face_inds,"zd",PT::PATCH));

  return M;
}

Mesh Mesh::make_2D_cartesian_mesh(int xdim, int ydim) {
  Mesh M = Mesh::make_3D_cartesian_mesh(xdim,ydim,1);
  M.patches[4].patch_type=PT::EMPTY;
  M.patches[4].face_inds.insert(M.patches[4].face_inds.end(),M.patches[5].face_inds.begin(),M.patches[5].face_inds.end());
  M.patches.pop_back();
  M.patches[4].patch_name="frontAndBackPlanes";
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
    points_out << boost::format("(%f %f %f)\n") % this->points[ind].x % this->points[ind].y % this->points[ind].z;
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
    internalfaces -= A.face_inds.size();
  }
  
  for(int i = 0; i < 17; i++) {
    getline(owner_template,line);
    if(i==12) {
      owner_out << boost::format("\t\tnote        \"nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d\";\n") % (int)this->points.size() % this->ncells() % (int)this->faces.size() % internalfaces;
    } else {
      owner_out << line << "\n";
    }
    getline(neighbour_template,line);
    if(i==12) {
      neighbour_out << boost::format("\t\tnote        \"nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d\";\n") % (int)this->points.size() % this->ncells() % (int)this->faces.size() % internalfaces;
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

    switch(P.patch_type) {
    case PT::PATCH:
      boundary_out << boost::format("\t%s\n\t{\n\t\ttype patch;\n\t\tnFaces %d;\n\t\tstartFace %d;\n\t}\n") %P.patch_name % P.face_inds.size() % start_face;
      break;
    case PT::EMPTY:
      boundary_out << boost::format("\t%s\n\t{\n\t\ttype empty;\n\t\tinGroups 1(empty);\n\t\tnFaces %d;\n\t\tstartFace %d;\n\t}\n") % P.patch_name % P.face_inds.size() % start_face;
      break;
    case PT::WALL:
      boundary_out << boost::format("\t%s\n\t{\n\t\ttype wall;\n\t\tinGroups 1(wall);\n\t\tnFaces %d;\n\t\tstartFace %d;\n\t}\n") % P.patch_name % P.face_inds.size() % start_face;
      break;
    default:
      printf("Other patch types not implemented\n");
      break;
    }

    start_face+=P.face_inds.size();
    patch_count++;
    for(int ind : P.face_inds) { //make so works with arbitrary numbers of points
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

  //inital condition
  std::ifstream U_template("initial_cond_templates/U");
  std::ifstream p_template("initial_cond_templates/p");
  std::ifstream k_template("initial_cond_templates/k");
  std::ifstream nut_template("initial_cond_templates/nut");
  std::ifstream omega_template("initial_cond_templates/omega");

  std::ofstream U_out("0/U");
  std::ofstream p_out("0/p");
  std::ofstream k_out("0/k");
  std::ofstream nut_out("0/nut");
  std::ofstream omega_out("0/omega");
  
  while(getline(U_template,line)) {
    U_out << line << '\n';
  }
  while(getline(p_template,line)) {
    p_out << line << '\n';
  }
  while(getline(k_template,line)) {
    k_out << line << '\n';
  }
  while(getline(nut_template,line)) {
    nut_out << line << '\n';
  }
  while(getline(omega_template,line)) {
    omega_out << line << '\n';
  }

  for(patch p : this->patches) {
    std::string bc;
    if(p.patch_type==PT::EMPTY) {
      bc="type empty;";
    } else {
      bc="type zeroGradient;";
    }
    U_out << p.patch_name << " { " << bc << " }\n\n";
    p_out << p.patch_name << " { " << bc << " }\n\n";
    k_out << p.patch_name << " { " << bc << " }\n\n";
    nut_out << p.patch_name << " { " << bc << " }\n\n";
    omega_out << p.patch_name << " { " << bc << " }\n\n";
  }
  
  U_out << "}";
  p_out << "}";
  k_out << "}";
  nut_out << "}";
  omega_out << "}";
  
  U_template.close();
  p_template.close();
  k_template.close();
  nut_template.close();
  omega_template.close();

  U_out.close();
  p_out.close();
  k_out.close();
  nut_out.close();
  omega_out.close();
}
