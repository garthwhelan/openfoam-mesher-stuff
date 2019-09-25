#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<boost/format.hpp>
#include<algorithm>
#include "point.hpp"
#include "mesh.hpp"
#include<random>
//#include<stdexcept>

//use more assert

int patch::patch_count = 0;

const double relax_mesh_const = 0.0005;//replace with 0.1*cell size

void Mesh::relax_mesh_2D(point (*v_dir)(point)) {
  std::random_device r;
  std::default_random_engine e1(r());
  std::uniform_real_distribution<double> uniform_dist(0,2*M_PI);

  for(int j = 0; j < 10; j++) {
    for(int i = 0; i < this->points.size(); i++) {
      if(this->points[i].z==0) {
        if(not pointind_in_boundary(i)) {
          for(int k = 0; k < 10; k++) {
            
            std::vector<int> face_inds = this->get_faces_by_point(i);

            point orig_point = this->points[i];
            double orig_skewness = this->skewness(i,face_inds);
            double orig_alignment = this->flow_alignment(i,face_inds,v_dir);
            double orig_badness = orig_skewness+orig_alignment;
          
            double rand_theta = uniform_dist(e1);
            this->points[i]+=point{relax_mesh_const*cos(rand_theta),relax_mesh_const*sin(rand_theta),0};
            double new_badness = this->skewness(i,face_inds)+this->flow_alignment(i,face_inds,v_dir);
            if(new_badness>orig_badness) {
              this->points[i]=orig_point;
            } else {
              bool found = false;
              for(point &p : this->points) {
                if(p==orig_point+point{0,0,1}) {p = this->points[i]+point{0,0,1}; found=true; break;}
              }
              assert(found);
            }
          }
        }
      }
    }
  }
}

int gen_ind(int i, int j, int k, int nx,int ny, int nz) {
  return i*(ny*nz)+j*(nz)+k;
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

//for combining meshes or patches
void Mesh::remove_duplicate_points() {  
  std::vector<point> point_array;
  for(int i = 0; i < this->points.size(); i++) {
    if(not this->is_pointind_used(i)) continue;
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

void Mesh::correctly_orient_face(int f_ind) {//todo refactor
  std::vector<int>& points = this->faces[f_ind].point_inds;
  
  point face_com = this->face_CoM(f_ind);
  point normal_vec = this->face_normal(f_ind);
  assert(point::len(normal_vec)==1);
  
  point p = this->points[points[0]]-face_com;
  p-=normal_vec*point::dot(normal_vec,p);
  p/=point::len(p);
  point second_normal = point::cross(normal_vec,p);
  second_normal/=point::len(second_normal);
  
  std::sort(points.begin()+1,points.end(),
            [&](int ind1, int ind2){
              point p1 = this->points[ind1]-face_com;
              point p2 = this->points[ind2]-face_com;
              p1-=normal_vec*point::dot(normal_vec,p1);
              p2-=normal_vec*point::dot(normal_vec,p2);
              assert(point::dot(p1,normal_vec)==0);
              assert(point::dot(p2,normal_vec)==0);
              p1/=point::len(p1);
              p2/=point::len(p2);              
              float adjacent = point::dot(p1,p);

              float opposite = point::dot(p1,second_normal);
              float theta1 = atan2(opposite,adjacent);
              if(theta1<0) theta1+=2*M_PI;
              adjacent = point::dot(p2,p);
              opposite = point::dot(p2,second_normal);
              float theta2 = atan2(opposite,adjacent);
              if(theta2<0) theta2+=2*M_PI;
              return theta1 < theta2;
            }
            );
  for(int i = 0; i < points.size(); i++) {
    point p1 = this->points[points[i]]-face_com;
              p1-=normal_vec*point::dot(normal_vec,p1);
              p1/=point::len(p1);
              float adjacent = point::dot(p1,p);
              float opposite = point::dot(p1,second_normal);
              float theta = atan2(opposite,adjacent);
              if(theta<0) theta+=2*M_PI;
              std::cout<<theta<<" ";
  }
  std::cout << "\n";
  if(point::dot(this->cell_CoM(this->faces[f_ind].owner)-this->face_CoM(f_ind),this->face_normal(f_ind))>0) {
   std::reverse(this->faces[f_ind].point_inds.begin(),this->faces[f_ind].point_inds.end());
  }
}


void Mesh::correctly_orient_faces() {
  for(int i = 0; i < this->faces.size(); i++) {
    this->correctly_orient_face(i);
  }
}

void Mesh::remove_invalid_faces() {
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
}

void Mesh::remove_invalid_cells() {
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
}

void Mesh::cleanup_boundary() {
  //put faces outside of patches into a new patch
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
      }}}}


void Mesh::cleanup() {

  this->remove_duplicate_points();
  
  this->remove_invalid_faces();
  
  //remove empty patches
  this->patches.erase(std::remove_if(this->patches.begin(),this->patches.end(),
      [](patch p){return p.face_inds.size()==0;}),this->patches.end());
  
  this->remove_invalid_cells();  
  this->correctly_orient_faces();
  this->cleanup_boundary(); 
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
          owner=gen_ind(i-1,j,k,xdim,ydim,zdim);
          neighbour=-1;
        } else {
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
          owner=gen_ind(i,j-1,k,xdim,ydim,zdim);
          neighbour=-1;
        } else {
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
          owner=gen_ind(i,j,k-1,xdim,ydim,zdim);
          neighbour=-1;
        } else {
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

  M.correctly_orient_faces();
  
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


void remove_duplicates_in_array(std::vector<int> &arr) {
  std::sort(arr.begin(),arr.end());
  arr.erase(std::unique(arr.begin(),arr.end()),arr.end());
}

Mesh Mesh::make_wedge_mesh(int xdim, int ydim) {
  Mesh M = Mesh::make_3D_cartesian_mesh(xdim,ydim,1);
  M.patches[4].patch_type=PT::WEDGE;
  M.patches[4].patch_name="backwedge";
  M.patches[5].patch_type=PT::WEDGE;
  M.patches[5].patch_name="frontwedge";
  for(point &p : M.points) {
    p.z-=0.5;
    p.z*=0.02*(p.x);
  }
  M.remove_duplicate_points();
  for(int i = 0; i < M.faces.size(); i++) {
    remove_duplicates_in_array(M.faces[i].point_inds);
    if(M.faces[i].point_inds.size()==2) {
      M.faces[i].owner=-1;
      M.faces[i].neighbour=-1;
    }
  }
  M.remove_invalid_faces();
  M.cleanup();
  return M;
}

void Mesh::write_mesh() {
  std::string line;
  
  std::ifstream points_template("outputfile_templates/points_template");
  std::ifstream faces_template("outputfile_templates/faces_template");
  std::ifstream owner_template("outputfile_templates/owner_template");
  std::ifstream neighbour_template("outputfile_templates/neighbour_template");
  std::ifstream boundary_template("outputfile_templates/boundary_template");
  std::ifstream *mesh_template_files[5] = {&points_template,&faces_template,&owner_template,&neighbour_template,&boundary_template};

  std::ofstream points_out("polyMesh/points");
  std::ofstream faces_out("polyMesh/faces");
  std::ofstream owner_out("polyMesh/owner");
  std::ofstream neighbour_out("polyMesh/neighbour");
  std::ofstream boundary_out("polyMesh/boundary");
  std::ofstream *mesh_out_files[5] = {&points_out,&faces_out,&owner_out,&neighbour_out,&boundary_out};

  while(getline(points_template,line)) {
    points_out << line << '\n';
  }
  points_out << boost::format("%d\n(\n") % (int)this->points.size();
  for(int ind = 0; ind < this->points.size(); ind++) {
    points_out << boost::format("(%f %f %f)\n") % this->points[ind].x % this->points[ind].y % this->points[ind].z;
  }
  
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
      switch(this->faces[ind].point_inds.size()) {
      case 3:
        faces_out << boost::format("3(%d %d %d)\n") % A.point_inds[0] % A.point_inds[1] % A.point_inds[2];
        break;
      case 4:
        faces_out << boost::format("4(%d %d %d %d)\n") % A.point_inds[0] % A.point_inds[1] % A.point_inds[2] % A.point_inds[3];
        break;
      default:
        assert(false);
      }
      owner_out << boost::format("%d\n") % A.owner;
      neighbour_out << boost::format("%d\n") % A.neighbour;
      start_face++;
    }
  }

  boundary_out << boost::format("%d\n(\n") % this->patches.size();
  int patch_count = 1;
  int wedges_out = 0;
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
    case PT::WEDGE:
      assert(wedges_out!=2);
      boundary_out << boost::format("\t%s\n\t{\n\t\ttype wedge;\n\t\tinGroups 1(wedge);\n\t\tnFaces %d;\n\t\tstartFace %d;\n\t}\n") % P.patch_name % P.face_inds.size() % start_face;
      wedges_out++;//need two?
      break;
    default:
      assert(false);
      break;
    }

    start_face+=P.face_inds.size();
    patch_count++;
    for(int ind : P.face_inds) { //make so works with arbitrary numbers of points
      switch(this->faces[ind].point_inds.size()) {
      case 3:
        faces_out << boost::format("3(%d %d %d)\n") % this->faces[ind].point_inds[0] % this->faces[ind].point_inds[1] % this->faces[ind].point_inds[2];
        break;
      case 4:
        faces_out << boost::format("4(%d %d %d %d)\n") % this->faces[ind].point_inds[0] % this->faces[ind].point_inds[1] % this->faces[ind].point_inds[2] % this->faces[ind].point_inds[3];
        break;
      default:
        assert(false);
      }
      owner_out << boost::format("%d\n") % this->faces[ind].owner;
    }
  }

  for(auto e : mesh_template_files) {
    e->close();
  }
  
  for(auto e : mesh_out_files) {
    *e << ")";
    e->close();
  }
  
  //inital condition
  std::ifstream U_template("initial_cond_templates/U");
  std::ifstream p_template("initial_cond_templates/p");
  std::ifstream k_template("initial_cond_templates/k");
  std::ifstream nut_template("initial_cond_templates/nut");
  std::ifstream omega_template("initial_cond_templates/omega");
  std::ifstream *template_files[5] = {&U_template,&p_template,&k_template,&nut_template,&omega_template};

  std::ofstream U_out("0/U");
  std::ofstream p_out("0/p");
  std::ofstream k_out("0/k");
  std::ofstream nut_out("0/nut");
  std::ofstream omega_out("0/omega");
  std::ofstream *out_files[5] = {&U_out,&p_out,&k_out,&nut_out,&omega_out};  

  for(int i = 0; i < 5; i++) {
    while(getline(*template_files[i],line)) {
      *out_files[i] << line << '\n';
    }
  }

  for(patch p : this->patches) {
    std::string bc;
    if(p.patch_type==PT::EMPTY) {
      bc="type empty;";
    } else {
      bc="type zeroGradient;";
    }
    for(auto e : out_files) {
      *e << p.patch_name << " { " << bc << " }\n\n";
    }
  }

  for(auto e : out_files) {
    *e << "}";
    e->close();
  }
  for(auto e : template_files) {
    e->close();
  }
}

