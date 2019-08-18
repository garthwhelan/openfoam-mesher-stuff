#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include <boost/format.hpp>
#include <algorithm>
//TODO: other patch types, parallelization

int gen_ind(int i, int j, int k, int nx,int ny, int nz) {
  return i*(ny*nz)+j*(nz)+k;
}

struct point {
  float p[3];
};

point make_point(float x, float y, float z) {
  point po;
  po.p[0]=x;
  po.p[1]=y;
  po.p[2]=z;
  return po;
}

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

bool face_order(face f1, face f2) {
  return f1.owner>f2.owner;
}

//TODO: refactor
void order_mesh(Mesh M) {
  //for uppertriangularity faces should be in order of owner cell IDs
  //remove faces associated with patches temporarily
  std::vector<face> internal_faces;
  for(face f : M.faces) {
    if(f.neighbour!=-1) internal_faces.push_back(f);
  }
  std::sort(internal_faces.begin(),internal_faces.end(),face_order);
          
  std::vector<std::vector<face>> patch_face_lists;
  for(patch p : M.patches) {
    std::cout << "patch:" << p.name << "\n";
    std::vector<face> face_list;
    for(int face_ind : p.faces) {
      face_list.push_back(M.faces[face_ind]);
    }
    std::sort(face_list.begin(),face_list.end(),face_order);
    patch_face_lists.push_back(face_list);
  }
  
  M.faces.clear();
  M.faces.insert(M.faces.begin(),internal_faces.begin(),internal_faces.end());

  for(int i = 0; i < M.patches.size(); i++) {
    std::vector<face> face_list = patch_face_lists[i];
    int initial_size = M.faces.size();
    M.faces.insert(M.faces.end(),face_list.begin(),face_list.end());
    int end_size = M.faces.size();
    M.patches[i].faces.clear();
    for(int j = initial_size; j < end_size; j++) M.patches[i].faces.push_back(j);
  }
}

Mesh make_cartesian_mesh(int xdim, int ydim, int zdim) {
  
  Mesh M;

  M.ncells=xdim*ydim*zdim;
  
  for(int i = 0; i < xdim+1; i++) {
    for(int j = 0; j < ydim+1; j++) {
      for(int k = 0; k < zdim+1; k++) {
        (M.points).push_back(make_point(i,j,k));
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
          std::reverse(A.point_inds.begin(),A.point_inds.end());
          A.owner=gen_ind(i,j,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else if(i==xdim){
          A.owner=gen_ind(i-1,j,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else {
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
          std::reverse(A.point_inds.begin(),A.point_inds.end());
          A.owner=gen_ind(i,j,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else if(j==ydim){
          A.owner=gen_ind(i,j-1,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else {
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
          std::reverse(A.point_inds.begin(),A.point_inds.end());
          A.owner=gen_ind(i,j,k,xdim,ydim,zdim);
          A.neighbour=-1;
        } else if(k==zdim){
          A.owner=gen_ind(i,j,k-1,xdim,ydim,zdim);
          A.neighbour=-1;
        } else {
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
  
  order_mesh(M);
  
  return M;
}

void write_mesh(Mesh M) {
  std::string line;
  //points
  std::ifstream points_template("points_template");
  std::ofstream points_out("polyMesh/points");
  while(getline(points_template,line)) {
    points_out << line << '\n';
  }
  points_out << boost::format("%d\n(\n") % (int)M.points.size();
  for(int ind = 0; ind < M.points.size(); ind++) {
    points_out << boost::format("(%f %f %f)\n") % M.points[ind].p[0] % M.points[ind].p[1] % M.points[ind].p[2];
  }
  points_out << ")\n";
  points_template.close();
  points_out.close();

  std::ifstream faces_template("faces_template");
  std::ifstream owner_template("owner_template");
  std::ifstream neighbour_template("neighbour_template");
  std::ifstream boundary_template("boundary_template");

  std::ofstream faces_out("polyMesh/faces");
  std::ofstream owner_out("polyMesh/owner");
  std::ofstream neighbour_out("polyMesh/neighbour");
  std::ofstream boundary_out("polyMesh/boundary");

  while(getline(faces_template,line)) {
    faces_out << line << '\n';
  }

  int internalfaces = M.faces.size();
  for(patch A : M.patches) {
    internalfaces -= A.faces.size();
  }
  
  for(int i = 0; i < 17; i++) {
    getline(owner_template,line);
    if(i==12) {
      owner_out << boost::format("\t\tnote        \"nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d\";\n") % (int)M.points.size() % M.ncells % (int)M.faces.size() % internalfaces;
    } else {
      owner_out << line << "\n";
    }
    getline(neighbour_template,line);
    if(i==12) {
      neighbour_out << boost::format("\t\tnote        \"nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d\";\n") % (int)M.points.size() % M.ncells % (int)M.faces.size() % internalfaces;
    } else {
      neighbour_out << line << "\n";
    }
  }
  
  faces_out << boost::format("%d\n(\n") % (int)M.faces.size();
  owner_out << boost::format("%d\n(\n") % (int)M.faces.size();
  neighbour_out << boost::format("%d\n(\n") % internalfaces;
    
  while(getline(boundary_template,line)) {
    boundary_out << line << '\n';
  }
    
  int start_face = 0;
    
  for(int ind = 0; ind < M.faces.size(); ind++) {
    face A = M.faces[ind];
    if(A.neighbour!=-1) {
      faces_out << boost::format("4(%d %d %d %d)\n") % A.point_inds[0] % A.point_inds[1] % A.point_inds[2] % A.point_inds[3];
      owner_out << boost::format("%d\n") % A.owner;
      neighbour_out << boost::format("%d\n") % A.neighbour;
      start_face++;
    }
  }

  boundary_out << boost::format("%d\n(\n") % M.patches.size();
  int patch_count = 1;
  for(patch P : M.patches) {

    switch(P.pt) {
    case PATCH:
      boundary_out << boost::format("\t%s\n\t{\n\t\ttype patch;\n\t\tnFaces %d;\n\t\tstartFace %d;\n\t}\n") %P.name % P.faces.size() % start_face;
      break;
    case EMPTY:
      boundary_out << boost::format("\t%s\n\t{\n\t\ttype empty;\n\t\tinGroups 1(empty);\n\t\tnFaces %d;\n\t\tstartFace %d;\n\t}\n") % P.name % P.faces.size() % start_face;
      break;
    default:
      printf("Other patch types not implemented\n");
      break;
    }

    start_face+=P.faces.size();
    patch_count++;
    for(int ind : P.faces) {
      faces_out << boost::format("4(%d %d %d %d)\n") % M.faces[ind].point_inds[0] % M.faces[ind].point_inds[1] % M.faces[ind].point_inds[2] % M.faces[ind].point_inds[3];
      owner_out << boost::format("%d\n") % M.faces[ind].owner;
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

int main() {
  Mesh M = make_cartesian_mesh(100,100,1);
  for(point p : M.points) {
    p.p[0]*=0.1;
    p.p[1]*=0.1;
    p.p[2]*=0.1;
  }
  M.patches[4].pt=EMPTY;
  //M.patches[5].pt=EMPTY;
  M.patches[4].faces.insert(M.patches[4].faces.end(),M.patches[5].faces.begin(),M.patches[5].faces.end());
  M.patches.pop_back();
  write_mesh(M);
  
  printf("here\n");
  return 0;
}
