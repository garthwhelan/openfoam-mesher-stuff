#include<iostream>
#include "mesh.hpp"
#define BOOST_TEST_MODULE MyTest
#include<boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(point_test) {

  point p1 = point(1.0,2.0,3.0);
  point p2 = p1;
  point p3 = point(2.0,1.0,3.0);
  point p4 = point(3.0,3.0,6.0);
  
  BOOST_CHECK(p1==p2);
  BOOST_CHECK(!(p1==p3));
  BOOST_CHECK((p1+p3)==p4);
  p1+=p2;
  BOOST_CHECK(!(p1==p2));
  p1/=2.0;
  BOOST_CHECK(p1==p2);
}

BOOST_AUTO_TEST_CASE(face_test) {
  face f1,f2,f3;
  f1 = face(1, -1, {0,1,2,3});
  f2 = face(2, -1, {3,1,2,0});
  f3 = face(5,  4, {1,2,3,4});

  BOOST_CHECK(f1==f2);
  BOOST_CHECK(!(f1==f3));
}

BOOST_AUTO_TEST_CASE(mesh_quality_tools_test) {
  int nx=5;
  int ny=5;
  Mesh M = Mesh::make_2D_cartesian_mesh(nx,ny);
  for(int i = 0; i < M.points.size(); i++) {
    std::vector<int> face_inds = M.get_faces_by_point(i);
    BOOST_CHECK(M.skewness(i,face_inds)==0);
  }
  for(point &p : M.points) {
    p.y+=p.x;
  }
  for(int i = 0; i < M.faces.size(); i++) {
    if(i < (nx+1)*ny+nx*(ny+1)) {
      BOOST_CHECK(M.face_skewness(i)==0);
    } else {
      BOOST_CHECK(point::near(pow(M.face_skewness(i),2),0.5));
    }
  }
}

void Test_Mesh(Mesh M, bool cartesian = true) {

  std::vector<int> all_f_inds;
  for(patch p : M.patches) {
    for(int f_ind : p.face_inds) {
      all_f_inds.push_back(f_ind);
      BOOST_CHECK(M.faces[f_ind].owner>=0);
      BOOST_CHECK(M.faces[f_ind].neighbour==-1);
      BOOST_CHECK((f_ind<M.faces.size())&&(f_ind>=0));
    }
  }
  std::sort(all_f_inds.begin(),all_f_inds.end());
  int patch_length = all_f_inds.size();
  std::unique(all_f_inds.begin(),all_f_inds.end());
  BOOST_CHECK(all_f_inds.size()==patch_length);//no duplicate faces in patches

  for(int i = 0; i < M.faces.size(); i++) {
    bool in_arr = std::find(all_f_inds.begin(),all_f_inds.end(),i)<all_f_inds.end();
    if(!in_arr) {
      BOOST_CHECK(M.faces[i].owner>=0);
      BOOST_CHECK(M.faces[i].neighbour>=0);
      BOOST_CHECK(M.faces[i].owner!=M.faces[i].neighbour);
    }
  }

  if(cartesian) {
    for(face f : M.faces) BOOST_CHECK(f.point_inds.size()==4);
    for(cell c : M.cells) BOOST_CHECK((c.owns.size()+c.neighbours.size())==6);
  }
  
  //checking consistency between face owners/neighbours and cells
  for(int i = 0; i < M.cells.size(); i++) {
    for(int f_ind : M.cells[i].owns) {
      BOOST_CHECK_MESSAGE((M.faces[f_ind].owner)==i, "cell_ind: " << i << " f->owner: " << M.faces[f_ind].owner << "\n");
    }
    for(int f_ind : M.cells[i].neighbours) {
      BOOST_CHECK_MESSAGE((M.faces[f_ind].neighbour)==i, "cell_ind: " << i << " f->neighbour: " << M.faces[f_ind].neighbour << "\n");
    }
  }
}

BOOST_AUTO_TEST_CASE(mesh_test) {

  Mesh M = Mesh::make_3D_cartesian_mesh(1,1,1);
  BOOST_CHECK(M.cell_CoM(0)==point(0.5,0.5,0.5));
  
  int nx = 7;
  int ny = 9;
  int nz = 11;
  M = Mesh::make_3D_cartesian_mesh(nx,ny,nz);  
  BOOST_CHECK(M.points.size()==(nx+1)*(ny+1)*(nz+1));
  BOOST_CHECK(M.faces.size()==(nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1));
  Test_Mesh(M);

  M = Mesh::make_2D_cartesian_mesh(nx,ny);
  BOOST_CHECK(M.points.size()==(nx+1)*(ny+1)*2);
  BOOST_CHECK(M.faces.size()==(nx+1)*ny+nx*(ny+1)+nx*ny*2);
  Test_Mesh(M);
  
  Mesh M2 = Mesh::make_3D_cartesian_mesh(nx,ny,nz);
  M2.remove_duplicate_points();//shouldn't do anything
  M2.remove_duplicate_faces();//shouldn't do anything either
  BOOST_CHECK(M2.points.size()==(nx+1)*(ny+1)*(nz+1));
  BOOST_CHECK(M2.faces.size()==(nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1));
  Test_Mesh(M2);
  
  //BOOST_CHECK(M.points==M2.points); properly compare std::vector
  // //BOOST_CHECK(M.faces==M2.face_inds);
  // //BOOST_CHECK(M.patches==M2.patches);

  nx = 6;
  ny = 7;
  Mesh M3 = Mesh::make_2D_cartesian_mesh(nx,ny);
  float r,theta;
  for(point &p : M3.points) {
    int i = p.x;
    int j = p.y;
    int k = p.z;
    r=j+1.0;
    theta=2.0*M_PI*((float)i)/((float)nx);//so theta=0 and theta=2*PI are both points
    p.x=r*sin(theta);
    p.y=r*cos(theta);
  }
  Test_Mesh(M3,false);
  M3.remove_duplicate_points();
  std::cout<< M3.points.size() << " " << M3.faces.size() << "\n";
  BOOST_CHECK(M3.points.size()==(nx)*(ny+1)*(2));
  Test_Mesh(M3,false);
  M3.remove_duplicate_faces();
  BOOST_CHECK(M3.faces.size()==(nx+1)*ny+nx*(ny+1)+nx*ny*2-ny);
  Test_Mesh(M3,false);

  M = Mesh::make_2D_cartesian_mesh(5,5);
  M.patches.push_back(patch({},"removed",PT::PATCH));
  std::vector<int> A {7,11,12,13,17};
  M.remove_cells(A,M.patches[5]);
  M.write_mesh();
}
