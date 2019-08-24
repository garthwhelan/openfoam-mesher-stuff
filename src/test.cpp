#include<iostream>
#include "mesh.hpp"
#define BOOST_TEST_MODULE MyTest
#include<boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(point_test) {
  //my_class test_object("abc");
  
  //BOOST_CHECK(test_object.is_valid());
  point p1 = point(1.0,2.0,3.0);
  point p2 = p1;
  point p3 = point(2.0,1.0,3.0);
  
  BOOST_CHECK(p1==p2);
  BOOST_CHECK(!(p1==p3));
}

BOOST_AUTO_TEST_CASE(face_test) {
  face f1,f2,f3;
  f1 = face(1, -1, {0,1,2,3});
  f2 = face(2, -1, {3,1,2,0});
  f3 = face(5,  4, {1,2,3,4});

  BOOST_CHECK(f1==f2);
  BOOST_CHECK(!(f1==f3));
}

void Test_Mesh(Mesh M) {

  std::vector<int> all_f_inds;
  for(patch p : M.patches) {
    for(int f_ind : p.faces) {
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
}

BOOST_AUTO_TEST_CASE(mesh_test) {
  int nx = 7;
  int ny = 9;
  int nz = 11;
  Mesh M = make_cartesian_mesh(nx,ny,nz);
  BOOST_CHECK(M.points.size()==(nx+1)*(ny+1)*(nz+1));
  BOOST_CHECK(M.faces.size()==(nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1));
  Test_Mesh(M);
  
  Mesh M2 = make_cartesian_mesh(nx,ny,nz);
  M2.remove_duplicate_points();//shouldn't do anything
  M2.remove_duplicate_faces();//see above
  BOOST_CHECK(M2.points.size()==(nx+1)*(ny+1)*(nz+1));
  BOOST_CHECK(M2.faces.size()==(nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1));
  Test_Mesh(M2);
  
  //BOOST_CHECK(M.points==M2.points); properly compare std::vector
  //BOOST_CHECK(M.faces==M2.faces);
  //BOOST_CHECK(M.patches==M2.patches);

  nx = 6;
  ny = 7;
  Mesh M3 = make_cartesian_mesh(nx,ny,1);
  float r,theta;
  for(point &p : M3.points) {
    int i = p.p[0];
    int j = p.p[1];
    int k = p.p[2];
    r=j+1.0;
    theta=2.0*M_PI*((float)i)/((float)nx);//so theta=0 and theta=2*PI are both points
    p.p[0]=r*sin(theta);
    p.p[1]=r*cos(theta);
  }
  Test_Mesh(M3);
  M3.remove_duplicate_points();
  BOOST_CHECK(M3.points.size()==(nx)*(ny+1)*(2));
  Test_Mesh(M3);
  M3.remove_duplicate_faces();  
  //BOOST_CHECK(M3.faces.size()==(nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1)-ny);
  Test_Mesh(M3);
  
}
