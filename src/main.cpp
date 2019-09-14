#include<iostream>
#include<cmath>
#include "mesh.hpp"
#include "NACA4.hpp"
#include<chrono>
#include<algorithm>

#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now()
#define START_TIMER  start = std::chrono::high_resolution_clock::now()
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now()-start ).count() << " ms " << std::endl; 

int main() {
  INIT_TIMER;
  START_TIMER;
  Mesh M1 = Mesh::make_2D_cartesian_mesh(25,25);
  Mesh M2 = Mesh::make_2D_cartesian_mesh(25,25);

  for(point &p : M2.points) {
    p.y+=25;
  }

  Mesh M = Mesh::combine_meshes(M1,M2);
    
  // STOP_TIMER("here\n");
  std::cout << "owns:\n";
  for(int ind :  M.cells[0].owns) {
    std::cout << "[";
    for(int ind2 : M.faces[ind].point_inds) {
      M.points[ind2].print();
      std::cout << ",";
    } std::cout << "]";std::cout << "\n";
  }
  std::cout << "neighbours:\n";
  for(int ind :  M.cells[0].neighbours) {
    std::cout << "[";
    for(int ind2 : M.faces[ind].point_inds) {
      M.points[ind2].print();
      std::cout << ",";
    } std::cout << "]";std::cout << "\n";
  }

  M.write_mesh();
  // //Don't forget to call 'renumberMesh -overwrite' !!!!
  // std::cout<<"Done!\n";
  return 0;
}
