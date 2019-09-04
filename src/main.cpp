#include<iostream>
#include<cmath>
#include "mesh.hpp"
#include "NACA4.hpp"
#include<chrono>
#include<algorithm>

#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now()
#define START_TIMER  start = std::chrono::high_resolution_clock::now()
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now()-start ).count() << " ms " << std::endl; 


//best to keep increase in size per adjacent cells 20%
int main() {
  Mesh M1 = make_cartesian_mesh(10,1,10);
  Mesh M2 = make_cartesian_mesh(10,1,10);
  std::transform(M2.points.begin(),M2.points.end(),M2.points.begin(),[](point &p){p.p[0]+=10;return p;});
  Mesh M = Mesh::combine_meshes(M1,M2);
  M.write_mesh();
  //Don't forget to call 'renumberMesh -overwrite' !!!!
  std::cout<<"Done!\n";
  return 0;
}
