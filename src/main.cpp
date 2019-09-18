#include<iostream>
#include<cmath>
//#include "point.hpp"
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
  Mesh M = test();//Mesh_w_flap();
  M.cleanup();
  M.write_mesh();
  // //Don't forget to call 'renumberMesh -overwrite' !!!!
  std::cout<<"Done!\n";
  return 0;
}
