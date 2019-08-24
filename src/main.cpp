#include<iostream>
#include<cmath>
#include "mesh.hpp"
#include<chrono>

/* Example making a C mesh around a NACA4 (2408) airfoil
http://airfoiltools.com/airfoil/details?airfoil=naca2408-il
 */
const float T = 0.08;
const float M = 0.02;
const float P = 0.4;

#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now()
#define START_TIMER  start = std::chrono::high_resolution_clock::now()
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now()-start ).count() << " ms " << std::endl; 

float y_t(float x) {
  float a0=0.2969;
  float a1=-0.126;
  float a2=-0.3516;
  float a3=0.2843;
  float a4=-0.1036;
  return (T/0.2)*(a0*pow(x,0.5)+a1*x+a2*pow(x,2)+a3*pow(x,3)+a4*pow(x,4));
}

float theta(float x) {
	if (x<P) return atan(2.0*M/pow(P,2)*(P-x));
  return atan(2.0*M/pow(1.0-P,2)*(P-x));
}

//chord coordinate
float y_c(float x) {
	if (x<P) return M/(pow(P,2))*(2.0*P*x-pow(x,2));
  return M/(pow(1.0-P,2))*(1.0-2.0*P+2.0*P*x-pow(x,2));
}

//x and y coordinates as functions of
//x (chord x) from 0 to 1
float x_u(float x) {
	return x-y_t(x)*sin(theta(x));
}
float x_l(float x) {
	return x+y_t(x)*sin(theta(x));
}
float y_u(float x) {
	return y_c(x)+y_t(x)*cos(theta(x));
}
float y_l(float x) {
	return y_c(x)-y_t(x)*cos(theta(x));
}

int main() {
  int nx = 70;
  int ny = 50;

  int ntail = 20;
  
  Mesh M = make_cartesian_mesh(nx,ny,1);
  float x,y;
  float theta,r;
  for(point &p : M.points) {

    int i = p.p[0];
    int j = p.p[1];
    int k = p.p[2];
    if(i < ntail) {//trailing lower half mesh
      x=1.05+0.075*(ntail-1-i);
      x+=0.5*pow(x-1,2);
      y=-0.005*j-0.0005*j*j;
    } else if (i < nx-ntail+1) {//c mesh
      r=0.005*j+0.0005*j*j;
      theta=M_PI*(i-ntail)/(nx-2*ntail);
      x=-r*sin(theta)+1-sin(theta);
      y=-r*cos(theta);
      if(theta>(0.5*M_PI)) {
        y+=y_u(1-sin(theta));
      } else {
        y+=y_l(1-sin(theta));
      }

    } else {//trailing upper mesh
      x=1.05+0.075*(i-(nx-ntail+1));
      x+=0.5*pow(x-1,2);
      y=0.005*j+0.0005*j*j;
    }
    p.p[0]=x;
    p.p[1]=p.p[2];
    p.p[2]=y;
  }
  M.patches[4].pt=EMPTY;
  M.patches[4].faces.insert(M.patches[4].faces.end(),M.patches[5].faces.begin(),M.patches[5].faces.end());  
  M.patches.pop_back();
  M.patches[4].name="frontAndBackPlanes";
  M.patches[0].faces.insert(M.patches[0].faces.end(),M.patches[1].faces.begin(),M.patches[1].faces.end());
  M.patches.erase(M.patches.begin()+1);
  M.patches[0].name="back";
  M.patches[1].name="wall";
  M.patches[2].name="front";
  //M.patches[1].pt=WALL;//TODO
  std::cout << M.points.size() << ' ' << M.faces.size() << '\n';
  M.remove_duplicate_points();
  INIT_TIMER;
  START_TIMER;
  M.remove_duplicate_faces();
  STOP_TIMER("ABC");
  std::cout << M.points.size() << ' ' << M.faces.size() << '\n';
  M.write_mesh();
  //Don't forget to call 'renumberMesh -overwrite' !!!!
  std::cout<<"Done!\n";
  return 0;
}
