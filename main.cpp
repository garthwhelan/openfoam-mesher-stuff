#include<iostream>
#include<cmath>
#include "mesh.hpp"

const float M = 0.02;
const float P = 0.4;

float y_t(float x) {
  float a0=0.2969;
  float a1=-0.126;
  float a2=-0.3516;
  float a3=0.2843;
  float a4=-0.1036;
  float T=0.1;
  return (T/0.2)*(a0*pow(x,0.5)+a1*x+a2*pow(x,2)+a3*pow(x,3)+a4*pow(x,4));
}

float theta(float x) {//check atan==arctan
	if (x<P) return atan(2.0*M/pow(P,2)*(P-x));
  return atan(2.0*M/pow(1.0-P,2)*(P-x));
}

float y_c(float x) {//chord coordinate
	if (x<P) return M/(pow(P,2))*(2.0*P*x-pow(x,2));
  return M/(pow(1.0-P,2))*(1.0-2.0*P+2.0*P*x-pow(x,2));
}

float x_u(float x) {// and y coors for upper and lower surfaces
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
  Mesh M = make_cartesian_mesh(100,100,1);
  // for(point &p : M.points) {
  //   int i = p.p[0];
  //   int j = p.p[1];
  //   int k = p.p[2];
  //   if(i < 30) {
  //     return 0;
  //   } else if (i < 45) {
  //     return 0;
  //   } else if (i < 60) {
  //     return 0;
  //   } else {
  //     return 0;
  //   }
  // }
  // M.patches[4].pt=EMPTY;
  // //M.patches[5].pt=EMPTY;
  // M.patches[4].faces.insert(M.patches[4].faces.end(),M.patches[5].faces.begin(),M.patches[5].faces.end());
  // M.patches.pop_back();
  write_mesh(M);

  std::cout<<"Done!\n";
  return 0;
}
