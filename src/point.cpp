#include<iostream>
#include "point.hpp"
#include<cmath>

float point::len(point a) {
  return pow(pow(a.x,2)+pow(a.y,2)+pow(a.z,2),0.5);
}

bool point::near(float a, float b) {
  return (std::abs(a-b)<small_number);
}

void point::print() {
  std::cout << "[" << this->x << "," << this->y << ","<< this->z << "]";
}

bool point::operator==(const point& rhs) const {
  return near(this->x,rhs.x)&&near(this->y,rhs.y)&&near(this->z,rhs.z);
}

point point::operator+(const point& rhs) const {
  point p = (*this);
  p.x+=rhs.x;
  p.y+=rhs.y;
  p.z+=rhs.z;
  return p;
}

point point::operator-(const point& rhs) const {
  point p = (*this);
  p.x-=rhs.x;
  p.y-=rhs.y;
  p.z-=rhs.z;
  return p;
}

void point::operator+=(const point& rhs) {
  this->x+=rhs.x;
  this->y+=rhs.y;
  this->z+=rhs.z;
}

void point::operator/=(const float& rhs) {
  this->x/=rhs;
  this->y/=rhs;
  this->z/=rhs;
}

point point::cross(point A, point B) {
  point res{0,0,0};
  res.x=A.y*B.z-A.z*B.y;
  res.y=A.z*B.x-A.x*B.z;
  res.z=A.x*B.y-A.y*B.x;
  return res;
}

float point::dot(point A, point B) {
  return A.x*B.x+A.y*B.y+A.z*B.z;
}
