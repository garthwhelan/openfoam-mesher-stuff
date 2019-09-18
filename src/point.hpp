#pragma once

const float small_number = 1e-4;

struct point {
  float x,y,z;
  bool operator==(const point&) const;
  point operator+(const point&) const;
  point operator-(const point&) const;
  void operator+=(const point&);
  void operator/=(const float&);
  static float dot(point,point);
  static point cross(point,point);
  static bool near(float,float);
  static float len(point);
  point(float x, float y, float z) : x{x},y{y},z{z} {}
  void print();
};
