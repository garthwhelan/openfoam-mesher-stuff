#pragma once

const double small_number = 1e-6;

struct point {
  double x,y,z;
  bool operator==(const point&) const;
  point operator+(const point&) const;
  point operator-(const point&) const;
  void operator+=(const point&);
  void operator-=(const point&);
  void operator/=(const double&);
  point operator/(const double&) const;
  void operator*=(const double&);
  point operator*(const double&) const;
  static double dot(point,point);
  static point cross(point,point);
  static bool near(double,double);
  static double len(point);
  point(double x, double y, double z) : x{x},y{y},z{z} {}
  void print();
};
