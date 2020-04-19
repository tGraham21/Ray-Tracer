// Vector Lib - 5607
// All vector math functions that are used in assignment 1

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <typeinfo>
#include "a1_struct.h"


Vector normalize(Vector v){
 double mag = sqrt((v.x * v.x) + (v.y* v.y) + (v.z * v.z));
 v.x = v.x /mag;
 v.y = v.y /mag;
 v.z = v.z / mag;
 return v;
}

Vector add(Vector u, Vector v){
  Vector vect = Vector();
  vect.x = u.x + v.x;
  vect.y = u.y + v.y;
  vect.z = u.z + v.z;

  return vect;

}

Point add(Vector v, Point p){
  Point point = Point();
  point.x = p.x + v.x;
  point.y = p.y + v.y;
  point.z = p.z + v.z;

  return point;
}
Vector sub(Vector u, Vector v){
  Vector vect = Vector();
  vect.x = u.x - v.x;
  vect.y = u.y - v.y;
  vect.z = u.z - v.z;

  return vect;

}

Vector sub(Point p, Point q){
  Vector vect = Vector();
  vect.x = p.x - q.x;
  vect.y = p.y - q.y;
  vect.z = p.z - q.z;

  return vect;
}

Vector mult(double s, Vector v){
  Vector vect = Vector();
  vect.x = s * v.x;
  vect.y = s * v.y;
  vect.z = s * v.z;

  return vect;
}


Vector cross(Vector u, Vector v){
  Vector vect = Vector();
  vect.x = (u.y * v.z) - (u.z * v.y);
  vect.y = (u.z * v.x) - (u.x * v.z);
  vect.z = (u.x * v.y) - (u.y * v.x);

  return vect;

}

Vector div (Vector v, double d){
  Vector vect = Vector();
  vect.x = v.x/d;
  vect.y = v.y/d;
  vect.z = v.z/d;

  return vect;

}

double mag (Vector v){
  double mag = sqrt((v.x * v.x) + (v.y* v.y) + (v.z * v.z));
  return mag;
}
double dot (Vector u, Vector v){
  return (u.x * v.x) + (u.y * v.y) + (u.z*v.z);
}
