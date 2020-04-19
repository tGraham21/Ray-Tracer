// a1_struct.h - 5607
// Definitions of all the data structures used in assignment1.cc
typedef struct{
  double x, y, z;
}Point;

typedef struct{
  double x,y,z;
}Vector;

typedef struct{
  double cx, cy, cz;
  double radius;
  int mat_index;
  int flag; // 1 - text
  int text_index;
} Sphere;

typedef struct{
  double cx, cy, cz;
  double dx, dy, dz;
  double radius;
  double length;
  int mat_index;
}Cylinder;

typedef struct {
  double r,g,b;
}Color;

typedef struct{
  double odr,odg,odb;
  double osr, osg, osb;
  double ka,kd,ks;
  int n;
  double alpha, eta;
} Mtlcolor;

typedef struct{
  double x,y,z,w;
  double r,g,b;
}Light;

typedef struct{
  double x, y, z;
  double dx, dy, dz;
} Ray;

typedef struct{
  Point v1, v2, v3 ;
  Vector n1, n2, n3;
  Point t1, t2, t3;
  int mat_index;
  int flag; // 1 - text
  int text_index;
  float alpha,beta,gamma;
} Face;
