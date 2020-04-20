// a1_struct.h - 5607
// Definitions of all the data structures used in assignment1.cc
 struct Point{
  double x, y, z;
  Point(double xVal = 0, double yVal = 0, double zVal = 0){
    x = xVal;
    y = yVal;
    z = zVal;
  }
};

struct Vector{
  double x,y,z;
  Vector(double xVal = 0, double yVal = 0, double zVal = 0){
    x = xVal;
    y = yVal;
    z = zVal;
  }
};

typedef struct {
  double cx, cy, cz;
  double radius;
  int mat_index;
  int flag; // 1 - text
  int text_index;
}Sphere;

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

struct Ray{
  double x, y, z;
  double dx, dy, dz;
  Ray(double xVal, double yVal, double zVal, double dirX, double dirY, double dirZ){
    x = xVal; y = yVal; z = zVal;
    dx = dirX; dy = dirY; dz = dirZ;
  }
};

typedef struct{
  Point v1, v2, v3 ;
  Vector n1, n2, n3;
  Point t1, t2, t3;
  int mat_index;
  int flag; // 1 - text
  int text_index;
  float alpha,beta,gamma;
} Face;
