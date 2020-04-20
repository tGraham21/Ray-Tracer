// Assignment 1c Ray Caster - 5607
// Trevor Graham graha662


#define _USE_MATH_DEFINES
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
#include <regex>
#include <cmath>
#include "vector.cc"


using namespace std;

// Global declorations
Point eye;
Vector viewdir;
Vector updir;
double hfov;
int width, height;
int mat_index = -1;
int texture_index = -1;
int text_flag;
Color bkgcolor;
std::vector <Mtlcolor> mats;
std::vector <Sphere> spheres;
std::vector <Cylinder> cylinders;
std::vector <Light> lights;
std::vector <Point> verticies;
std::vector <Vector> surfaceNormals;
std::vector <Point> textureCoords;
std::vector <Face> faces;
std::vector <std::vector<std::vector<Color>>> textures;
// values contains color values of PPM of input textures
std::vector<std::vector<Color>> values;


void parseFace(std::stringstream& ss, Face& f){
  // 2 spaces between values 
  std::regex verts("\\d+\\s\\s\\d+\\s\\s\\d+"); // 3 verticies
  std::regex norms("\\d+\\/\\/\\d+\\s\\s\\d+\\/\\/\\d+\\s\\s\\d+\\/\\/\\d+");// verticies and normals v1//vn1 v2//vn2
  std::regex texts("\\d+\\/\\d+\\s\\s\\d+\\/\\d+\\s\\s\\d+\\/\\d+"); // verticies and texture coords v1/vt1 v2/vt2
  std::regex all("\\d+\\/\\d+\\/\\d+\\s\\s\\d+\\/\\d+\\/\\d+\\s\\s\\d+\\/\\d+\\/\\d+"); // verticies, texture coords and normals v1/vt1/vn1 v2/vt2/vn2
  
  // get values from verticies
  string in = ss.str();
  in = in.substr(3,in.length()-4); 

  if(std::regex_match(in,verts)){
    int x,y,z;
    ss >> x >> y >> z;
    f.v1 = verticies[x-1];
    f.v2 = verticies[y-1];
    f.v3 = verticies[z-1];
  }

  else if(std::regex_match(in,norms)){
    int x, y, z, nx, ny, nz;
    int pos = 0;

    while(pos < (int)(in.length() - 1)){
      if(in[pos] == '/' && in[pos+1] == '/'){
        // modify in for stringstream parsing
        in.replace(pos,1," ");
        in.erase(pos+1,1);
      }
      pos++;
    }

    ss = stringstream(in);
    ss >> x >> nx >> y >> ny >> z >> nz;

    f.v1 = verticies[x-1];
    f.v2 = verticies[y-1];
    f.v3 = verticies[z-1];

    f.n1 = surfaceNormals[nx-1];
    f.n2 = surfaceNormals[ny-1];
    f.n3 = surfaceNormals[nz-1];
  }

  else if(std::regex_match(in,texts)){
    int x, y, z, tx, ty, tz;
    int pos = 0;

    while(pos < (int)in.length()){
      if( in[pos] == '/' ){
        // modify in for stringstream parsing
        in.replace(pos, 1, " ");
      }
      pos++;
    }

    ss = stringstream(in);
    ss >> x >> tx >> y >> ty >> z >> tz;
    f.v1 = verticies[x-1];
    f.v2 = verticies[y-1];
    f.v3 = verticies[z-1];

    f.t1 = textureCoords[tx-1];
    f.t2 = textureCoords[ty-1];
    f.t3 = textureCoords[tz-1];
  }

  else if(std::regex_match(in, all)){
    int x, y, z, tx, ty, tz, nx, ny, nz;
    int i = 0;

    while(i < (int)in.length()){
      // modify for stringstream parsing
      if(in[i] == '/'){
        in.replace(i, 1, " ");
      }
      i++;
    }

    ss = stringstream(in);
    ss >> x >> tx >> nx >> y >> ty >> ny >>
    z >> tz >> nz;
    f.v1 = verticies[x-1];
    f.v2 = verticies[y-1];
    f.v3 = verticies[z-1];

    f.t1 = textureCoords[tx-1];
    f.t2 = textureCoords[ty-1];
    f.t3 = textureCoords[tz-1];

    f.n1 = surfaceNormals[nx-1];
    f.n2 = surfaceNormals[ny-1];
    f.n3 = surfaceNormals[nz-1];

  }
}

void readTexture(string text){

  std::ifstream inputTexture;
  inputTexture.open(text);

  // read in first line of input texture
  stringstream ss;
  string line;
  getline(inputTexture, line);

  // start contains the PPM identifier
  string start;
  int width, height, scale;

  ss = stringstream(line);
  ss >> start >> width >> height >> scale;

  for(int j = 0; j < height; j ++){
   std::vector<Color>row ;
    for(int i = 0; i < width; i++){
      getline(inputTexture, line);
      ss = stringstream(line);
      Color new_col = Color();
      ss >> new_col.r >> new_col.g >> new_col.b;
      row.push_back(new_col);
    }
    values.push_back(row);
  }

  inputTexture.close();
}

int parseFile(ifstream& in){
  string line;
  stringstream ss;
  string id;
  
  while(getline(in,line)){
    ss = stringstream(line);
    ss >> id;
    check:
    if(ss){
      if(!id.compare("eye")){
        ss >> eye.x;
        ss >> eye.y;
        ss >> eye.z;
      }
      else if(!id.compare("viewdir")){
        ss >> viewdir.x;
        ss >> viewdir.y;
        ss >> viewdir.z;
      }
      else if(!id.compare("updir")){
        ss >> updir.x;
        ss >> updir.y;
        ss >> updir.z;
      }
      else if(!id.compare("hfov")){
        ss >> hfov;
      }
      else if(!id.compare("imsize")){
        ss >> width;
        ss >> height;
      }
      else if(!id.compare("bkgcolor")){
        ss >> bkgcolor.r;
        ss >> bkgcolor.g;
        ss >> bkgcolor.b;
      }
      else if(!id.compare("mtlcolor")){
        Mtlcolor mtlcolor = Mtlcolor();
        ss >> mtlcolor.odr;
        ss >> mtlcolor.odg;
        ss >> mtlcolor.odb;
        ss >> mtlcolor.osr;
        ss >> mtlcolor.osg;
        ss >> mtlcolor.osb;
        ss >> mtlcolor.ka;
        ss >> mtlcolor.kd;
        ss >> mtlcolor.ks;
        ss >> mtlcolor.n;
		    ss >> mtlcolor.alpha;
		    ss >> mtlcolor.eta;
        mat_index++;
        text_flag = 0;
        mats.push_back(mtlcolor);
      }
      else if(!id.compare("sphere")){
        Sphere sphere = Sphere();
        ss >> sphere.cx;
        ss >> sphere.cy;
        ss >> sphere.cz;
        ss >> sphere.radius;
        sphere.mat_index = mat_index;
        sphere.text_index = texture_index;
        if(text_flag){
          sphere.flag = 1;
        }
        else{
          sphere.flag = 0;
        }
        spheres.push_back(sphere);
      }
      else if(!id.compare("cylinder")){
        Cylinder cylinder = Cylinder();
        ss >> cylinder.cx;
        ss >> cylinder.cy;
        ss >> cylinder.cz;
        ss >> cylinder.dx;
        ss >> cylinder.dy;
        ss >> cylinder.dz;
        ss >> cylinder.radius;
        ss >> cylinder.length;
      }
      else if(!id.compare("light")){
        Light light = Light();
        ss >> light.x;
        ss >> light.y;
        ss >> light.z;
        ss >> light.w;
        ss >> light.r;
        ss >> light.g;
        ss >> light.b;
        lights.push_back(light);
      }
      else if(!id.compare("v")){  
         Point vertex = Point();
         ss >> vertex.x;
         ss >> vertex.y;
         ss >> vertex.z;
         verticies.push_back(vertex);
      }
      else if(!id.compare("vn")){
        Vector norm = Vector();
        ss >> norm.x;
        ss >> norm.y;
        ss >> norm.z;
        surfaceNormals.push_back(norm);
      }
      else if(!id.compare("vt")){
        Point text = Point();
        ss >> text.x;
        ss >> text.y;
        ss >> text.z;
        textureCoords.push_back(text);
      }
      else if(!id.compare("f")){
        Face face =  Face();
        face.mat_index = mat_index;
        face.text_index = texture_index;
        if(text_flag){
          face.flag = 1;
        }
        else{
          face.flag = 0;
        }
        parseFace(ss, face);
        faces.push_back(face);
      }
      else if(!id.compare("texture")){
        string text;
        ss >> text;
    
        readTexture(text);
        // add each texture to textures
        textures.push_back(values);
        texture_index ++;
        text_flag = 1;
      }
    }
    // check to see if there is another line
    else if (getline(in,line)){
      ss = stringstream(line);
      ss >> id;
      // jump to check to continue parsing
      goto check; 
    }
    else{
      return 1;
    }
  }
  return 0;
}

// function determines if the ray intersects the light source or needs a shadow
bool traceShadowRay(Ray ray, double dist){
  double min_distance = 10000;
  double b, c;
  double pos_point, neg_point;
  Sphere inter_sphere = Sphere();

  for(int i = 0; i< (int) spheres.size(); i++){
    // check all spheres for intersection
    b = 2 * ((ray.dx * (ray.x - spheres[i].cx)) + (ray.dy * (ray.y - spheres[i].cy)) +
    (ray.dz * (ray.z - spheres[i].cz)));
    c = pow((ray.x - spheres[i].cx),2) + pow((ray.y - spheres[i].cy),2) +
     pow((ray.z - spheres[i].cz),2) - pow(spheres[i].radius,2);

    double disc = pow(b,2) - (4 * c);
    pos_point = ((-1 * b) + (sqrt(disc))) /(2);
    neg_point = ((-1 * b) - (sqrt(disc))) /(2);

    if(disc >= 0){
      if(neg_point <= pos_point && neg_point > 0.001){
        if(neg_point < min_distance){
          min_distance = neg_point;
          inter_sphere = spheres[i];
        }
      }
      else if(pos_point > 0.001){
        if(pos_point < min_distance){
          min_distance = pos_point;
          inter_sphere = spheres[i];
        }
      }
    }
  }

  if(min_distance > -1 && min_distance < dist){
    // return true if there is an intersection and should be a shadow
    return true;
  }

  return false;
}

Color traceRay(Ray ray, int count, float eta_I){
  Vector N;
  Vector n, n_temp;
  float t = 0;
  Color pos_col;
  double min_distance = 10000;
  int min_index = -1;
  double b, c;
  double pos_point, neg_point;
  Point inter_point = Point();
  Sphere inter_sphere = Sphere();
  Face inter_face = Face();
  float alpha, beta, gamma;
  std::vector <Vector> L_list;

  // set default values for col if none are specified
  Mtlcolor col = Mtlcolor();
  col.osr = 1 ; col.osg = 1; col.osb = 1;
  col.ka = .3; col.kd = .5; col.ks = .2;
  col.n = 2; col.alpha = 1; col.eta = 1;
  
  Point rOrigin = Point();
  rOrigin.x = ray.x; rOrigin.y = ray.y; rOrigin.z =ray.z;

  for(int i = 0; i< (int)spheres.size(); i++){ 
    // check spheres for intersection
    b = 2 * ((ray.dx * (ray.x - spheres[i].cx)) + (ray.dy * (ray.y - spheres[i].cy)) +
    (ray.dz * (ray.z - spheres[i].cz)));
    c = pow((ray.x - spheres[i].cx),2) + pow((ray.y - spheres[i].cy),2) +
     pow((ray.z - spheres[i].cz),2) - pow(spheres[i].radius,2);
    double disc = pow(b,2) - (4 * c);
    pos_point = ((-1 * b) + (sqrt(disc))) /(2);
    neg_point = ((-1 * b) - (sqrt(disc))) /(2);

    if(disc >= 0){
      if(neg_point <= pos_point && neg_point > 0.01){
        if(neg_point < min_distance){
          min_distance = neg_point;
          min_index = i;
          inter_sphere = spheres[i];
        }
      }
      else if(pos_point > 0.01){
        if(pos_point < min_distance){
          min_distance = pos_point;
          min_index = i;
          inter_sphere = spheres[i];
        }
      }
    }
  }

    // define curr_light vector depending on what type of light source and add to L_list
    for(int i = 0; i < (int) lights.size(); i++){
      Vector curr_light = Vector();
      if(lights[i].w == 0){ 
        // directional light 
        curr_light.x =  lights[i].x; curr_light.y = lights[i].y; curr_light.z =lights[i].z;
      }
      else if(lights[i].w == 1){
         // point light
        curr_light.x = lights[i].x - inter_point.x; curr_light.y = lights[i].y - inter_point.y; curr_light.z = lights[i].z - inter_point.z;
      }
      L_list.push_back(normalize(curr_light));
    }
    
    if(min_index > -1){

      Vector dir = Vector(); 
      Point p = Point();
      dir.x = ray.dx; dir.y = ray.dy; dir.z = ray.dz;
      p.x = ray.x ;  p.y = ray.y ; p.z = ray.z ;
      // calculate intersection point 
      inter_point = add(mult(min_distance,dir), p);

      Point center = Point();
      center.x = inter_sphere.cx; center.y = inter_sphere.cy; center.z = inter_sphere.cz;
      N = normalize(sub(inter_point,center));
      Vector E = normalize(sub(eye, inter_point));

      if(inter_sphere.mat_index >= 0){
       // if there is a material defined apply that, otherwise default
       col =  mats[inter_sphere.mat_index];
    }
      if(inter_sphere.flag){
        // if there is a texture defined for sphere read texture to determine color
        float theta = std::atan2(-1 * N.z, - 1 * N.y);
        float phi = acos(- 1 * N.x);
        float sv = phi/ M_PI;
        float su = .5 + theta/(2 * M_PI);
        int sx = floor(su * ((textures[inter_sphere.text_index].size()) - 1 ));
        int sy = floor(sv * ((textures[inter_sphere.text_index][0].size()) - 1 ));
          col.odr = textures[inter_sphere.text_index][sx][sy].r;
          col.odg = textures[inter_sphere.text_index][sx][sy].g;
          col.odb = textures[inter_sphere.text_index][sx][sy].b;
      }
      pos_col = Color();
      // start with ambient light
      pos_col.r = (col.ka * col.odr); 
      pos_col.g = (col.ka * col.odg);
      pos_col.b = (col.ka * col.odb);

      for(int i = 0; i < (int)lights.size(); i++){
        int s = 1;
        Vector shadow_dir = L_list[i];
        Ray shadow_ray = Ray(inter_point.x, inter_point.y, inter_point.z, shadow_dir.x, shadow_dir.y, shadow_dir.z );
        double dist = sqrt(pow(lights[i].x,2) + pow(lights[i].y,2) + pow(lights[i].z,2));

        if(traceShadowRay(shadow_ray,dist) == true){
          // set to 0 if point is determined to be within a shadow
          s = 1 ;
        }

        // halfway vector for calculations
        Vector H = normalize(add(L_list[i], E)); 

        // add color calulations for diffuse and specular light
        // account for shadow
        pos_col.r +=  s *(lights[i].r * ((col.kd * col.odr * (max(0.0,dot(N,L_list[i])))) + (col.ks * col.osr * (pow(max(0.0,dot(N,H)),col.n)))));
        pos_col.g +=  s * (lights[i].g * ((col.kd * col.odg * (max(0.0,dot(N,L_list[i])))) + (col.ks * col.osg * (pow(max(0.0,dot(N,H)),col.n)))));
        pos_col.b +=  s * (lights[i].b * ((col.kd * col.odb * (max(0.0,dot(N,L_list[i])))) + (col.ks * col.osb * (pow(max(0.0,dot(N,H)),col.n)))));
      }
     
      
	  // reflection and refraction calculations
    if (count < 6){
      count ++;
      Vector refractN = N;
      Vector iDir = normalize(sub(rOrigin, inter_point));
      iDir.x = ray.dx * -1; iDir.y = ray.dy * -1; iDir.z = ray.dz * -1;  

      // Direction for R
      Vector rDirSphere = normalize(sub(mult(2 * max(0.0,dot(N,iDir)),N),iDir));
      Ray R = Ray(inter_point.x, inter_point.y, inter_point.z, rDirSphere.x, rDirSphere.y, rDirSphere.z);

      // accumulated color of reflections
      Color reflectColor = traceRay(R, count, 1);
        
      // Fresnel Coefficient Calculations 
      float fZero = pow(((col.eta - eta_I) / (col.eta + eta_I)),2);
      float fR = fZero + ((1 - fZero) * pow((1 - max(0.0,(dot(N,iDir)))),5));
      
      // weight R by Fresnel Coefficient
      pos_col.r += (fR * reflectColor.r);   
      pos_col.g += (fR * reflectColor.g);
      pos_col.b += (fR * reflectColor.b);

      // flip N if dot(N,I) is negative
      // intersection is inside of the sphere
      if(dot(refractN, iDir) < 0){
        refractN = mult(-1, refractN);
      }

      Vector negN = mult(-1, refractN);
      float eta_T = col.eta;

      // calculations for refracted ray T
      float tFront = sqrt(1 - (pow((eta_I/eta_T), 2) * (1 - pow(dot(refractN,iDir),2))));
      Vector tValue1 = mult(tFront, negN);
      Vector tValue2  = sub(mult(dot(refractN,iDir),refractN),iDir);
      tValue2 = mult((eta_I/eta_T), tValue2);
      Vector tDirSphere = add(tValue1, tValue2);

      Ray T = Ray(inter_point.x, inter_point.y, inter_point.z, tDirSphere.x, tDirSphere.y, tDirSphere.z);
      Color transparencyColor = traceRay(T, count, eta_T);
        
      pos_col.r += ((1 - fR) * (1 - col.alpha) * transparencyColor.r);
      pos_col.g += ((1 - fR) * (1 - col.alpha) * transparencyColor.g);
      pos_col.b += ((1 - fR) * (1 - col.alpha) * transparencyColor.b);

      // cap values at 1
      if(pos_col.r > .999){
        pos_col.r = 1;
      }
      if (pos_col.g > .999){
        pos_col.g = 1;
      }
      if (pos_col.b > .999){
        pos_col.b = 1;
      }       
        return pos_col;
      }
    else{
      if(pos_col.r > 1){
        pos_col.r = 1;
      }
      if (pos_col.g > 1){
        pos_col.g = 1;
      }
      if (pos_col.b > 1){
        pos_col.b = 1;
      }
        return pos_col;
    }

    }
  
  Point intersect = Point();
  float min_dist_face = 10000;

  // check for intersection with faces
  for(int i = 0; i < (int)faces.size(); i++){ 
    Point p0 = faces[i].v1; Point p1 = faces[i].v2; Point p2 = faces[i].v3;
      Vector e1 = sub(p1,p0); Vector e2 = sub(p2, p0);
      n_temp = cross(e1, e2); 
      n_temp = normalize(n_temp);
      float d = (n_temp.x * p0.x) +(n_temp.y *p0.y) + (n_temp.z * p0.z);
      d = -1 * d;
     float denom = (n_temp.x * ray.dx) +(n_temp.y *ray.dy) + (n_temp.z * ray.dz); 
     
     if(denom != 0 ){ 
      t =  (-1* ((n_temp.x * ray.x) + (n_temp.y * ray.y) + (n_temp.z *ray.z) + d)/denom);
      Vector rDir = Vector();
      rDir.x = ray.dx; rDir.y = ray.dy; rDir.z = ray.dz;
      // calculate intersection point
      intersect = add(mult(t, rDir), rOrigin);
      if( t > 0.00001 && t < min_dist_face){

       // check if intersection is within face
       Vector ep = sub(intersect, p0);
       float d11 = dot(e1,e1); float d12 = dot(e1,e2); float d22 = dot(e2,e2);
       float dp1 = dot(ep,e1); float dp2 = dot(ep,e2);
       float D_base = (d11 * d22) - (d12 * d12);
       float D_beta = (d22 * dp1) - (d12 * dp2);
       float D_gamma = (d11 * dp2) - (d12 * dp1);
       beta = D_beta/D_base;
       gamma = D_gamma/ D_base;
       alpha = 1 - (beta + gamma);
       if((alpha >= 0 && alpha <=1) && (beta >= 0 && beta <=1) && (gamma >= 0 && gamma <= 1)){
         min_dist_face = t;
         inter_face = faces[i];
         inter_face.alpha = alpha; inter_face.beta = beta; inter_face.gamma = gamma;
         
         // set n if there is a vertex normal
         if(mag(faces[i].n1) + mag(faces[i].n2) + mag(faces[i].n3) != 0){ 
           n = add(add(mult(alpha, faces[i].n1), mult(beta, faces[i].n2)), mult(gamma, faces[i].n3));
           n = normalize(n);
         }
         else{
           n = n_temp;
         }

       }
      }
    }
  }
    float u = 0; float v = 0;
    float x = 0; float y = 0;
    if( min_dist_face <  min_distance && min_dist_face > 0.0001){
      if(inter_face.mat_index >= 0){
      col =  mats[inter_face.mat_index];
      }

      // look up nearest pixel in texture to determine the color 
      if(inter_face.flag && (inter_face.alpha >= 0 && inter_face.alpha <=1) && (inter_face.beta >= 0 && inter_face.beta <=1) && (inter_face.gamma > 0 && inter_face.gamma <= 1)){
        u = (inter_face.alpha * inter_face.t1.x) + (inter_face.beta * inter_face.t2.x) + (inter_face.gamma * inter_face.t3.x);
        v =  (alpha * inter_face.t1.y) + (beta * inter_face.t2.y) + (gamma * inter_face.t3.y);
        x = floor(u * ((textures[inter_face.text_index].size()) - 1 ));
        y = floor(v * ((textures[inter_face.text_index][0].size()) - 1 ));
        col.odr = textures[inter_face.text_index][x][y].r;
        col.odg = textures[inter_face.text_index][x][y].g;
        col.odb = textures[inter_face.text_index][x][y].b;
      }

      // start with ambient light
      pos_col.r = (col.ka * col.odr); 
      pos_col.g = (col.ka * col.odg);
      pos_col.b = (col.ka * col.odb);
      for(int i = 0; i < (int)lights.size(); i++){
        int s = 1;
        Vector shadow_dir = L_list[i];
        Ray shadow_ray = Ray(inter_point.x, inter_point.y, inter_point.z, shadow_dir.x, shadow_dir.y, shadow_dir.z );
        double dist = sqrt(pow(lights[i].x,2) + pow(lights[i].y,2) + pow(lights[i].z,2));

        if(traceShadowRay(shadow_ray,dist) == true){
          s = 0; 
        }
        Vector E = normalize(sub(eye, intersect));
        // halfway vector for use in light calculations
        Vector H = normalize(add(L_list[i], E)); 
        // add in diffuse and specular calculations
        pos_col.r +=  s *(lights[i].r * ((col.kd * col.odr * (max(0.0,dot(n,L_list[i])))) + (col.ks * col.osr * (pow(max(0.0,dot(n,H)),col.n)))));
        pos_col.g +=  s * (lights[i].g * ((col.kd * col.odg * (max(0.0,dot(n,L_list[i])))) + (col.ks * col.osg * (pow(max(0.0,dot(n,H)),col.n)))));
        pos_col.b +=  s * (lights[i].b * ((col.kd * col.odb * (max(0.0,dot(n,L_list[i])))) + (col.ks * col.osb * (pow(max(0.0,dot(n,H)),col.n)))));
      }

      // reflection and refraction calculations for faces
      if (count < 6){
        count ++;
        Vector refractN = n;
        Vector iDir = normalize(sub(rOrigin, intersect));   

        // calculate direction of R vector
        Vector rDirFace = normalize(sub(mult(2 * max(0.0,dot(n,iDir)),n),iDir));
        Ray R = Ray(intersect.x, intersect.y, intersect.z, rDirFace.x, rDirFace.y, rDirFace.z);
        
        // accumulation of reflected color
        Color reflectColor = traceRay(R, count, 1);
        
        // Fresnel Coefficient Calculations 
        float fZero = pow(((col.eta - eta_I) / (col.eta + eta_I)),2);
        float fR = fZero + ((1 - fZero) * pow((1 - max(0.0,(dot(n,iDir)))),5));

        
        // weight R by Fresnel Coefficient
        pos_col.r += (fR * reflectColor.r);   
        pos_col.g += (fR * reflectColor.g);
        pos_col.b += (fR * reflectColor.b);

        // flip N if dot(N,I) is negative
        // intersection is inside of the sphere
        if(dot(refractN, iDir) < 0){
          refractN = mult(-1, refractN);
        }

        // calculations for refraction ray T
        Vector negN = mult(-1, refractN);
        float eta_T = col.eta;
      
        float tFront = sqrt(1 - (pow((eta_I/eta_T), 2) * (1 - pow(dot(refractN,iDir),2))));
        Vector tValue1 = mult(tFront, negN);
        Vector tValue2  = sub(mult(dot(refractN,iDir),refractN),iDir);
        tValue2 = mult((eta_I/eta_T), tValue2);
        Vector tDirFace = add(tValue1, tValue2);
        Ray T = Ray(inter_point.x, inter_point.y, inter_point.z, tDirFace.x, tDirFace.y, tDirFace.z);
        
        // shift slightly to avoid ray collisions and spotted appearance
        T.x = intersect.x + .001 ; T.y = intersect.y + .001; T.z = intersect.z + .001; 
        Color transparencyColor = traceRay(T, count, eta_T);
        
        pos_col.r += ((1 - fR) * (1 - col.alpha) * transparencyColor.r);
        pos_col.g += ((1 - fR) * (1 - col.alpha) * transparencyColor.g);
        pos_col.b += ((1 - fR) * (1 - col.alpha) * transparencyColor.b);



        // cap values at 1
        if(pos_col.r > 1){
          pos_col.r = 1;
        }
        if (pos_col.g > 1){
          pos_col.g = 1;
        }
        if (pos_col.b > 1){
          pos_col.b = 1;
        }
        
        
        return pos_col;
      }
      else{
         if(pos_col.r > 1){
          pos_col.r = 1;
        }
        if (pos_col.g > 1){
          pos_col.g = 1;
        }
        if (pos_col.b > 1) {
          pos_col.b = 1;
        }

        return pos_col;
      }
    }
 return bkgcolor;
  
}

int main (int argc, char **argv){
  std::string fileName = argv[1];
  std::ifstream inputFile;
  inputFile.open(fileName);
 
  // argument error handling
  if (argc != 2){
    std::cout << "* Incorrect args * USAGE: ./a.out inputfile" << std::endl;
    exit(0);
  }
  else if (!inputFile){
   std::cout << "Could not open file" << std::endl;
   exit(0);
  }
  else if(parseFile(inputFile)){
   std::cout << "Error parsing input file" << std::endl;
   exit(0);

  }

  int d = 2;
  double aspect = width / height;
  Color pixels[width][height];
  hfov =  hfov * (M_PI/180);

  int window_width = round((2 * d) * tan((hfov/2)));
  int window_height = window_width / aspect;

  Vector u = cross(viewdir, updir);
  u = normalize(u);

  Vector v = cross(u, viewdir);
  v= normalize(v);

  Vector w = mult(-1, viewdir);
  w = normalize(w);


  Vector n = normalize(viewdir);

 // Define the All the corners in the viewing window ul, ur, ll, lr
  Vector ulVector = add(sub(mult(d,n), (mult((window_width/2),u))),mult(window_height/2, v));
  Point ul = eye;
  ul.x = ul.x + ulVector.x;
  ul.y = ul.y + ulVector.y;
  ul.z = ul.z + ulVector.z;

  Vector urVector = add(add(mult(d,n), (mult((window_width/2),u))),mult(window_height/2, v));
  Point ur = eye;
  ur.x = ur.x + urVector.x;
  ur.y = ur.y + urVector.y;
  ur.z = ur.z + urVector.z;

  Vector llVector = sub(sub(mult(d,n), (mult((window_width/2),u))),mult(window_height/2, v));
  Point ll = eye;
  ll.x = ll.x + llVector.x;
  ll.y = ll.y + llVector.y;
  ll.z = ll.z + llVector.z;

  Vector lrVector = sub(add(mult(d,n), (mult((window_width/2),u))),mult(window_height/2, v));
  Point lr = eye;
  lr.x = lr.x + lrVector.x;
  lr.y = lr.y + lrVector.y;
  lr.z = lr.z + lrVector.z;

// Define offset values for viewing window
 Vector ch_offset = div((sub(ur,ul)), (2 * width));
 Vector cv_offset = div((sub(ll,ul)), (2 * height));
 Vector h_offset = div((sub(ur,ul)), (width));
 Vector v_offset = div((sub(ll,ul)), (height));

  for(int j = 0; j < height; j++){
    for(int i = 0; i < width; i++){
      // create ray through pixel
      Point ray_point = add(add(add(add(mult(i,h_offset),mult(j,v_offset)), ch_offset), cv_offset)  ,ul);
     
      // eye is origin
      Vector norm_dir = Vector();
      norm_dir = normalize(sub(ray_point, eye));
      Ray ray = Ray(eye.x, eye.y, eye.z, norm_dir.x, norm_dir.y, norm_dir.z);
      
      pixels[i][j] = traceRay(ray, 0, 1);
    }
  }

  // output to PPM file
  std::string output = "output_image.ppm";
  std::ofstream outputFile;
  outputFile.open(output);
  outputFile << "P3" << std::endl;
  outputFile <<width << " "<< height << std::endl;
  outputFile <<255<< std::endl;
  for(int j = 0; j < height; j++){
    for(int i = 0; i < width; i++){
      char buff[100];
      int r,g,b;
      // scale 0-255 instead of 0-1
      r = (int)255 *pixels[i][j].r;
      g = (int)255 *pixels[i][j].g;
      b = (int)255 * pixels[i][j].b;
      snprintf(buff, sizeof(buff), "%d %d %d", r, g, b);
      std::string output = buff;
      outputFile << output << std::endl;
    }

  }
  //close files
  outputFile.close();
  inputFile.close();

}
