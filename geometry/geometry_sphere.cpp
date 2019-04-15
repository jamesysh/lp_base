#include <iostream>
#include <cmath>
#include <algorithm>
#include "geometry_sphere.h"

Sphere::Sphere(): xcen(0),ycen(0),zcen(0),radius(1){}

bool Sphere::operator()(double x, double y, double z) const{
	double r2=(x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen);
	return ( r2<radius*radius);
}

void Sphere::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax){
	xmin = xcen-radius;
	xmax = xcen+radius;
	ymin = ycen-radius;
	ymax = ycen+radius;
	zmin = zcen-radius;
	zmax = zcen+radius;
}

