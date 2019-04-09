#include "geometry_cylinder.h"
#include <iostream>
#include <cmath>


Cylinder::Cylinder():radius(3/2.),length(8),xCen(0),yCen(0),zCen(){}


bool Cylinder::operator()(double x, double y, double z) const {	
    
    if(x <= (xCen+length/2.) && x >= (xCen-length/2.) && (y*y+z*z)<radius*radius )

        return true;
        
    else

        return false;

}

void Cylinder::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-length/2.;
	xmax = xCen+length/2.;
	ymin = yCen-radius;
	ymax = yCen+radius;
	zmin = zCen-radius;
	zmax = zCen+radius;
}


