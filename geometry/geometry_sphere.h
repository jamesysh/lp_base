#ifndef __GEOMETRY_SPHERE_H
#define __GEOMETRY_SPHERE_H

#include "geometry.h"
#include "initializer.h"
#include <math.h>




class Sphere: public Geometry {
public:
	Sphere();
	virtual ~Sphere() {}
	virtual bool operator()(double x, double y, double z) const;
	virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double&zmax);
private:
	double xcen;
	double ycen;
	double zcen;
	double radius;
};

#endif
