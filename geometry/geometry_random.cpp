#include "geometry_random.h"
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

Uniform3D::Uniform3D():length(1.0){xCen={0,3,0,3}, yCen={0,0,3.5,3.5}, zCen={0,0,0,0};}

bool Uniform3D::operator()(double x, double y, double z) const {
        return (x<*max_element(xCen.begin(),xCen.end())+length/2 && x>*min_element(xCen.begin(),xCen.end())-length/2 && y<*max_element(yCen.begin(),yCen.end())+length/2 && y>*min_element(yCen.begin(),yCen.end())-length/2 && z<*max_element(zCen.begin(),zCen.end())+length/2 && z>*min_element(zCen.begin(),zCen.end())-length/2);
}

void Uniform3D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = *min_element(xCen.begin(),xCen.end())-length/2;
        xmax = *max_element(xCen.begin(),xCen.end())+length/2;
        ymin = *min_element(yCen.begin(),yCen.end())-length/2;
        ymax = *max_element(yCen.begin(),yCen.end())+length/2;
        zmin = *min_element(zCen.begin(),zCen.end())-length/2;
        zmax = *max_element(zCen.begin(),zCen.end())+length/2;
}

void Uniform3D::randomlocation(double& x, double& y, double& z){
	
	int i=rand()%xCen.size();

	x = xCen[i] - length/2 + length*((double)rand()/(double)RAND_MAX);
        y = yCen[i] - length/2 + length*((double)rand()/(double)RAND_MAX);
        z = zCen[i] - length/2 + length*((double)rand()/(double)RAND_MAX);
}

Gaussian3D::Gaussian3D():sigma(1),radius(3),xCen(0),yCen(0),zCen(0){}

bool Gaussian3D::operator()(double x, double y, double z) const {
        return ((x-xCen)*(x-xCen)+(y-yCen)*(y-yCen)+(z-zCen)*(z-zCen)<radius*radius);
}

void Gaussian3D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = xCen-radius;
        xmax = xCen+radius;
        ymin = yCen-radius;
        ymax = yCen+radius;
        zmin = zCen-radius;
        zmax = zCen+radius;
}

void Gaussian3D::randomlocation(double& x, double& y, double& z){

	double r1,r2,w=10;
	while(w>=1.0){
		r1=2.0*((double)rand()/(double)RAND_MAX)-1.0;
		r2=2.0*((double)rand()/(double)RAND_MAX)-1.0;
		w=r1*r1+r2*r2;
	}
	w=sqrt((-2.0*log(w))/w);
	x=r1*w*sigma+xCen;
	y=r2*w*sigma+yCen;
	w=10;
        while(w>=1.0){
                r1=2.0*((double)rand()/(double)RAND_MAX)-1.0;
                r2=2.0*((double)rand()/(double)RAND_MAX)-1.0;
                w=r1*r1+r2*r2;
        }
        w=sqrt((-2.0*log(w))/w);
        z=r1*w*sigma+zCen;
}

MultiGaussian3D::MultiGaussian3D():sigma(0.1),radius(4){xCen={0.,0.,2.,2.,1.},yCen={0.,2.,0.,2.,1.},zCen={0.,0.,0.,0.,0};}

bool MultiGaussian3D::operator()(double x, double y, double z) const {
        return ((x-xCen[0])*(x-xCen[0])+(y-yCen[0])*(y-yCen[0])+(z-zCen[0])*(z-zCen[0])<radius*radius);
}

void MultiGaussian3D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = xCen[0]-radius;
        xmax = xCen[0]+radius;
        ymin = yCen[0]-radius;
        ymax = yCen[0]+radius;
        zmin = zCen[0]-radius;
        zmax = zCen[0]+radius;
}

void MultiGaussian3D::randomlocation(double& x, double& y, double& z){

	int i=rand()%xCen.size();

        double r1,r2,w=10;
        while(w>=1.0){
                r1=2.0*((double)rand()/(double)RAND_MAX)-1.0;
                r2=2.0*((double)rand()/(double)RAND_MAX)-1.0;
                w=r1*r1+r2*r2;
        }
        w=sqrt((-2.0*log(w))/w);
        x=r1*w*sigma+xCen[i];
        y=r2*w*sigma+yCen[i];
        w=10;
        while(w>=1.0){
                r1=2.0*((double)rand()/(double)RAND_MAX)-1.0;
                r2=2.0*((double)rand()/(double)RAND_MAX)-1.0;
                w=r1*r1+r2*r2;
        }
        w=sqrt((-2.0*log(w))/w);
        z=r1*w*sigma+zCen[i];
}
