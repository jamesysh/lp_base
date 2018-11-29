#include <iostream>
#include <cmath>
#include <algorithm>
#include "geometry_pellet.h"

PelletLayer::PelletLayer(): xcen(0),ycen(0),zcen(0),innerradius(0.2), outerradius(0.24){}

bool PelletLayer::operator()(double x, double y, double z) const{
	double r2=(x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen);
	return (innerradius*innerradius<r2 && r2<outerradius*outerradius);
}

void PelletLayer::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax){
	xmin = xcen-outerradius;
	xmax = xcen+outerradius;
	ymin = ycen-outerradius;
	ymax = ycen+outerradius;
	zmin = zcen-outerradius;
	zmax = zcen+outerradius;
}

MultiPelletLayer::~MultiPelletLayer(){
    	delete[] xcen;
        delete[] ycen;
        delete[] zcen;
        delete[] innerradius;
        delete[] outerradius;
}

bool MultiPelletLayer::operator()(double x, double y, double z) const{
	for(int i=0;i<NumberofPellet;i++){
        	double r2=(x-xcen[i])*(x-xcen[i])+(y-ycen[i])*(y-ycen[i])+(z-zcen[i])*(z-zcen[i]);
        	if(innerradius[i]*innerradius[i]<r2 && r2<outerradius[i]*outerradius[i])
		return true;
	}
	return false;
}

void MultiPelletLayer::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax){
        xmin = xcen[0]-outerradius[0];
        xmax = xcen[0]+outerradius[0];
        ymin = ycen[0]-outerradius[0];
        ymax = ycen[0]+outerradius[0];
        zmin = zcen[0]-outerradius[0];
        zmax = zcen[0]+outerradius[0];

	for(int i=1;i<NumberofPellet;i++){
		        xmin=std::min(xmin,xcen[i]-outerradius[i]);
                xmax=std::max(xmax,xcen[i]+outerradius[i]);
                ymin=std::min(ymin,ycen[i]-outerradius[i]);
                ymax=std::max(ymax,ycen[i]+outerradius[i]);
                zmin=std::min(zmin,zcen[i]-outerradius[i]);
                zmax=std::max(zmax,zcen[i]+outerradius[i]);
	}
}



void MultiPelletLayer::initPellet(Initializer* init){
   
    NumberofPellet = init->getNumberofPellet();
    std::cout<<"Num"<<NumberofPellet<<std::endl;
    xcen = init->getPelletPositionX();
    ycen = init->getPelletPositionY();
    zcen = init->getPelletPositionZ();
    innerradius = init->getPelletInnerRadius();
    outerradius = init->getPelletRadius();
    }

