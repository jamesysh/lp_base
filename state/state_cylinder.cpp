#include "state_cylinder.h"
#include <iostream>

CylinderState::CylinderState(){}

double CylinderState::pressure(double x, double y, double z) {
    double R = 83.14/2.014;
    double T = 10*11604.525;
    double den = density(x,y,z);
    m_fPressure = den*R*T;
    
    return m_fPressure;
}

double CylinderState::density(double x, double y, double z){
	if(x<-1)
        m_fDen = -1.044e-6*(x+1)*(x+5) + 1.e-10;
    else if(x>1)
        m_fDen = -1.044e-6*(x-1)*(x-5) +1.e-10;
    else
        m_fDen = 1.e-10;
    return m_fDen;
}

void CylinderState::velocity(double x, double y, double z, double& vX, double& vY,double& vZ){
	vX=0;
	vY=0;
	vZ=0;
}
