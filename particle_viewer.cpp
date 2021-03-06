#include "particle_viewer.h"
#include "particle_data.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip> //setw

//#define DEBUG_LW_V
//#define DEBUG_LW_S
//#define NEIGHOFPARTICLE
//#define DEBUG_SWITCH

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
// Start of ParticleViewer
////////////////////////////////////////////////////////////////////////////////////////


string ParticleViewer::rightFlush(size_t writeStep, size_t numDigits) {
	
	assert(pow(10,numDigits) > writeStep);

	string result;

	if(writeStep == 0) numDigits--;
	for(size_t i=writeStep; i!=0; i /= 10) numDigits--;

	for( ; numDigits>0; numDigits--) result.push_back('0');
	
	result += to_string(writeStep); 
	
	return result;

}


////////////////////////////////////////////////////////////////////////////////////////
// End of ParticleViewer
////////////////////////////////////////////////////////////////////////////////////////


double VTKParticleViewer::exact_density(double x, double y, double z, double time)
{
if(exactsolution=="simplewave"){
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
        double v_post=0.92745;
        double time1=3.0+time;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);
        double x1=-c_left*time1;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time1;
	if(x<x1)
		return rho_left;
	else if (x<x2)
	{
                        double c=mu*mu*(-x/time1)+(1-mu*mu)*c_left;
                        return rho_left*pow((c/c_left),(2/(gamma-1)));
	}
	else
	{
                        double c=mu*mu*(-x2/time1)+(1-mu*mu)*c_left;
                        return rho_left*pow((c/c_left),(2/(gamma-1)));
	}
}
else if(exactsolution=="gresho"){
	return 1;
}
else if(exactsolution=="shocktube"){
//	x=x-0.005;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
//        double u_left=0.0;

        double rho_right=1.0;

//        double rho_right=0.125;
//        double p_right=0.1;
//        double u_right=0.0;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);

        double rho_middle=0.62847;
        double rho_post=2.8803;
//        double p_post=0.52191;
        double v_post=0.52481;
        double v_shock=0.80392;


//        double rho_middle=0.42632;
//        double rho_post=0.26557;
//        double p_post=0.30313;
//        double v_post=0.92745;
//        double v_shock=1.75216;

        double x1=-c_left*time;
        double x3=v_post*time;
        double x4=v_shock*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;

	if(x<x1)
		return rho_left;
	if(x<x2)
	{
                double c=mu*mu*(-x/time)+(1-mu*mu)*c_left;
                return rho_left*pow((c/c_left),(2/(gamma-1)));
	}
	if(x<x3)
		return rho_middle;
	if(x<x4)
		return rho_post;
	return rho_right;
}
else if(exactsolution=="sod"){
//      x=x-0.005;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
        double u_left=0.0;

        double rho_right=0.125;
        double p_right=0.1;
        double u_right=0.0;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);

        double rho_middle=0.42632;
        double rho_post=0.26557;
        double p_post=0.30313;
        double v_post=0.92745;
        double v_shock=1.75216;

        double x1=-c_left*time;
        double x3=v_post*time;
        double x4=v_shock*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;
        if(x<x1)
                return rho_left;
        if(x<x2)
        {
                double c=mu*mu*(-x/time)+(1-mu*mu)*c_left;
                return rho_left*pow((c/c_left),(2/(gamma-1)));
        }
        if(x<x3)
                return rho_middle;
        if(x<x4)
                return rho_post;
        return rho_right;
}
else if(exactsolution=="noh"){
	double r=sqrt(x*x+y*y);
	if(r<=(time/3))
		return 16;
	else
		return (1+time/r);
}
else if(exactsolution=="yee"){
        double val=0;
        double r2=x*x+y*y;
        double gamma=1.4;
        double beta=5.0;
        double T=1-(gamma-1)*beta*beta/8.0/gamma/M_PI/M_PI*exp(1.0-r2);
        val=pow(T,1.0/(gamma-1.0));

        return val;
}
	return 0;
}

double VTKParticleViewer::exact_volume(double x, double y, double z, double time)
{
	if(exact_density(x,y,z,time))
		return 1.0/exact_density(x,y,z,time);
	else
		return 0.0;
}

double VTKParticleViewer::exact_pressure(double x, double y, double z, double time)
{
if(exactsolution=="simplewave"){
//	x=x-0.0025;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
        double v_post=0.92745;
        double time1=3.0+time;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);
        double x1=-c_left*time1;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time1;
        if(x<x1)
                return p_left;
        else if (x<x2)
        {
                        double c=mu*mu*(-x/time1)+(1-mu*mu)*c_left;
                        double rho=rho_left*pow((c/c_left),(2/(gamma-1)));
                        return p_left*pow(rho/rho_left,gamma);
        }
        else
        {
                        double c=mu*mu*(-x2/time1)+(1-mu*mu)*c_left;
                        double rho=rho_left*pow((c/c_left),(2/(gamma-1)));
                        return p_left*pow(rho/rho_left,gamma);
        }
}else if(exactsolution=="gresho"){
        double val;
        double r = sqrt(x*x + y*y);

        if (r <= 0.2) val = 5.0 + 12.5*r*r;
        else if  (r > 0.2 && r <= 0.4) val = 9.0 - 4*log(0.2) + 12.5*r*r - 20*r + 4*log(r);
        else val = 3.0 + 4*log(2.0);
        return val;	
}else if(exactsolution=="shocktube"){
//        x=x-0.005;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
//        double u_left=0.0;

//        double rho_right=0.125;
        double p_right=0.1;
//        double u_right=0.0;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);

//        double rho_middle=0.62847;
//        double rho_post=2.8803;
        double p_post=0.52191;
        double v_post=0.52481;
        double v_shock=0.80392;


//        double rho_middle=0.42632;
//        double rho_post=0.26557;
//        double p_post=0.30313;
//        double v_post=0.92745;
//        double v_shock=1.75216;

        double x1=-c_left*time;
        double x3=v_post*time;
        double x4=v_shock*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;

	if(x<x1)
		return p_left;
	if(x<x2)
	{
                double c=mu*mu*(-x/time)+(1-mu*mu)*c_left;
                double rho=rho_left*pow((c/c_left),(2/(gamma-1)));
                return p_left*pow(rho/rho_left,gamma);
	}
	if(x<x3)
		return p_post;
	if(x<x4)
		return p_post;
	return p_right;
}else if(exactsolution=="sod"){

//      x=x-0.005;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
        double u_left=0.0;

        double rho_right=0.125;
        double p_right=0.1;
        double u_right=0.0;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);

        double rho_middle=0.42632;
        double rho_post=0.26557;
        double p_post=0.30313;
        double v_post=0.92745;
        double v_shock=1.75216;

        double x1=-c_left*time;
        double x3=v_post*time;
        double x4=v_shock*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;

        if(x<x1)
                return p_left;
        if(x<x2)
        {
                double c=mu*mu*(-x/time)+(1-mu*mu)*c_left;
                double rho=rho_left*pow((c/c_left),(2/(gamma-1)));
                return p_left*pow(rho/rho_left,gamma);
        }
        if(x<x3)
                return p_post;
        if(x<x4)
                return p_post;
        return p_right;
}else if(exactsolution=="noh"){
        double r=sqrt(x*x+y*y);
        if(r<=(time/3))
                return 16/3;
        else
                return 1e-6;	
}else if(exactsolution=="yee"){
        double val=0;
        double r2=x*x+y*y;
        double gamma=1.4;
        double beta=5.0;
        double T=1-(gamma-1)*beta*beta/8.0/gamma/M_PI/M_PI*exp(1.0-r2);
        val=pow(T,gamma/(gamma-1.0));

        return val;
}
        return 0;
}

double VTKParticleViewer::exact_velocityu(double x, double y, double z, double time)
{
if(exactsolution=="simple"){
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
        double v_post=0.92745;
        double time1=3.0+time;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);
        double x1=-c_left*time1;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time1;
        if(x<x1)
                return 0.0;
        else if (x<x2)
        {
                        return (1-mu*mu)*(x/time1+c_left);
        }
        else
        {
                        return (1-mu*mu)*(x2/time1+c_left);
        }
}else if(exactsolution=="gresho"){
        double r = sqrt(x*x + y*y);
        double th = atan2(y,x);
        double u;
	double vX;

        if (r <= 0.2)
          {
            u = 5*r;
            vX = -u*sin(th);
//            vY = u*cos(th);
          }
        else if  (r > 0.2 && r <= 0.4)
          {
            u = 2.0 - 5*r;
            vX = -u*sin(th);
//            vY = u*cos(th);
          }
        else
          {
            vX = 0;
//            vY = 0;
	  }
	return vX;
}else if(exactsolution=="shocktube"){
//        x=x-0.005;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
        double u_left=0.0;

//        double rho_right=0.125;
//        double p_right=0.1;
        double u_right=0.0;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);

//        double rho_middle=0.62847;
//        double rho_post=2.8803;
//        double p_post=0.52191;
        double v_post=0.52481;
        double v_shock=0.80392;


//        double rho_middle=0.42632;
//        double rho_post=0.26557;
//        double p_post=0.30313;
//        double v_post=0.92745;
//        double v_shock=1.75216;

        double x1=-c_left*time;
        double x3=v_post*time;
        double x4=v_shock*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;

	if(x<x1)
		return u_left;
	if(x<x2)
		return (1-mu*mu)*(x/time+c_left);
	if(x<x3)
		return v_post;
	if(x<x4)
		return v_post;
	return u_right;
}else if(exactsolution=="sod"){

//      x=x-0.005;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
        double u_left=0.0;

        double rho_right=0.125;
        double p_right=0.1;
        double u_right=0.0;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);

        double rho_middle=0.42632;
        double rho_post=0.26557;
        double p_post=0.30313;
        double v_post=0.92745;
        double v_shock=1.75216;

        double x1=-c_left*time;
        double x3=v_post*time;
        double x4=v_shock*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;
        if(x<x1)
                return u_left;
        if(x<x2)
                return (1-mu*mu)*(x/time+c_left);
        if(x<x3)
                return v_post;
        if(x<x4)
                return v_post;
        return u_right;
}else if(exactsolution=="noh"){
        double r=sqrt(x*x+y*y);
        if(r<=(time/3))
                return 0;
        else
                return (-x/r);
}else if(exactsolution=="yee"){
        double r2=x*x+y*y;
        double beta=5.0;
        return (beta/2.0/M_PI*exp((1.0-r2)/2.0)*(-y));
}

        return 0;
}

double VTKParticleViewer::exact_velocityv(double x, double y, double z, double time)
{
if(exactsolution=="gresho"){
        double r = sqrt(x*x + y*y);
        double th = atan2(y,x);
        double u;
        double vY;

        if (r <= 0.2)
          {
            u = 5*r;
//            vX = -u*sin(th);
            vY = u*cos(th);
          }
        else if  (r > 0.2 && r <= 0.4)
          {
            u = 2.0 - 5*r;
//            vX = -u*sin(th);
            vY = u*cos(th);
          }
        else
          {
//            vX = 0;
            vY = 0;
          }
        return vY;
}else if(exactsolution=="noh"){
        double r=sqrt(x*x+y*y);
        if(r<=(time/3))
                return 0;
        else
                return (-y/r);
}else if(exactsolution=="yee"){
        double r2=x*x+y*y;
        double beta=5.0;
        return (beta/2.0/M_PI*exp((1.0-r2)/2.0)*x);
}
	return 0;
}

double VTKParticleViewer::exact_velocity_magnitude(double x, double y, double z, double time)
{
	return exact_velocityu(x,y,z,time)*exact_velocityu(x,y,z,time)+exact_velocityv(x,y,z,time)*exact_velocityv(x,y,z,time);
}



////////////////////////////////////////////////////////////////////////////////////////
// Start of VTKParticleViewer
////////////////////////////////////////////////////////////////////////////////////////

VTKParticleViewer::VTKParticleViewer(ParticleData* data, const std::string& particleType, 
const string& outputfileName, int numDigits) {
	m_pParticleData = data;
	m_sOutputfileName = outputfileName;
	m_iNumDigits = numDigits; 
	m_sParticleType = particleType;
	outputerror = false;
	exactsolution = "";
	auxiliaryoutput = "";
}


int VTKParticleViewer::writeResult(double time, size_t writeStep) {
	 
	// alias pointers
	const double* positionX = m_pParticleData->getPositionX();
	const double* positionY = m_pParticleData->getPositionY();
	const double* positionZ = m_pParticleData->getPositionZ();
	if(m_pParticleData->getDimension()==2)
		positionZ=positionY;
	const double* velocityU = m_pParticleData->getVelocityU();
	const double* velocityV = m_pParticleData->getVelocityV();
	const double* velocityW = m_pParticleData->getVelocityW();
    const double* temperature = m_pParticleData->getTemperature();
	const double* volume = m_pParticleData->getVolume();
	const double* mass = m_pParticleData->getMass();
	const double* pressure = m_pParticleData->getPressure();
	const double* soundSpeed = m_pParticleData->getSoundSpeed();
        const double* voronoi_volume = m_pParticleData->getVolumeVoronoi();
	const double* leftintegral = m_pParticleData->getLeftIntegral();
	const double* rightintegral = m_pParticleData->getRightIntegral();
	const double* Deltaq = m_pParticleData->getDeltaQ();
	const double* Qplusminus = m_pParticleData->getQplusminus();
    
    int* timetrack = m_pParticleData->getTimeTrack();

#ifdef LW_DEBUG        
        
//	const double* phi = m_pParticleData->getPhi();
        const double* perror0 = m_pParticleData->getPError0();
        const double* perror1 = m_pParticleData->getPError1();
        const double* velerror0 = m_pParticleData->getVelError0();
        const double* velerror1 = m_pParticleData->getVelError1();        
	    const int* IfSPHDensity = m_pParticleData->getIfSPHDensity();
#endif
#ifdef DEBUG_LW_S
        const double* pxl = m_pParticleData->getPxl();
        const double* pxr = m_pParticleData->getPxr();
        const double* vxl = m_pParticleData->getVxl();
        const double* vxr = m_pParticleData->getVxr();
        const double* pyl = m_pParticleData->getPyl();
        const double* pyr = m_pParticleData->getPyr();
        const double* vyl = m_pParticleData->getVyl();
        const double* vyr = m_pParticleData->getVyr();
        const double* vtx = m_pParticleData->getVtx();
        const double* vty = m_pParticleData->getVty();
        const double* ptx = m_pParticleData->getPtx();
        const double* pty = m_pParticleData->getPty();
        const double* volumetx = m_pParticleData->getVolumetx();
        const double* volumety = m_pParticleData->getVolumety();
#endif
#ifdef DEBUG_SWITCH
	const double* volume_x = m_pParticleData->getVolume_x();
        const double* volume_y = m_pParticleData->getVolume_y();
        const double* volume_z = m_pParticleData->getVolume_z();
#endif
	const double* localParSpacing = m_pParticleData->getLocalParSpacing();

	const int* objectTag = m_pParticleData->getObjectTag();

	// Create an output file the name "filename"
	string filename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".vtk";
	FILE *outfile;
	outfile = fopen(filename.c_str(), "w");
	if(outfile==nullptr) {
		printf("Unable to open file: %s\n",filename.c_str()); 
		return 1;
	}
	size_t startIndex, numParticle;
	if(m_sParticleType=="all") {
		//startIndex = 0;
		startIndex = m_pParticleData->getFluidStartIndex()+m_pParticleData->getFluidNum();
		numParticle = m_pParticleData->getTotalNum()-m_pParticleData->getFluidNum();
	}
	else if(m_sParticleType=="fluid") {
		startIndex = m_pParticleData->getFluidStartIndex();
		numParticle = m_pParticleData->getFluidNum();
	}

	size_t endIndex = startIndex + numParticle;


	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.16g\n",time);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	
	fprintf(outfile,"POINTS %ld double\n",numParticle);

    if(m_pParticleData->getDimension()==2) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",positionX[i], positionY[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",positionX[i], positionY[i], positionZ[i]);
	}

    fprintf(outfile,"POINT_DATA %ld\n",numParticle);
if(m_pParticleData->m_iPrintVelocity){
	
	fprintf(outfile,"VECTORS Velocity double\n");
	if(m_pParticleData->getDimension()==2) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",velocityU[i], velocityV[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",velocityU[i], velocityV[i], velocityW[i]);
	}

}
if(m_pParticleData->m_iPrintPressure){
	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",pressure[i]);
}	

if(m_pParticleData->m_iPrintVolume){
	fprintf(outfile,"SCALARS volume double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",volume[i]);
}
/*        fprintf(outfile,"SCALARS volume_voronoi double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",voronoi_volume[i]);
*/

if(m_pParticleData->m_iPrintDensity){
        fprintf(outfile,"SCALARS density double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
	{
		if(fabs(volume[i])>1e-10)
                	fprintf(outfile,"%.16g\n",1.0/volume[i]);
		else
			fprintf(outfile,"%.16g\n",0.0);
	}	
}
if(m_pParticleData->m_iPrintSoundSpeed){
	fprintf(outfile,"SCALARS sound_speed double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",soundSpeed[i]);
}

if(m_pParticleData->m_iPrintMass){
        fprintf(outfile,"SCALARS mass double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",mass[i]);
}
/*	fprintf(outfile,"SCALARS object_tag int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",objectTag[i]);
*/
if(m_pParticleData->m_iPrintLocalSpacing){    
	fprintf(outfile,"SCALARS local_par_spacing double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",localParSpacing[i]);
}

if(m_pParticleData->m_iPrintVelocityU){
	fprintf(outfile,"SCALARS velocity_u double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",velocityU[i]);
	}
if(m_pParticleData->m_iPrintVelocityV){
	fprintf(outfile,"SCALARS velocity_v double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",velocityV[i]);
	}
	if(m_pParticleData->getDimension()==3) {
        if(m_pParticleData->m_iPrintVelocityW){
		fprintf(outfile,"SCALARS velocity_w double\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g\n",velocityW[i]);
    }
}	
if(m_pParticleData->m_iPrintTemperature){
	fprintf(outfile,"SCALARS temperature double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",temperature[i]);
	}


/*	fprintf(outfile,"SCALARS index int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%zu\n",i);
       
        fprintf(outfile,"SCALARS order double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",phi[i]);
        fprintf(outfile,"SCALARS density_modification_count int\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%zu\n",IfSPHDensity[i]);
*/
if(m_pParticleData->getNumberofPellet()){
	if(m_pParticleData->getDimension()==3){
        if(m_pParticleData->m_iPrintLeftIntegral){
            fprintf(outfile,"SCALARS left_integral double\n");
        	fprintf(outfile,"LOOKUP_TABLE default\n");
	        for(size_t i=startIndex; i<endIndex; i++)
        	        fprintf(outfile,"%.16g\n",leftintegral[i]);
        }
        if(m_pParticleData->m_iPrintRightIntegral){
	        fprintf(outfile,"SCALARS right_integral double\n");
        	fprintf(outfile,"LOOKUP_TABLE default\n");
	        for(size_t i=startIndex; i<endIndex; i++)
        	        fprintf(outfile,"%.16g\n",rightintegral[i]);
        }
        if(m_pParticleData->m_iPrintDeltaq){
                fprintf(outfile,"SCALARS delta_q double\n");
                fprintf(outfile,"LOOKUP_TABLE default\n");
                for(size_t i=startIndex; i<endIndex; i++)
                        fprintf(outfile,"%.16g\n",Deltaq[i]);
        }
        if(m_pParticleData->m_iPrintQplusminus){
                fprintf(outfile,"SCALARS q_plusminus double\n");
                fprintf(outfile,"LOOKUP_TABLE default\n");
                for(size_t i=startIndex; i<endIndex; i++)
                        fprintf(outfile,"%.16g\n",Qplusminus[i]);
	    }
    }
}
        double Delta_T = 0.0005;

        int timelabel =  (int)(time/Delta_T);
        if(m_pParticleData->m_iPrintTimeTrack){
                fprintf(outfile,"SCALARS timetrack int\n");
                fprintf(outfile,"LOOKUP_TABLE default\n");
                for(size_t i=startIndex; i<endIndex; i++){
                    double t = temperature[i];
                    if (t>1&&timetrack[i]==0){
                        timetrack[i] = timelabel;
                        fprintf(outfile,"%d\n",timelabel);

                    }
                    else
                        fprintf(outfile,"%d\n",0);
              }
        }
#ifdef DEBUG_LW_V
//        fprintf(outfile,"SCALARS phi double\n");
        fprintf(outfile,"SCALARS u_x double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",phi[i]);

//        fprintf(outfile,"SCALARS perror0 double\n");
        fprintf(outfile,"SCALARS u_y double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",perror0[i]);

//        fprintf(outfile,"SCALARS perror1 double\n");
        fprintf(outfile,"SCALARS u_xx double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",perror1[i]);

//        fprintf(outfile,"SCALARS velerror0 double\n");
        fprintf(outfile,"SCALARS u_tt double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",velerror0[i]);

//        fprintf(outfile,"SCALARS velerror1 double\n");
        fprintf(outfile,"SCALARS p_tt double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",velerror1[i]);


#endif
#ifdef DEBUG_SWITCH
        fprintf(outfile,"SCALARS volume_x double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",volume_x[i]);
        fprintf(outfile,"SCALARS volume_y double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",volume_y[i]);
        fprintf(outfile,"SCALARS volume_z double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",volume_z[i]);

        fprintf(outfile,"VECTORS Volume_gradient double\n");
        if(m_pParticleData->getDimension()==2) {
                for(size_t i=startIndex; i<endIndex; i++)
                        fprintf(outfile,"%.16g %.16g %.16g\n",volume_x[i], volume_y[i], 0.0);
        }
        else if(m_pParticleData->getDimension()==3) {
                for(size_t i=startIndex; i<endIndex; i++)
                        fprintf(outfile,"%.16g %.16g %.16g\n",volume_x[i], volume_y[i], volume_z[i]);
        }
#endif
#ifdef DEBUG_LW_S
        fprintf(outfile,"SCALARS p_xl double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",pxl[i]);
        fprintf(outfile,"SCALARS p_xr double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",pxr[i]);
        fprintf(outfile,"SCALARS v_xl double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",vxl[i]);
        fprintf(outfile,"SCALARS v_xr double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",vxr[i]);
        fprintf(outfile,"SCALARS p_yl double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",pyl[i]);
        fprintf(outfile,"SCALARS p_yr double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",pyr[i]);
        fprintf(outfile,"SCALARS v_yl double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",vyl[i]);
        fprintf(outfile,"SCALARS v_yr double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",vyr[i]);
        fprintf(outfile,"SCALARS p_tx double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",ptx[i]);
        fprintf(outfile,"SCALARS p_ty double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",pty[i]);
        fprintf(outfile,"SCALARS v_tx double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",vtx[i]);
        fprintf(outfile,"SCALARS v_ty double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",vty[i]);
        fprintf(outfile,"SCALARS volume_tx double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",volumetx[i]);
        fprintf(outfile,"SCALARS volume_ty double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",volumety[i]);
#endif
if(outputerror){
        fprintf(outfile,"SCALARS exact_density double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                fprintf(outfile,"%.16g\n",exact_density(positionX[i],positionY[i],positionZ[i],time));
        }

        double l_1_error_density=0;
        double l_inf_error_density=0;
        fprintf(outfile,"SCALARS error_density double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
		if(fabs(volume[i])>1e-10)
		{
			l_1_error_density=l_1_error_density+fabs(1.0/volume[i]-exact_density(positionX[i],positionY[i],positionZ[i],time));
			if (fabs(1.0/volume[i]-exact_density(positionX[i],positionY[i],positionZ[i],time))>l_inf_error_density)
				l_inf_error_density=fabs(1.0/volume[i]-exact_density(positionX[i],positionY[i],positionZ[i],time));
	                fprintf(outfile,"%.16g\n",fabs(1.0/volume[i]-exact_density(positionX[i],positionY[i],positionZ[i],time)));
		}
		else
			fprintf(outfile,"%.16g\n",0.0);
        }
        if(endIndex>startIndex)
                l_1_error_density=l_1_error_density/(endIndex-startIndex);
        printf("L1 error of density: %.2e\n",l_1_error_density);
        printf("L_inf error of density: %.2e\n",l_inf_error_density);


        fprintf(outfile,"SCALARS exact_pressure double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                fprintf(outfile,"%.16g\n",exact_pressure(positionX[i],positionY[i],positionZ[i],time));
        }

	double l_1_error_pressure=0;
	double l_inf_error_pressure=0;
        fprintf(outfile,"SCALARS error_pressure double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
		l_1_error_pressure=l_1_error_pressure+fabs(pressure[i]-exact_pressure(positionX[i],positionY[i],positionZ[i],time));
		if (fabs(pressure[i]-exact_pressure(positionX[i],positionY[i],positionZ[i],time))>l_inf_error_pressure)
			l_inf_error_pressure=fabs(pressure[i]-exact_pressure(positionX[i],positionY[i],positionZ[i],time));
                fprintf(outfile,"%.16g\n",fabs(pressure[i]-exact_pressure(positionX[i],positionY[i],positionZ[i],time)));
        }
	if(endIndex>startIndex)
		l_1_error_pressure=l_1_error_pressure/(endIndex-startIndex);
	printf("L1 error of pressure: %.2e\n",l_1_error_pressure);
        printf("L_inf error of pressure: %.2e\n",l_inf_error_pressure);

        fprintf(outfile,"SCALARS error_pressure_raw double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                fprintf(outfile,"%.16g\n",pressure[i]-exact_pressure(positionX[i],positionY[i],positionZ[i],time));
        }


        fprintf(outfile,"SCALARS exact_velocity_u double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                fprintf(outfile,"%.16g\n",exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
        }

        fprintf(outfile,"SCALARS error_velocity_u double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                fprintf(outfile,"%.16g\n",fabs(velocityU[i]-exact_velocityu(positionX[i],positionY[i],positionZ[i],time)));
        }

        double l_1_error_velocity=0;
        double l_inf_error_velocity=0;
        fprintf(outfile,"SCALARS error_velocity_magnitude double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
		l_1_error_velocity=l_1_error_velocity+fabs(velocityU[i]*velocityU[i]+velocityV[i]*velocityV[i]-exact_velocity_magnitude(positionX[i],positionY[i],positionZ[i],time));
		if(fabs(velocityU[i]*velocityU[i]+velocityV[i]*velocityV[i]-exact_velocity_magnitude(positionX[i],positionY[i],positionZ[i],time))>l_inf_error_velocity)
			l_inf_error_velocity=fabs(velocityU[i]*velocityU[i]+velocityV[i]*velocityV[i]-exact_velocity_magnitude(positionX[i],positionY[i],positionZ[i],time));
                fprintf(outfile,"%.16g\n",fabs(velocityU[i]*velocityU[i]+velocityV[i]*velocityV[i]-exact_velocity_magnitude(positionX[i],positionY[i],positionZ[i],time)));
        }
        if(endIndex>startIndex)
                l_1_error_velocity=l_1_error_velocity/(endIndex-startIndex);
        printf("L1 error of velocity magnitude: %.2e\n",l_1_error_velocity);
        printf("L_inf error of velocity magnitude: %.2e\n",l_inf_error_velocity);


#ifdef DEBUG_LW_V
        fprintf(outfile,"SCALARS exact_u_x double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",(exact_volume(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_volume(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2);
	}

        fprintf(outfile,"SCALARS error_u_x double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",fabs(phi[i]-(exact_volume(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_volume(positionX[i],positionY[i],positionZ[i]-dx,time))/dx/2));
        }

        fprintf(outfile,"SCALARS exact_u_xx double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",(exact_volume(positionX[i]+dx,positionY[i],positionZ[i],time)-2*exact_volume(positionX[i],positionY[i],positionZ[i],time)+exact_volume(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/dx);
        }

        fprintf(outfile,"SCALARS error_u_xx double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",fabs(perror1[i]-(exact_volume(positionX[i]+dx,positionY[i],positionZ[i],time)-2*exact_volume(positionX[i],positionY[i],positionZ[i],time)+exact_volume(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/dx));
	}

        fprintf(outfile,"SCALARS exact_u_tt double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
		double dx=0.00001;
		double dt=0.00001;
		double ftt=(exact_velocityu(positionX[i],positionY[i],positionZ[i],time+dt)-2*exact_velocityu(positionX[i],positionY[i],positionZ[i],time)+exact_velocityu(positionX[i],positionY[i],positionZ[i],time-dt))/dt/dt;
		double fxt=(exact_velocityu(positionX[i]+dx,positionY[i],positionZ[i],time+dt)+exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time-dt)-exact_velocityu(positionX[i],positionY[i],positionZ[i]+dx,time-dt)-exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time+dt))/dx/dt/4;
		double ut=(exact_velocityu(positionX[i],positionY[i],positionZ[i],time+dt)-exact_velocityu(positionX[i],positionY[i],positionZ[i],time-dt))/dt/2;
		double fx=(exact_velocityu(positionX[i],positionY[i],positionZ[i]+dx,time)-exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2;
                double fxx=(exact_velocityu(positionX[i]+dx,positionY[i],positionZ[i],time)-2*exact_velocityu(positionX[i],positionY[i],positionZ[i],time)+exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/dx;
		double ux=(exact_velocityu(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2;
		double f_tt=ftt+2*exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*fxt+ut*fx+exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*ux*fx+exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*fxx;
                fprintf(outfile,"%.16g\n",f_tt);

 
	}

        fprintf(outfile,"SCALARS exact_p_tt double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                double dt=0.00001;
                double ftt=(exact_pressure(positionX[i],positionY[i],positionZ[i],time+dt)-2*exact_pressure(positionX[i],positionY[i],positionZ[i],time)+exact_pressure(positionX[i],positionY[i],positionZ[i],time-dt))/dt/dt;
                double fxt=(exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time+dt)+exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time-dt)-exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time-dt)-exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time+dt))/dx/dt/4;
                double ut=(exact_velocityu(positionX[i],positionY[i],positionZ[i],time+dt)-exact_velocityu(positionX[i],positionY[i],positionZ[i],time-dt))/dt/2;
                double fx=(exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2;
                double fxx=(exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time)-2*exact_pressure(positionX[i],positionY[i],positionZ[i],time)+exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/dx;
                double ux=(exact_velocityu(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2;
                double f_tt=ftt+2*exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*fxt+ut*fx+exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*ux*fx+exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*fxx;
                fprintf(outfile,"%.16g\n",f_tt);

        }

        fprintf(outfile,"SCALARS error_p_tt double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                double dt=0.00001;
                double ftt=(exact_pressure(positionX[i],positionY[i],positionZ[i],time+dt)-2*exact_pressure(positionX[i],positionY[i],positionZ[i],time)+exact_pressure(positionX[i],positionY[i],positionZ[i],time-dt))/dt/dt;
                double fxt=(exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time+dt)+exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time-dt)-exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time-dt)-exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time+dt))/dx/dt/4;
                double ut=(exact_velocityu(positionX[i],positionY[i],positionZ[i],time+dt)-exact_velocityu(positionX[i],positionY[i],positionZ[i],time-dt))/dt/2;
                double fx=(exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2;
                double fxx=(exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time)-2*exact_pressure(positionX[i],positionY[i],positionZ[i],time)+exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/dx;
                double ux=(exact_velocityu(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2;
                double f_tt=ftt+2*exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*fxt+ut*fx+exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*ux*fx+exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*exact_velocityu(positionX[i],positionY[i],positionZ[i],time)*fxx;
                fprintf(outfile,"%.16g\n",fabs(velerror1[i]-f_tt));
	}
#endif
#ifdef DEBUG_LW_S
        fprintf(outfile,"SCALARS exact_p_x double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",(exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2);
        }

        fprintf(outfile,"SCALARS error_p_xl double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",fabs(pxl[i]-(exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2));
        }


        fprintf(outfile,"SCALARS error_p_xr double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",fabs(pxr[i]-(exact_pressure(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_pressure(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2));
        }

        fprintf(outfile,"SCALARS exact_u_x double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",(exact_velocityu(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2);
        }

        fprintf(outfile,"SCALARS error_u_xl double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",fabs(vxl[i]-(exact_velocityu(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2));
        }


        fprintf(outfile,"SCALARS error_u_xr double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dx=0.00001;
                fprintf(outfile,"%.16g\n",fabs(vxr[i]-(exact_velocityu(positionX[i]+dx,positionY[i],positionZ[i],time)-exact_velocityu(positionX[i]-dx,positionY[i],positionZ[i],time))/dx/2));
        }

        fprintf(outfile,"SCALARS exact_u_t double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dt=0.00001;
		double u1=exact_velocityu(positionX[i],positionY[i],positionZ[i],time);
		double v1=exact_velocityv(positionX[i],positionY[i],positionZ[i],time);

                double ut=(exact_velocityu(positionX[i]+u1*dt,positionY[i]+v1*dt,positionZ[i],time+dt)-exact_velocityu(positionX[i]-u1*dt,positionY[i]-v1*dt,positionZ[i],time-dt))/dt/2;
                fprintf(outfile,"%.16g\n",ut);
        }

        fprintf(outfile,"SCALARS error_utx double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dt=0.00001;
                double u1=exact_velocityu(positionX[i],positionY[i],positionZ[i],time);
                double v1=exact_velocityv(positionX[i],positionY[i],positionZ[i],time);

                double ut=(exact_velocityu(positionX[i]+u1*dt,positionY[i]+v1*dt,positionZ[i],time+dt)-exact_velocityu(positionX[i]-u1*dt,positionY[i]-v1*dt,positionZ[i],time-dt))/dt/2;
                fprintf(outfile,"%.16g\n",fabs(vtx[i]-ut));
        }

        fprintf(outfile,"SCALARS exact_p_t double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dt=0.00001;
                double u1=exact_velocityu(positionX[i],positionY[i],positionZ[i],time);
                double v1=exact_velocityv(positionX[i],positionY[i],positionZ[i],time);

                double pt=(exact_pressure(positionX[i]+u1*dt,positionY[i]+v1*dt,positionZ[i],time+dt)-exact_pressure(positionX[i]-u1*dt,positionY[i]-v1*dt,positionZ[i],time-dt))/dt/2;
                fprintf(outfile,"%.16g\n",pt);
        }

        fprintf(outfile,"SCALARS error_p_t double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dt=0.00001;
                double u1=exact_velocityu(positionX[i],positionY[i],positionZ[i],time);
                double v1=exact_velocityv(positionX[i],positionY[i],positionZ[i],time);

                double pt=(exact_pressure(positionX[i]+u1*dt,positionY[i]+v1*dt,positionZ[i],time+dt)-exact_pressure(positionX[i]-u1*dt,positionY[i]-v1*dt,positionZ[i],time-dt))/dt/2;
                fprintf(outfile,"%.16g\n",fabs(ptx[i]+pty[i]-pt));
        }

        fprintf(outfile,"SCALARS exact_volume_t double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dt=0.00001;
                double u1=exact_velocityu(positionX[i],positionY[i],positionZ[i],time);
                double v1=exact_velocityv(positionX[i],positionY[i],positionZ[i],time);

                double volumet=(exact_volume(positionX[i]+u1*dt,positionY[i]+v1*dt,positionZ[i],time+dt)-exact_volume(positionX[i]-u1*dt,positionY[i]-v1*dt,positionZ[i],time-dt))/dt/2;
                fprintf(outfile,"%.16g\n",volumet);
        }

        fprintf(outfile,"SCALARS error_volume_t double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
        {
                double dt=0.00001;
                double u1=exact_velocityu(positionX[i],positionY[i],positionZ[i],time);
                double v1=exact_velocityv(positionX[i],positionY[i],positionZ[i],time);

                double volumet=(exact_volume(positionX[i]+u1*dt,positionY[i]+v1*dt,positionZ[i],time+dt)-exact_volume(positionX[i]-u1*dt,positionY[i]-v1*dt,positionZ[i],time-dt))/dt/2;
                fprintf(outfile,"%.16g\n",fabs(volumetx[i]+volumety[i]-volumet));
        }

#endif
}
	fclose(outfile);

	if(auxiliaryoutput=="onedplot"){
        string onedfilename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".txt";
        FILE *onedoutfile;
        onedoutfile = fopen(onedfilename.c_str(), "w");
        if(onedoutfile==nullptr) {
                printf("Unable to open file: %s\n",onedfilename.c_str());
                return 1;
        }
        size_t i0=startIndex;
        while(positionY[i0]<0)  i0++;
        while(positionX[i0]>positionX[i0-1]) i0++;
        size_t i=i0;
        while(positionX[i]<positionX[i+1])
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
	i0=i;
        while(positionY[i0]>-0.1)  i0++;
        while(positionY[i0]<0)  i0++;
        while(positionX[i0]>positionX[i0-1]) i0++;
        i=i0;
        while(positionX[i]<positionX[i+1])
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        fclose(onedoutfile);
	}

	if(auxiliaryoutput=="onedplottwoline"){
        string onedfilename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".txt";
        FILE *onedoutfile;
        onedoutfile = fopen(onedfilename.c_str(), "w");
        if(onedoutfile==nullptr) {
                printf("Unable to open file: %s\n",onedfilename.c_str());
                return 1;
        }
        size_t i0=startIndex;
        while(positionY[i0]<0)  i0++;
        while(positionX[i0]<-0.1) i0++;
        while(positionX[i0]>=-0.1) i0++;
        size_t i=i0;
        while(positionX[i]<-0.1)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionX[i]>=-0.1)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionX[i]<-0.1)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionX[i]>=-0.1)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }

        fclose(onedoutfile);
	}

	if(auxiliaryoutput=="onedplotfourline"){
        string onedfilename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".txt";
        FILE *onedoutfile;
        onedoutfile = fopen(onedfilename.c_str(), "w");
        if(onedoutfile==nullptr) {
                printf("Unable to open file: %s\n",onedfilename.c_str());
                return 1;
        }
        size_t i0=startIndex;
        while((positionY[i0]<0) || (positionZ[i0]<0))  i0++;
        while(positionY[i0]>=0) i0++;
        while(positionY[i0]<0) i0++;
        while(positionX[i0]<0) i0++;
        while(positionX[i0]>=0) i0++;
        size_t i=i0;
        while(positionX[i]<0)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionX[i]>=0)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionX[i]<0)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionX[i]>=0)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionY[i0]>=0) i0++;
        while(positionY[i0]<0) i0++;
        while(positionX[i0]<0) i0++;
        while(positionX[i0]>=0) i0++;
        i=i0;
        while(positionX[i]<0)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionX[i]>=0)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionX[i]<0)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        while(positionX[i]>=0)
        {
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
                i=i+1;
        }
        fclose(onedoutfile);
	}

	if(auxiliaryoutput=="nozzle"){
        string onedfilename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".txt";
        FILE *onedoutfile;
        onedoutfile = fopen(onedfilename.c_str(), "w");
        if(onedoutfile==nullptr) {
                printf("Unable to open file: %s\n",onedfilename.c_str());
                return 1;
        }
	for(size_t i0=startIndex;i0<endIndex;i0++)
	{
		double temp=positionY[i0]*positionY[i0];
		if(m_pParticleData->getDimension()==3)
			temp+=positionZ[i0]*positionZ[i0];
		temp=sqrt(temp);
//		if(temp<positionX[i0]*0.02)
		if((positionX[i0]-0.0079)*(positionX[i0]-0.008)<0)
			fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g\n",temp, 1.0/volume[i0],  pressure[i0], velocityU[i0], mass[i0], positionX[i0]);
	}
        fclose(onedoutfile);


        string onedfilename2 = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + "_begin_.txt";
        FILE *onedoutfile2;
        onedoutfile2 = fopen(onedfilename2.c_str(), "w");
        if(onedoutfile2==nullptr) {
                printf("Unable to open file: %s\n",onedfilename2.c_str());
                return 1;
        }
        for(size_t i0=startIndex;i0<endIndex;i0++)
        {
                double temp=positionY[i0]*positionY[i0];
                if(m_pParticleData->getDimension()==3)
                        temp+=positionZ[i0]*positionZ[i0];
                temp=sqrt(temp);
                if((positionX[i0]-0)*(positionX[i0]-0.00001)<0)
                        fprintf(onedoutfile2,"%.16g %.16g %.16g %.16g %.16g %.16g\n",temp, 1.0/volume[i0],  pressure[i0], velocityU[i0], mass[i0], positionX[i0]);
        }
        fclose(onedoutfile2);
	}

	if(auxiliaryoutput=="region"){
        string onedfilename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".txt";
        FILE *onedoutfile;
        onedoutfile = fopen(onedfilename.c_str(), "w");
        if(onedoutfile==nullptr) {
                printf("Unable to open file: %s\n",onedfilename.c_str());
                return 1;
        }
        for(size_t i=startIndex;i<endIndex;i++)
        {
                if((positionX[i]+2)*(positionX[i]-2)<0&&(positionY[i]+0.1)*(positionY[i]-0.1)<0&&(positionZ[i]+0.1)*(positionZ[i]-0.1)<0)
                fprintf(onedoutfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",positionX[i], 1.0/volume[i], exact_density(positionX[i],positionY[i],positionZ[i],time),pressure[i], exact_pressure(positionX[i],positionY[i],positionZ[i],time), velocityU[i], exact_velocityu(positionX[i],positionY[i],positionZ[i],time));
        }
        fclose(onedoutfile);
	}

	if(auxiliaryoutput=="radial"){
        string radialfilename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + "_radial_plot.dat";
        FILE *radialoutfile;
        radialoutfile = fopen(radialfilename.c_str(), "w");
        if(radialoutfile==nullptr) {
                printf("Unable to open file: %s\n",radialfilename.c_str());
                return 1;
        }
        for(size_t i=startIndex; i<endIndex; i++)
        {
//		if((positionX[i]*positionX[i]+positionY[i]*positionY[i])<0.36)
			fprintf(radialoutfile,"%.16g %.16g %.16g %.16g\n",sqrt(positionX[i]*positionX[i]+positionY[i]*positionY[i]),sqrt(velocityU[i]*velocityU[i]+velocityV[i]*velocityV[i]),pressure[i],1.0/volume[i]);
        }
        fclose(radialoutfile);
	}

	if(auxiliaryoutput=="kelvinhelmholtz"){
	double max_energy=0;
	double big_energy[10];
	for(size_t j=0;j<10;j++) big_energy[j]=0;
	double s_sum=0;
	double c_sum=0;
	double d_sum=0;
	double mode=0;

        for(size_t i=startIndex; i<endIndex; i++)
        {	
		if(volume[i]>0)
		{
//			if ((0.5*velocityV[i]*velocityV[i]/volume[i])>max_energy)
//				max_energy=0.5*velocityV[i]*velocityV[i]/volume[i];
			double energy_y=0.5*velocityV[i]*velocityV[i]/volume[i];
			for(size_t j=0;j<10;j++)
			{
				if(energy_y>big_energy[j])
				{
					for(size_t k=9;k>j;k--)
						big_energy[k]=big_energy[k-1];
					big_energy[j]=energy_y;
					break;
				}
			}

		}		
		if(positionY[i]<0.5)
		{
			s_sum=s_sum+velocityV[i]*volume[i]*sin(4.0*M_PI*positionX[i])*exp(-4.0*M_PI*fabs(positionY[i]-0.25));
			c_sum=c_sum+velocityV[i]*volume[i]*cos(4.0*M_PI*positionX[i])*exp(-4.0*M_PI*fabs(positionY[i]-0.25));
			d_sum=d_sum+volume[i]*exp(-4.0*M_PI*fabs(positionY[i]-0.25));
		}
		else
		{
			s_sum=s_sum+velocityV[i]*volume[i]*sin(4.0*M_PI*positionX[i])*exp(-4.0*M_PI*fabs(1.0-positionY[i]-0.25));
			c_sum=c_sum+velocityV[i]*volume[i]*cos(4.0*M_PI*positionX[i])*exp(-4.0*M_PI*fabs(1.0-positionY[i]-0.25));
			d_sum=d_sum+volume[i]*exp(-4.0*M_PI*fabs(1.0-positionY[i]-0.25));
		}
        }

	max_energy=big_energy[9];
	mode=2.0*sqrt((s_sum*s_sum+c_sum*c_sum)/d_sum/d_sum);
        string growthfilename = m_sOutputfileName + "_growth.txt";
        FILE *growthoutfile;
	if(writeStep>0)
	        growthoutfile = fopen(growthfilename.c_str(), "a");
	else
                growthoutfile = fopen(growthfilename.c_str(), "w");
        if(growthoutfile==nullptr) {
                printf("Unable to open file: %s\n",growthfilename.c_str());
                return 1;
        }
	fprintf(growthoutfile,"%.16g %.16g %.16g\n",time,max_energy,mode);
        fclose(growthoutfile);
	}

	if(auxiliaryoutput=="tvplot"){
        string radialfilename = m_sOutputfileName + "_tvplot.txt";
        FILE *radialoutfile;
	if(writeStep==0)
		radialoutfile = fopen(radialfilename.c_str(), "w");
	else
        	radialoutfile = fopen(radialfilename.c_str(), "a");
        if(radialoutfile==nullptr) {
                printf("Unable to open file: %s\n",radialfilename.c_str());
                return 1;
        }
#ifdef LW_DEBUG
    double tv_p=0,tv_rho=0,tv_u=0;
        for(size_t i=startIndex; i<endIndex; i++)
        {
		tv_p+=perror0[i]*volume[i];
		tv_rho+=perror1[i]*volume[i];
		tv_u+=velerror0[i]*volume[i];
        }
	fprintf(radialoutfile,"%.16g %.16g %.16g %.16g\n",time,tv_p/(endIndex-startIndex),tv_rho/(endIndex-startIndex),tv_u/(endIndex-startIndex));	
        fclose(radialoutfile);
	
    
#endif
    }
	if(auxiliaryoutput=="pelletablation"){
        
       /* double E = 0;
        for(size_t i = startIndex; i<endIndex;i++){
           double e = m_fDt*mass[i]*volume[i]*Deltaq[i];
           E += e;
        }*/
    	string coupledquantity  = m_sOutputfileName + "_test_" + rightFlush(writeStep, m_iNumDigits) + ".txt";
	    FILE *outfile;
	    outfile = fopen(coupledquantity.c_str(), "w");
	    if(outfile==nullptr) {
		printf("Unable to open file: %s\n",filename.c_str()); 
		return 1;
	    }
        for(size_t i = startIndex; i<endIndex;i++){
           
           double x = positionX[i];
           double y = positionY[i];
           double z = positionZ[i];
           
           double vx = velocityU[i];
           double vy = velocityV[i];
           double vz = velocityW[i];
           

           double r_cord = sqrt(y*y+z*z);
           double phi = atan(z/y);
           double v_phi = vy*(-z)/r_cord + vz*y/r_cord;
           double v_r = vy*y/r_cord + vz*z/r_cord;
           double T = temperature[i];
           double V = volume[i];
           double P = pressure[i];
          
           fprintf(outfile,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",r_cord,phi,x,v_r,v_phi,vx,T,V,P,1./V);
        
        } 


        string mfrfilename = m_sOutputfileName + "_massflowrate.txt";
        FILE *mfroutfile;
        if(writeStep>0)
                mfroutfile = fopen(mfrfilename.c_str(), "a");
        else{
                mfroutfile = fopen(mfrfilename.c_str(), "w");
		fprintf(mfroutfile,"Time  Massflowrate around pellets.\n");
	}
        if(mfroutfile==nullptr) {
                printf("Unable to open file: %s\n",mfrfilename.c_str());
                return 1;
        }

    double *massflowrate = m_pParticleData->getMassFlowRate();


	double r,dis,dx,dr,mfr;
	int NumberofPellet = m_pParticleData->getNumberofPellet();
	const double* px = m_pParticleData->getPelletPositionX();
	const double* py = m_pParticleData->getPelletPositionY();
	const double* pz = m_pParticleData->getPelletPositionZ();

	fprintf(mfroutfile,"%.16g ",time);

	r=1;
	dis=10;

	for(int k=0;k<NumberofPellet;k++)
	{
		
	    fprintf(mfroutfile,"%.16g ",massflowrate[k]);
        mfr = 0;

		for(size_t i=startIndex;i<endIndex;i++)
		{
			double tr=(positionX[i]-px[k])*(positionX[i]-px[k])+(positionY[i]-py[k])*(positionY[i]-py[k])+(positionZ[i]-pz[k])*(positionZ[i]-pz[k]);
			tr=fabs(sqrt(tr)-r);
			if(tr<dis)
			{	
				dis=tr;
				dx=localParSpacing[i];
			}
		}
		dr=5*dx;
		for(size_t i=startIndex;i<endIndex;i++)
		{	
        	double tr=(positionX[i]-px[k])*(positionX[i]-px[k])+(positionY[i]-py[k])*(positionY[i]-py[k])+(positionZ[i]-pz[k])*(positionZ[i]-pz[k]);
        	tr=fabs(sqrt(tr)-r);
			if(tr<dr)
			{
				double vr=velocityU[i]*(positionX[i]-px[k])+velocityV[i]*(positionY[i]-py[k])+velocityW[i]*(positionZ[i]-pz[k]);
				vr=vr/sqrt((positionX[i]-px[k])*(positionX[i]-px[k])+(positionY[i]-py[k])*(positionY[i]-py[k])+(positionZ[i]-pz[k])*(positionZ[i]-pz[k]));
				mfr+=mass[i]*vr;

			}
		}
		mfr=mfr/2/dr;
		fprintf(mfroutfile,"%.16g ", mfr);
		if((k+1)%NumberofPellet == 0)
			fprintf(mfroutfile,"\n");
	}
	


/*
	double mfr_1=mfr;
    
        r = 1;
        dis=10;
        mfr=0;
        for(size_t i=startIndex;i<endIndex;i++)
        {
                double tr=(positionX[i]-4)*(positionX[i]-4)+(positionY[i])*(positionY[i])+positionZ[i]*positionZ[i];
                tr=fabs(sqrt(tr)-r);
                if(tr<dis)
                {
                        dis=tr;
                        dx=localParSpacing[i];
                }
        }
        dr=5*dx;
        for(size_t i=startIndex;i<endIndex;i++)
        {
                double tr=(positionX[i]-4)*(positionX[i]-4)+(positionY[i])*(positionY[i])+positionZ[i]*positionZ[i];
                tr=fabs(sqrt(tr)-r);
                if(tr<dr)
                {
                        double vr=velocityU[i]*(positionX[i]-4)+velocityV[i]*(positionY[i])+velocityW[i]*positionZ[i];
                        vr=vr/sqrt((positionX[i]-4)*(positionX[i]-4)+(positionY[i])*(positionY[i])+positionZ[i]*positionZ[i]);
                        mfr+=mass[i]*vr;
                }
        }
        mfr=mfr/2/dr;
        double mfr_2=mfr;


        r=1;
        dis=10;
        mfr=0;
        for(size_t i=startIndex;i<endIndex;i++)
        {
                double tr=(positionX[i])*(positionX[i])+(positionY[i]-4)*(positionY[i]-4)+positionZ[i]*positionZ[i];
                tr=fabs(sqrt(tr)-r);
                if(tr<dis)
                {
                        dis=tr;
                        dx=localParSpacing[i];
                }
        }
        dr=5*dx;
        for(size_t i=startIndex;i<endIndex;i++)
        {
                double tr=(positionX[i])*(positionX[i])+(positionY[i]-4)*(positionY[i]-4)+positionZ[i]*positionZ[i];
                tr=fabs(sqrt(tr)-r);
                if(tr<dr)
                {
                        double vr=velocityU[i]*(positionX[i])+velocityV[i]*(positionY[i]-4)+velocityW[i]*positionZ[i];
                        vr=vr/sqrt((positionX[i])*(positionX[i])+(positionY[i]-4)*(positionY[i]-4)+positionZ[i]*positionZ[i]);
                        mfr+=mass[i]*vr;
                }
        }
        mfr=mfr/2/dr;
        double mfr_3=mfr;
	
        fprintf(mfroutfile,"%.16g %.16g %.16g %.16g %.16g\n",time,m_pParticleData->getMassFlowRate(), mfr_1, mfr_2, mfr_3);
     
  */  
        fclose(mfroutfile);
	}
	return 0;
}

/*

int VTKParticleViewer::writeResult(double time, size_t writeStep, size_t startIndex, size_t numParticle) {
	
	// alias pointers
	const double* positionX = m_pParticleData->getPositionX();
	const double* positionY = m_pParticleData->getPositionY();
	const double* positionZ = m_pParticleData->getPositionZ();

	const double* velocityU = m_pParticleData->getVelocityU();
	const double* velocityV = m_pParticleData->getVelocityV();
	const double* velocityW = m_pParticleData->getVelocityW();

	const double* volume = m_pParticleData->getVolume();
	const double* pressure = m_pParticleData->getPressure();
	const double* soundSpeed = m_pParticleData->getSoundSpeed();

	const int* objectTag = m_pParticleData->getObjectTag();

	// Create an output file the name "filename"
	string filename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".vtk";
	FILE *outfile;
	outfile = fopen(filename.c_str(), "w");
	if(outfile==nullptr) {
		printf("Error opening file: %s\n",filename.c_str()); 
		return 1;
	}
	size_t endIndex = startIndex + numParticle;


	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.16g\n",time);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	
	fprintf(outfile,"POINTS %ld double\n",numParticle);
	if(m_pParticleData->getDimension()==2) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",positionX[i], positionY[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",positionX[i], positionY[i], positionZ[i]);
	}
	fprintf(outfile,"POINT_DATA %ld\n",numParticle);
	
	fprintf(outfile,"VECTORS Velocity double\n");
	if(m_pParticleData->getDimension()==2) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",velocityU[i], velocityV[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",velocityU[i], velocityV[i], velocityW[i]);
	}

	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",pressure[i]);
		
	fprintf(outfile,"SCALARS volume double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",volume[i]);
	
	fprintf(outfile,"SCALARS velocity_u double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",velocityU[i]);
	
	fprintf(outfile,"SCALARS velocity_v double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",velocityV[i]);
	
	if(m_pParticleData->getDimension()==3) {
		fprintf(outfile,"SCALARS velocity_w double\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g\n",velocityW[i]);
	}
	
	fprintf(outfile,"SCALARS object_tag int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",objectTag[i]);
	
	fprintf(outfile,"SCALARS sound_speed double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",soundSpeed[i]);

		

	fclose(outfile);
	
	return 0;
	
	
}
*/

////////////////////////////////////////////////////////////////////////////////////////
// End of VTKParticleViewer
////////////////////////////////////////////////////////////////////////////////////////







////////////////////////////////////////////////////////////////////////////////////////
// Start of TXTParticleViewer1D
////////////////////////////////////////////////////////////////////////////////////////

TXTParticleViewer1D::TXTParticleViewer1D(ParticleData* data, const std::string& particleType, 
const string& outputfileName, int numDigits) {
	m_pParticleData = data;
	m_sOutputfileName = outputfileName;
	m_iNumDigits = numDigits; 
	m_sParticleType = particleType;
}


int TXTParticleViewer1D::writeResult(double time, size_t writeStep) {
	 
	// alias pointers
	const double* positionX  = m_pParticleData->getPositionX();	
	const double* velocityU  = m_pParticleData->getVelocityU();	
	const double* volume     = m_pParticleData->getVolume();
	const double* pressure   = m_pParticleData->getPressure();
	const double* soundSpeed = m_pParticleData->getSoundSpeed();
	const double* dd1        = m_pParticleData->getDD1();
	const double* dd2_left   = m_pParticleData->getDD2Left();
	const double* dd2_right  = m_pParticleData->getDD2Right();
	const double* dd3_left   = m_pParticleData->getDD3Left();
	const double* dd3_right  = m_pParticleData->getDD3Right();
	
	
	// Create an output file with the name "filename"
	string filename = m_sOutputfileName + rightFlush(writeStep, m_iNumDigits) + ".txt";
	ofstream outfile;
	outfile.open(filename.c_str());
	
	size_t startIndex, numParticle;
	if(m_sParticleType=="all") {
		startIndex = 0;
		numParticle = m_pParticleData->getTotalNum();
	}
	else if(m_sParticleType=="fluid") {
		startIndex = 1;
		numParticle = m_pParticleData->getTotalNum()-1;
	}
	size_t endIndex = startIndex + numParticle;


	if(outfile.is_open()) {
		outfile<<"x"<<setw(24)
			   <<"V"<< setw(24) 
			   <<"vel"<< setw(24) 
			   <<"p"<<setw(24)
			   <<"cs"<<setw(24)
			   <<"dd1[i]"<<setw(24)
			   <<"dd2_left[i]"<<setw(24)
			   <<"dd2_right[i]"<<setw(24)
			   <<"dd3_left[i]"<<setw(24)
			   <<"dd3_right[i]"<<setw(24)
			   <<endl;
		for(size_t i=startIndex; i<endIndex; i++) {	
			outfile.precision(15);
			outfile<<left<<setw(24) 
			<<positionX[i]<<setw(24)
			<<volume[i] <<setw(24) 
			<<velocityU[i]<<setw(24) 
			<<pressure[i]<<setw(24)  
			<<soundSpeed[i]<<setw(24)
			<<dd1[i]<<setw(24)
			<<dd2_left[i]<<setw(24)
			<<dd2_right[i]<<setw(24)
			<<dd3_left[i]<<setw(24)
			<<dd3_right[i]<<setw(24)
			<<endl;
		}
		
	}
	else {
		cout<<"Unable to open file: "<<filename<<endl;
		return 1; // ERROR
	}

	outfile.close();	
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of TXTParticleViewer1D
////////////////////////////////////////////////////////////////////////////////////////
