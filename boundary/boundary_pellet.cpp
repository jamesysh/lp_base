#include "boundary_pellet.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

PelletInflowBoundary::PelletInflowBoundary():Pinflow(15), Uinflow(0), Vinflow(55){}

double calculateMassFlowRate(double energy){
	return energy;
}

int PelletInflowBoundary::UpdateInflowBoundary(ParticleData* m_pParticleData, EOS* m_pEOS, double dt, double dx){
        size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
        size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
        size_t inflowEndIndex = fluidEndIndex + m_pParticleData->m_iInflowNum;
        double *x = m_pParticleData->m_vPositionX;
        double *y = m_pParticleData->m_vPositionY;
        double *z = m_pParticleData->m_vPositionZ;
        double *pressure = m_pParticleData->m_vPressure;
        double *vx = m_pParticleData->m_vVelocityU;
        double *vy = m_pParticleData->m_vVelocityV;
        double *vz = m_pParticleData->m_vVelocityW;
        double *volume = m_pParticleData->m_vVolume;
        double *volumeold = m_pParticleData->m_vVolumeOld;
        double *localParSpacing = m_pParticleData->m_vLocalParSpacing;
        double *mass = m_pParticleData->m_vMass;
       double *sound = m_pParticleData->m_vSoundSpeed;

      double *m_vmassflowrate = m_pParticleData->m_vMassFlowRate; 
        int *pelletstate = m_pParticleData->m_vPelletState; 

	int *pelletid = m_pParticleData->m_vPelletID;
	size_t pelletn = m_pParticleData->m_iNumberofPellet;
	double *pelletx = m_pParticleData->m_vPelletPositionX;
	double *pellety = m_pParticleData->m_vPelletPositionY;
	double *pelletz = m_pParticleData->m_vPelletPositionZ;
	double *pelletr = m_pParticleData->m_vPelletRadius;
	double *pelletir = m_pParticleData->m_vPelletInnerRadius;
	double sublimationenergy = m_pParticleData->sublimationenergy;
	double *pelletvelocity = m_pParticleData->m_vPelletVelocity;
    double *volumeOnBoundary = m_pParticleData->volumeOnBoundary;
    double *ssOnBoundary = m_pParticleData->ssOnBoundary;
    double *pelletqsum = m_pParticleData->pelletqsum;
    double *pressureOnBoundary = m_pParticleData->pressureOnBoundary;
    double *uOnBoundary = m_pParticleData->uOnBoundary;
    static bool FIRST = true;    
    vector<double> m_voldv(pelletn,0);
    static double t_total = 0;
    t_total += dt; 
    	double R = 83.1446/20.1797;
    	double Ts = Vinflow*Pinflow/R;
    	double gamma0 = 1.67;
	double gamma1 = gamma0 - 1;


	for(int pi=0;pi<pelletn;pi++)
    {

	    	cout<<"average volume on boundary is "<<volumeOnBoundary[pi]<<endl;
        	cout<<"average pressure on boundary is "<<pressureOnBoundary[pi]<<endl;
        	cout<<"average soundspeed on boundary is "<<ssOnBoundary[pi]<<endl;
            cout<<"average radial velocity = "<<uOnBoundary[pi]<<endl;
        

        double pr=pelletr[pi];
       	pelletir[pi] = pr;
		double pir=pelletir[pi];
		

		double massflowrate = m_vmassflowrate[pi];
         
      
        m_voldv[pi] = pelletvelocity[pi];

        if(FIRST)
            {   
                
                FIRST = false;
                m_voldv[pi] = 0;
                
                }

		double B = (pressureOnBoundary[pi]+dt*gamma1*(pelletqsum[pi]))*volumeOnBoundary[pi]/ssOnBoundary[pi] - uOnBoundary[pi];
		double C = -(massflowrate/4/M_PI/pr/pr)*R*Ts*volumeOnBoundary[pi]/ssOnBoundary[pi];
	
		cout<<"B ="<<B<<" "<<"C ="<<C<<endl;
        if(massflowrate == 0)
		{
			pelletvelocity[pi] = uOnBoundary[pi];
			volumeOnBoundary[pi] = Vinflow;
		    pressureOnBoundary[pi] = Pinflow;
		}
		else
		{

      // volumeOnBoundary[pi] = Vinflow;
        //    pelletvelocity[pi]=massflowrate*Vinflow/4.0/M_PI/pr/pr;
       // pressureOnBoundary[pi] = Pinflow;

         pelletvelocity[pi] = (-B + sqrt(B*B - 4*C))/2;
            volumeOnBoundary[pi] = pelletvelocity[pi]/(massflowrate/4/M_PI/pr/pr);
			pressureOnBoundary[pi] = R*Ts/volumeOnBoundary[pi];
		}
           
           cout<<"volume on boundary is "<<volumeOnBoundary[pi]<<endl;
        	cout<<"pressure on boundary is "<<pressureOnBoundary[pi]<<endl;
          cout<<"pellet velocity is "<<pelletvelocity[pi]<<endl;
    cout<<"massflowrate is "<<massflowrate<<endl;            
	 

 
 
 } 
	    	inflowEndIndex = fluidEndIndex; 
 //generate new inflow particles in region 0<r<pir

    for(int pi=0;pi<pelletn;pi++){
    	double pr=pelletr[pi];
       	pelletir[pi] = pr;
		double pir=pelletir[pi];
		
    int n=4.0*3.1416*pr*pr*pr/dx/dx/dx*sqrt(2.0)*Vinflow/volumeOnBoundary[pi];
//		cout<<n<<endl;
		if(inflowEndIndex+n>=m_pParticleData->m_iCapacity) {
			cout<<"Error: too many inflow particles: n = "<<n<<endl;
			return 1;//too many
		}
		double newpir=pir-n*dx*dx*dx/4.0/3.1416/pr/pr/sqrt(2.0)*volumeOnBoundary[pi]/Vinflow;
        	for(int i=0;i<n;i++)
		{
			double tx=1,ty=1,tz=1,tr;
			while(tx*tx+ty*ty+tz*tz>1)
			{
				tx=2*(double)rand()/(double)RAND_MAX-1;
				ty=2*(double)rand()/(double)RAND_MAX-1;
				tz=2*(double)rand()/(double)RAND_MAX-1;
			}
			tr=sqrt(tx*tx+ty*ty+tz*tz);
			tx=tx/tr,ty=ty/tr,tz=tz/tr;
			tr=newpir+(double)rand()/(double)RAND_MAX*(pir-newpir);
			x[inflowEndIndex]=pelletx[pi]+tx*tr;
			y[inflowEndIndex]=pellety[pi]+ty*tr;
			z[inflowEndIndex]=pelletz[pi]+tz*tr;
            double d_x = tx*tr;
            double d_y = ty*tr;
            double d_z = tz*tr;
            double dr = sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
	       		vx[inflowEndIndex] = m_voldv[pi]*d_x/dr;
           		vy[inflowEndIndex] = m_voldv[pi]*d_y/dr;
           		vz[inflowEndIndex] = m_voldv[pi]*d_z/dr;
			pelletid[inflowEndIndex]=pi;
			inflowEndIndex++;
		}
	}
	
//	cout<<"move particles"<<endl;
	for(size_t index=fluidEndIndex;index<inflowEndIndex;index++)
	{
        	int pi = pelletid[index];
        	double d_x=x[index]-pelletx[pi];
		double d_y=y[index]-pellety[pi];
		double d_z=z[index]-pelletz[pi];
		double dr=sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
		double pir=pelletir[pi];
		double pr=pelletr[pi];
        double old_v=sqrt(vx[index]*vx[index]+vy[index]*vy[index]+vz[index]*vz[index]);
		double newv=pelletvelocity[pi];
    //		volumeold[index] = volume[index]=volumeOnBoundary[pi];
  //     	pressure[index] = pressureOnBoundary[pi];
    //   	localParSpacing[index]=dx;
	//	mass[index]=dx*dx*dx/Vinflow/sqrt(2);
     //  	sound[index]=m_pEOS->getSoundSpeed(pressure[index],1./volume[index]);
       	x[index]+=dt*0.5*(old_v+newv)*d_x/dr;
		y[index]+=dt*0.5*(old_v+newv)*d_y/dr;
		z[index]+=dt*0.5*(old_v+newv)*d_z/dr;
		vx[index]=newv*d_x/dr;
		vy[index]=newv*d_y/dr;
		vz[index]=newv*d_z/dr;
		dr+=dt*0.5*(old_v+newv);
		
        volumeold[index] = volume[index]=volumeOnBoundary[pi]*dr*dr/pr/pr;
       	pressure[index] = pressureOnBoundary[pi]*pr*pr/dr/dr;
       	localParSpacing[index]=dx;
		mass[index]=dx*dx*dx/Vinflow/sqrt(2);
       	sound[index]=m_pEOS->getSoundSpeed(pressure[index],1./volume[index]);


        if(dr>pelletr[pelletid[index]])//change label from inflow to fluid if r>pr
		{
		 	volumeold[index] = volume[index] = volumeOnBoundary[pi];
         	pressure[index] = pressureOnBoundary[pi];
       		sound[index]=m_pEOS->getSoundSpeed(pressure[index],1./volume[index]);
       		if(index>fluidEndIndex)
                        	m_pParticleData->swap(index,fluidEndIndex);
			fluidEndIndex++;
            pelletid[fluidEndIndex-1] = -1; 
		}
	}

//	cout<<"inflow finished"<<endl;
        m_pParticleData->m_iFluidNum=fluidEndIndex-fluidStartIndex;
        m_pParticleData->m_iBoundaryStartIndex=fluidEndIndex;
        m_pParticleData->m_iBoundaryNum=m_pParticleData->m_iInflowNum=inflowEndIndex-fluidEndIndex;
        m_pParticleData->m_iGhostStartIndex=inflowEndIndex;
        m_pParticleData->m_iTotalNum=inflowEndIndex-fluidStartIndex;
	return 0;

}

PelletOutflowBoundary::PelletOutflowBoundary():xmin(-16),xmax(16),ymin(-16),ymax(16),zmin(-16),zmax(16) {
}

int PelletOutflowBoundary::UpdateInflowBoundary(ParticleData *m_pParticleData, EOS* m_pEOS, double dt, double dx) {
        size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
        size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
        size_t inflowEndIndex = fluidEndIndex + m_pParticleData->m_iInflowNum;
        double *x = m_pParticleData->m_vPositionX;
        double *y = m_pParticleData->m_vPositionY;
        double *z = m_pParticleData->m_vPositionZ;
        int* timetrack = m_pParticleData->m_vTimeTrack;
//        double *volume = m_pParticleData->m_vVolume;
//        double volumeforvacuum = 1e6;
     
     
      
        for(size_t index=fluidStartIndex;index<fluidEndIndex;index++)
        {
                if(x[index]<xmin || x[index]>xmax || y[index]<ymin || y[index]>ymax || z[index]<zmin || z[index]>zmax)
                {
                        
                        if(index+1<fluidEndIndex)
                        {
                                m_pParticleData->swap(index,fluidEndIndex-1);
                          }
                        
			            timetrack[fluidEndIndex-1] = 0; //RESET TIMETRACK TO BE ZERO
                        if(fluidEndIndex<inflowEndIndex)
                        {
                                m_pParticleData->swap(fluidEndIndex-1,inflowEndIndex-1);
                        }
                        fluidEndIndex--;
                        inflowEndIndex--;
	             }
        }
        m_pParticleData->m_iFluidNum=fluidEndIndex-fluidStartIndex;
        m_pParticleData->m_iBoundaryStartIndex=fluidEndIndex;
        m_pParticleData->m_iGhostStartIndex=inflowEndIndex;
        m_pParticleData->m_iTotalNum=inflowEndIndex-fluidStartIndex;
        return 0;
}

