#include "boundary_pellet.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

PelletInflowBoundary::PelletInflowBoundary():Pinflow(30), Uinflow(0), Vinflow(100){}

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
	double *pellete = m_pParticleData->m_vPelletEnergy;
	double *qplusminus = m_pParticleData->m_vQplusminus;
	double sublimationenergy = m_pParticleData->sublimationenergy;
	double *pelletvelocity = m_pParticleData->m_vPelletVelocity;
	vector<double> pelletqsum(pelletn,0);
	vector<int> pelletneighbor(pelletn,0);
    vector<double> volumeOnBoundary(pelletn,0);
    vector<double> pressureOnBoundary(pelletn,0);
    vector<double> m_voldv(pelletn,0);
	vector<double> uOnBoundary(pelletn,0);
	vector<double> ssOnBoundary(pelletn,0);	
    static double t_total = 0;
    t_total += dt; 
    	double R = 83.1446/20.1797;
    	double Ts = 200;//Vinflow*Pinflow/R;
    	double gamma0 = 1.67;
	double gamma1 = gamma0 - 1;

//	cout<<"insert particles"<<endl;
//	cout<<inflowEndIndex<<endl;
//	cout<<m_pParticleData->m_iCapacity<<endl;
	for(int pi=0;pi<pelletn;pi++)
	{
//		cout<<pi<<endl;
		double pr=pelletr[pi];
       	pelletir[pi] = pr;
		double pir=pelletir[pi];
		//calculate ablation velocity
/*
    for(size_t index=fluidStartIndex;index<inflowEndIndex;index++)
		{
			
            double d_x=x[index]-pelletx[pi];
			double d_y=y[index]-pellety[pi];
			double d_z=z[index]-pelletz[pi];
			double r=d_x*d_x+d_y*d_y+d_z*d_z;

			
			if( r>pr*pr && r<(pr+dx)*(pr+dx))
            {
				pelletneighbor[pi]++;
			}
		}
		if(pelletneighbor[pi]==0)
		{
			cout<<"Error: cannot find neighbor for pellet"<<endl;
			assert(false);
			return 0;
		}
		cout<<"Number of neighbor for pellet = "<<pelletneighbor[pi]<<endl;

        pellete[pi]=pelletqsum[pi]/pelletneighbor[pi]*4*M_PI*pr*pr;
		cout<<"pellete is "<<pellete[pi]<<endl;
   
   double massflowrate=pellete[pi]/sublimationenergy;

		m_vmassflowrate[pi] = massflowrate;
        cout<<"massflowrate is "<<massflowrate<<endl;
        pelletneighbor[pi] = 0;
*/
    for(size_t index=fluidStartIndex;index<inflowEndIndex;index++)
		{
			
            double ss = sound[index];
            double shift = ss*dt;
            double d_x=x[index]-pelletx[pi];
			double d_y=y[index]-pellety[pi];
			double d_z=z[index]-pelletz[pi];
			double r=d_x*d_x+d_y*d_y+d_z*d_z;
            double r_shift = (sqrt(r)-shift)*(sqrt(r)-shift);

			if(r_shift<(pr+dx/5)*(pr+dx/5) && r_shift > (pr-dx/5)*(pr-dx/5) && r>(pr)*(pr) && r<(pr+dx)*(pr+dx))
			
		//	if( r>pr*pr && r<(pr+dx/2)*(pr+dx/2))
            {
                
				pelletqsum[pi]+=qplusminus[index];
                volumeOnBoundary[pi] += volume[index];
                pressureOnBoundary[pi] += pressure[index];
				ssOnBoundary[pi] += sound[index];
				double vr = (vx[index]*d_x + vy[index]*d_y + vz[index]*d_z)/sqrt(r);
				uOnBoundary[pi] += vr;
				pelletneighbor[pi]++;
			}
		}
		if(pelletneighbor[pi]==0)
		{
			cout<<"Error: cannot find neighbor for pellet"<<endl;
			assert(false);
			return 0;
		}
		cout<<"Number of neighbor for pellet = "<<pelletneighbor[pi]<<endl;
      pellete[pi]=pelletqsum[pi]/pelletneighbor[pi]*4*M_PI*pr*pr*2/M_PI;
		cout<<"pellete is "<<pellete[pi]<<endl;
   
   double massflowrate=pellete[pi]/sublimationenergy;

		m_vmassflowrate[pi] = massflowrate;
        cout<<"massflowrate is "<<massflowrate<<endl;
    
        volumeOnBoundary[pi] = volumeOnBoundary[pi]/pelletneighbor[pi];
		
        pressureOnBoundary[pi] = pressureOnBoundary[pi]/pelletneighbor[pi];
        uOnBoundary[pi] = uOnBoundary[pi]/pelletneighbor[pi];
        if(uOnBoundary[pi]<0) 
            uOnBoundary[pi] = 0.;
        //volumeOnBoundary[pi] = pow(pow(volumeOnBoundary[pi],gamma0)*pressureOnBoundary[pi]/R/Ts,1./(gamma0-1)); 
      
      m_voldv[pi] = pelletvelocity[pi];

		double B = (pressureOnBoundary[pi]+dt*gamma1*(pelletqsum[pi]/pelletneighbor[pi]))*volumeOnBoundary[pi]/ssOnBoundary[pi] - uOnBoundary[pi];
		double C = -(massflowrate/4/M_PI/pr/pr)*R*Ts*volumeOnBoundary[pi]/ssOnBoundary[pi];
	
		cout<<"B ="<<B<<" "<<"C ="<<C<<endl;
        cout<<"u average is "<<uOnBoundary[pi]<<endl;	
		if(massflowrate == 0)
		{
			pelletvelocity[pi] = 0.0;
			volumeOnBoundary[pi] = Vinflow;
			pressureOnBoundary[pi] = Pinflow;
		}
		else
		{
			// pelletvelocity[pi] =volumeOnBoundary[pi]*(massflowrate/4/M_PI/pr/pr);

           pelletvelocity[pi] = (-B + sqrt(B*B - 4*C))/2;
			volumeOnBoundary[pi] = pelletvelocity[pi]/(massflowrate/4/M_PI/pr/pr);
			pressureOnBoundary[pi] = R*Ts/volumeOnBoundary[pi];
		}
           double D =  (pressureOnBoundary[pi]+dt*gamma1*(pelletqsum[pi]/pelletneighbor[pi]))*volumeOnBoundary[pi]/ssOnBoundary[pi]; 
           cout<<"D is "<<D<<endl; 
        //pelletvelocity[pi]=massflowrate*volumeOnBoundary[pi]/4.0/M_PI/pr/pr;
	    	cout<<"volume on boundary is "<<volumeOnBoundary[pi]<<endl;
        	cout<<"pressure on boundary is "<<pressureOnBoundary[pi]<<endl;
        	cout<<"pellet ablation velocity = "<<pelletvelocity[pi]<<endl;

 } 
/* 
    for(size_t index=fluidEndIndex;index<inflowEndIndex;index++)
	    {
        	int pi = pelletid[index];
        	double d_x=x[index]-pelletx[pi];
		double d_y=y[index]-pellety[pi];
		double d_z=z[index]-pelletz[pi];
		double dr=sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
		double pir=pelletir[pi];
		double pr=pelletr[pi];
		if(dr>pelletr[pi])//change label from inflow to fluid if r>pr
		{
		 //	volumeold[index] = volume[index]=volumeOnBoundary[pi];
         //	pressure[index] = pressureOnBoundary[pi];
       	//	sound[index]=m_pEOS->getSoundSpeed(pressure[index],1./volume[index]);
       		if(index>fluidEndIndex)
                        	m_pParticleData->swap(index,fluidEndIndex);
			fluidEndIndex++;
            pelletid[fluidEndIndex-1] = -1; 
		}
	}

 */
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
//			if(i%1000==0)
//			cout<<i<<endl;
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
	       		vx[inflowEndIndex] = m_voldv[pi]*d_x/dr*dr*dr/pr/pr;
           		vy[inflowEndIndex] = m_voldv[pi]*d_y/dr*dr*dr/pr/pr;
           		vz[inflowEndIndex] = m_voldv[pi]*d_z/dr*dr*dr/pr/pr;
        volumeold[inflowEndIndex] = volume[inflowEndIndex]=volumeOnBoundary[pi];
       	pressure[inflowEndIndex] = pressureOnBoundary[pi];
       	localParSpacing[inflowEndIndex]=dx;
		mass[inflowEndIndex]=dx*dx*dx/Vinflow/sqrt(2);
       	sound[inflowEndIndex]=m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);

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
		if(dr>pelletr[pelletid[index]])//change label from inflow to fluid if r>pr
		{
		 	volumeold[index] = volume[index]=volumeOnBoundary[pi];
         	pressure[index] = pressureOnBoundary[pi];;
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

PelletOutflowBoundary::PelletOutflowBoundary():xmin(-20),xmax(20),ymin(-20),ymax(20),zmin(-20),zmax(20) {
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

