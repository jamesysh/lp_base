#include "boundary_pellet.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <time.h>
using namespace std;

float     Bessel_Kn( int n,float  x);
float     Bessel_I0(float  x);
float     Bessel_K1(float  x);
float     Bessel_I1(float  x);


PelletInflowBoundary::PelletInflowBoundary():Pinflow(30), Uinflow(0), Vinflow(100.0){}

double calculateMassFlowRate(double energy){
	return energy;
}

int PelletInflowBoundary::UpdateInflowBoundary(ParticleData* m_pParticleData, EOS* m_pEOS, double dt, double dx){
        size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
        size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
        size_t inflowEndIndex = fluidEndIndex;
        //+ m_pParticleData->m_iInflowNum;
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
	double masse = m_pParticleData->masse;
	double massNe = m_pParticleData->massNe;
	double teinf = m_pParticleData->teinf;
	double INe = m_pParticleData->INe;
	int ZNe = m_pParticleData->ZNe;
	double neinf = m_pParticleData->neinf;
	double heatK = m_pParticleData->heatK;
	double e = heatK*(2.99792458e7)/100;
    double lnLambda = log(2*teinf/INe*sqrt(exp(1)/2));
	const double* leftintegral = m_pParticleData->m_vLeftIntegral;
	const double* rightintegral = m_pParticleData->m_vRightIntegral;
	
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
//	cout<<"insert particles"<<endl;
//	cout<<inflowEndIndex<<endl;
//	cout<<m_pParticleData->m_iCapacity<<endl;
        
        double R = 83.1446/20.1797;
        double Ts = Vinflow*Pinflow/R;

        double gamma0 = 5.0/3;
   

    for(int pi=0;pi<pelletn;pi++){
		double pir=pelletir[pi];
		double pr=pelletr[pi];
           	for(size_t index=fluidStartIndex;index<fluidEndIndex;index++)
		{
			double d_x=x[index]-pelletx[pi];
			double d_y=y[index]-pellety[pi];
			double d_z=z[index]-pelletz[pi];
			double r=d_x*d_x+d_y*d_y+d_z*d_z;
			if(r<(pr+dx)*(pr+dx) && r>pr*pr)
			{   
                double tauleft = leftintegral[index]/massNe*ZNe;
                double tauright = rightintegral[index]/massNe*ZNe;
                double tauinf = heatK*heatK*teinf*teinf/(8.0*3.1416*e*e*e*e*lnLambda);
                double taueff = tauinf*sqrt(2.0/(1.0+ZNe));
                double uleft = tauleft/taueff;
                double uright = tauright/taueff;
                double qinf=sqrt(2.0/3.1416/masse)*neinf*pow(heatK*teinf,1.5);
    
 /*               if(d_x>0)
				    pelletqsum[pi] += qinf*0.5*uright*Bessel_Kn(2,sqrt(uright));
                
                else
                    pelletqsum[pi] += qinf*0.5*uleft*Bessel_Kn(2,sqrt(uleft));
   */          
                pelletqsum[pi] += qplusminus[index];

                volumeOnBoundary[pi] += volume[index];
                pressureOnBoundary[pi] += pressure[index];
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
		pellete[pi]=pelletqsum[pi]/pelletneighbor[pi]*4*3.1416*pr*pr;
        volumeOnBoundary[pi] = volumeOnBoundary[pi]/pelletneighbor[pi];
        pressureOnBoundary[pi] /= pelletneighbor[pi];
        volumeOnBoundary[pi] = pow(pow(volumeOnBoundary[pi],gamma0)*pressureOnBoundary[pi]/R/Ts,1./(gamma0-1)); 
        cout<<"volume on boundary is "<<volumeOnBoundary[pi]<<endl;
        cout<<"pressure on boundary is "<<pressureOnBoundary[pi]<<endl;
        double massflowrate=pellete[pi]/sublimationenergy;
		cout<<"Mass flow rate = "<<massflowrate<<endl;

		m_pParticleData->m_vMassFlowRate = massflowrate;

		double oldv=pelletvelocity[pi];
		pelletvelocity[pi]=massflowrate*Vinflow/4.0/3.1416/pr/pr;
		cout<<"pellet ablation velocity = "<<pelletvelocity[pi]<<endl;
        
//		cout<<pi<<endl;
    pir = pr;
		//generate new inflow particles in region 0<r<pir
   if(pir > 3*dx/2){
    int layer_n = floor((pir)/dx) ;
    cout<<"layer_n "<<layer_n<<endl;
    double layer_r[layer_n];
     if(layer_n != 1){
   for(int i=0;i<layer_n;i++){
        layer_r[i] = pir - dx/2 - i*dx;
        if (layer_r[i] < dx/2 ) 
        {   
            layer_r[i] = dx/2;
            }
        }
     }
     else{
     layer_r[0] = dx/2;
     }
    int numberOfParticleOnLayer = sqrt(2)*4*M_PI*(pr)*(pr)/dx/dx;//*Vinflow/volumeOnBoundary[pi];


    int n_tmp = layer_n*numberOfParticleOnLayer;
     
 	if(inflowEndIndex+n_tmp>=m_pParticleData->m_iCapacity) {
			cout<<"Error: too many inflow particles: n = "<<n_tmp<<endl;
			return 1;//too many
		}
	
    
    srand(time(NULL));
    for(int layer_id=0;layer_id<layer_n;layer_id++){

         
       int rnd = rand();
       double offset = 2./numberOfParticleOnLayer;
       double increment = M_PI*(3-sqrt(5));
       for(int i=0;i<numberOfParticleOnLayer;i++){
           double r_random = (2*(double)rand()/(double)RAND_MAX - 1)*dx*0.5; 
           double  y_tmp =  ((i*offset-1)+offset/2);
           y[inflowEndIndex] = y_tmp * (layer_r[layer_id]+r_random) + pellety[pi];
           double r = sqrt(1-y_tmp*y_tmp);
           double phi = ((i+rnd)%numberOfParticleOnLayer) * increment;

           x[inflowEndIndex] = cos(phi)*r*(layer_r[layer_id]+r_random) + pelletx[pi];
           z[inflowEndIndex] = sin(phi)*r*(layer_r[layer_id]+r_random) + pelletz[pi];
           double d_x = x[inflowEndIndex] - pelletx[pi];
           double d_y = y[inflowEndIndex] - pellety[pi];
           double d_z = z[inflowEndIndex] - pelletz[pi];
           double dr = sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
           vx[inflowEndIndex] = oldv*d_x/dr;
           vy[inflowEndIndex] = oldv*d_y/dr;
           vz[inflowEndIndex] = oldv*d_z/dr;
   
//			    volumeold[inflowEndIndex]=volume[inflowEndIndex]=Vinflow;
  //         pressure[inflowEndIndex]=Pinflow;
			    localParSpacing[inflowEndIndex]=dx;
			    mass[inflowEndIndex]=dx*dx*dx/sqrt(2.0)/Vinflow;
	//		    sound[inflowEndIndex]=m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
	
                pelletid[inflowEndIndex]=pi;
			    inflowEndIndex++;


                 }

        }

    }
	 

    }
//	cout<<"move particles"<<endl;

 

	for(size_t index=fluidEndIndex;index<inflowEndIndex;index++)
	{
		//move particle[index]
       
        int pelleti = pelletid[index];
        double d_x=x[index]-pelletx[pelletid[index]];
		double d_y=y[index]-pellety[pelletid[index]];
		double d_z=z[index]-pelletz[pelletid[index]];
		double dr=sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
		double oldv=sqrt(vx[index]*vx[index]+vy[index]*vy[index]+vz[index]*vz[index]);
		double newv=pelletvelocity[pelletid[index]];
  
        double Tb = Ts;//- (gamma0-1)/(2*gamma0*R)*newv*newv;
        double mach = newv/sqrt(gamma0*R*Ts);
        if(index == fluidEndIndex){cout<<"mach "<<mach<<endl;}
       if(mach<1){

           
           volume[index]  = volumeold[index] = Vinflow;//dr*dr/pelletr[pelleti]/pelletr[pelleti]*volumeOnBoundary[pelleti]; 
     
           pressure[index] = R*Tb/volume[index];
           sound[index] = m_pEOS->getSoundSpeed(pressure[index],1./volume[index]);

        } 
        
        else{
                pressure[index]=Pinflow;
			    volumeold[index]=volume[index]=Vinflow;
			    sound[index]=m_pEOS->getSoundSpeed(pressure[index],1./volume[index]);
	          
            newv = newv*Vinflow/volumeOnBoundary[pelleti];
        
        }
        
        x[index] += dt*0.5*(oldv+newv)*d_x/dr;
		y[index] += dt*0.5*(oldv+newv)*d_y/dr;
		z[index] += dt*0.5*(oldv+newv)*d_z/dr;
		vx[index] = newv*d_x/dr;
		vy[index] = newv*d_y/dr;
		vz[index] = newv*d_z/dr;
		dr+=dt*0.5*(oldv+newv);
		if(dr>pelletr[pelleti])//change label from inflow to fluid if r>pr
		{ 
       /* 
           volume[index]  = volumeold[index] = volumeOnBoundary[pelleti]; 
           pressure[index] = R*Tb/volume[index];
           sound[index] = m_pEOS->getSoundSpeed(pressure[index],1./volume[index]);
*/

			if(index>fluidEndIndex)
          {  m_pParticleData->swap(index,fluidEndIndex);
          }
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
                        if(fluidEndIndex<inflowEndIndex)
                        {
                                m_pParticleData->swap(fluidEndIndex-1,inflowEndIndex-1);
                        }
                        fluidEndIndex--;
                        inflowEndIndex--;
			            timetrack[fluidEndIndex] = 0; //RESET TIMETRACK TO BE ZERO
	             }
        }
        m_pParticleData->m_iFluidNum=fluidEndIndex-fluidStartIndex;
        m_pParticleData->m_iBoundaryStartIndex=fluidEndIndex;
        m_pParticleData->m_iGhostStartIndex=inflowEndIndex;
        m_pParticleData->m_iTotalNum=inflowEndIndex-fluidStartIndex;
        return 0;
}

