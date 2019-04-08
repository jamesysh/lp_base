#include "particle_data.h"
#include "initializer.h"
#include "boundary.h"
#include <algorithm>
#include <iostream>
//#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Start : ParticleData
////////////////////////////////////////////////////////////////////////////////


ParticleData::ParticleData(Initializer& init) {
	
	// get the scalar parameters
	m_iDimension = init.getDimension(); 
	m_iCapacity = init.getCapacity();
	m_iMaxNeighbourNum = init.getMaxNeighbourNum(); 
//	m_iMaxNeighbourNumInOneDir = init.getMaxNeighbourNumInOneDir();	
        m_iMaxNeighbourNumInOneDir = init.getMaxNeighbourNum(); 

//	m_iUseLimiter = init.getUseLimiter();	


	m_iFluidNum = init.getFluidNum();
//	m_iBoundaryNum = init.getBoundaryNum();
	m_iInflowNum = 0;
	m_iBoundaryNum = 0;
	m_iGhostNum = 0;
//	m_iTotalNum = m_iFluidNum + m_iBoundaryNum + m_iGhostNum;
	m_iTotalNum = m_iFluidNum;

	m_iFluidStartIndex = init.getFluidStartIndex();
	m_iBoundaryStartIndex = m_iFluidStartIndex + m_iFluidNum;
	m_iGhostStartIndex = m_iFluidStartIndex + m_iFluidNum;

	// get the resources allocated on heap by the Initializer class
	m_vPositionX = init.getPositionX();
	m_vPositionY = init.getPositionY();	
	m_vPositionZ = init.getPositionZ();
	m_vVelocityU = init.getVelocityU();
	m_vVelocityV = init.getVelocityV();
	m_vVelocityW = init.getVelocityW();
	m_vVolume = init.getVolume();
	m_vPressure = init.getPressure();
	m_vSoundSpeed = init.getSoundSpeed();	
	m_vTemperature = init.getTemperature();
    m_vTimeTrack = init.getTimeTrack();
    m_vObjectTag = init.getObjectTag();
	m_vLocalParSpacing = init.getLocalParSpacing();
	m_vMass=init.getMass();

	m_iNumberofPellet = init.getNumberofPellet();
	m_iMaxParticlePerCell=init.getMaxParticlePerCell();
	m_iQuadtreeResolution = init.getQuadtreeResolution();
    m_iBinarytreeResolution = init.getBinarytreeResolution();
    m_iHeatModel = init.getHeatModel();
    m_iMagneticField = init.getMagneticField();
    m_vPelletPositionX = init.getPelletPositionX();
	m_vPelletPositionY = init.getPelletPositionY();
	m_vPelletPositionZ = init.getPelletPositionZ();
	m_vPelletRadius = init.getPelletRadius();
	m_vPelletInnerRadius = init.getPelletInnerRadius();
    m_vPelletOuterRadius = init.getPelletOuterRadius();
	Magx = init.getMagx();
	Magy = init.getMagy();
	Magz = init.getMagz();
	masse = init.getmasse();
	massNe = init.getmassNe();
	teinf = init.getteinf();
	INe = init.getINe();
	ZNe = init.getZNe();
	neinf = init.getneinf();
	heatK = init.getheatK();
	one_plus_Zstar = init.getOnePlusZstar();
    conductivity = init.getconductivity();
	sublimationenergy = init.getsublimationenergy();

    //------------output option control--------------------
    m_iPrintVelocity = init.getPrintVelocity();
    m_iPrintVelocityU = init.getPrintVelocityU();
    m_iPrintVelocityV = init.getPrintVelocityV();
    m_iPrintVelocityW = init.getPrintVelocityW();
    m_iPrintVolume = init.getPrintVolume();
    m_iPrintDensity = init.getPrintDensity();
    m_iPrintPressure = init.getPrintPressure();
    m_iPrintSoundSpeed = init.getPrintSoundSpeed();
    m_iPrintMass = init.getPrintMass();
    m_iPrintLeftIntegral = init.getPrintLeftIntegral();
    m_iPrintRightIntegral = init.getPrintRightIntegral();
    m_iPrintDeltaq = init.getPrintDeltaq();
    m_iPrintQplusminus = init.getPrintQplusminus();
    m_iPrintLocalSpacing = init.getPrintLocalSpacing();
    m_iPrintTemperature = init.getPrintTemperature();
    m_iPrintTimeTrack = init.getPrintTimeTrack();

	try {
        /*
        m_vPositionXOld = new double[m_iCapacity];	
        for(size_t i=0; i<m_iCapacity; i++) m_vPositionXOld[i]=m_vPositionX[i];

        m_vPositionYOld = new double[m_iCapacity];	
        for(size_t i=0; i<m_iCapacity; i++) m_vPositionYOld[i]=m_vPositionY[i];

        m_vPositionZOld = new double[m_iCapacity];	
        for(size_t i=0; i<m_iCapacity; i++) m_vPositionZOld[i]=m_vPositionZ[i];

         m_vVelocityRadialOld = new double[m_iCapacity];
         fill_n(m_vVelocityRadialOld,m_iCapacity,0);

         m_vSoundSpeedOld = new double[m_iCapacity];	
		for(size_t i=0; i<m_iCapacity; i++) m_vSoundSpeedOld[i]=m_vSoundSpeed[i];

        m_vPressure = new double[m_iCapacity];	
		for(size_t i=0; i<m_iCapacity; i++) m_vPressureOld[i]=m_vPressure[i];

*/

        m_vVolumeOld = new double[m_iCapacity];	
		for(size_t i=0; i<m_iCapacity; i++) m_vVolumeOld[i]=m_vVolume[i];

                m_vVolumeVoronoi = new double[m_iCapacity];
                for(size_t i=0; i<m_iCapacity; i++) m_vVolumeVoronoi[i]=m_vVolume[i];
		
        		m_vFluidBoundingBox = init.getFluidBoundingBox();
		        m_vBoundaryObj = init.getBoundaryObj();
		        m_vBoundaryObjTypes = init.getBoundaryObjTypes();
                m_bLeftInflow = new bool[m_iCapacity];
                fill_n(m_bLeftInflow,m_iCapacity,false);
                    
		// allocate pellet memory 
                if(m_iNumberofPellet){
                    pelletqsum = new double[m_iNumberofPellet];
                    volumeOnBoundary = new double[m_iNumberofPellet];
                    pressureOnBoundary = new double[m_iNumberofPellet];
                    ssOnBoundary = new double[m_iNumberofPellet];
                    uOnBoundary = new double[m_iNumberofPellet];

                    m_vLeftIntegral = new double[m_iCapacity];
                    fill_n(m_vLeftIntegral, m_iCapacity, 0);
                    m_vRightIntegral = new double[m_iCapacity];
                    fill_n(m_vRightIntegral, m_iCapacity, 0);
                    m_vDeltaq = new double[m_iCapacity];
                    fill_n(m_vDeltaq, m_iCapacity, 0);
                    m_vQplusminus = new double[m_iCapacity];
                    fill_n(m_vQplusminus, m_iCapacity, 0);

                    m_vPelletEnergy = new double[m_iNumberofPellet];
                    fill_n(m_vPelletEnergy, m_iNumberofPellet, 0);
                    m_vPelletVelocity = new double[m_iNumberofPellet];
                    fill_n(m_vPelletVelocity, m_iNumberofPellet, 0);
                    m_vMassFlowRate = new double[m_iNumberofPellet];
                    fill_n(m_vMassFlowRate,m_iNumberofPellet, 0);
                    m_vPelletID = new int[m_iCapacity];
                    fill_n(m_vPelletID, m_iCapacity, -1);
                    m_vPelletState = new int[m_iNumberofPellet];
                    fill_n(m_vPelletState, m_iNumberofPellet, 1);

                }

                    m_vPhi = new double[m_iCapacity];
                    fill_n(m_vPhi,m_iCapacity,0);
                    m_vIfSPHDensity = new int[m_iCapacity];
                    fill_n(m_vIfSPHDensity,m_iCapacity,0);

#ifdef LW_DEBUG
                m_vPError0 = new double[m_iCapacity];
                fill_n(m_vPError0,m_iCapacity,0);
                m_vPError1 = new double[m_iCapacity];
                fill_n(m_vPError1,m_iCapacity,0);
                m_vVelError0 = new double[m_iCapacity];
                fill_n(m_vVelError0,m_iCapacity,0);
                m_vVelError1 = new double[m_iCapacity];
                fill_n(m_vVelError1,m_iCapacity,0);
                m_vPxl = new double[m_iCapacity];
                fill_n(m_vPxl,m_iCapacity,0);
                m_vPxr = new double[m_iCapacity];
                fill_n(m_vPxr,m_iCapacity,0);
                m_vVxl = new double[m_iCapacity];
                fill_n(m_vVxl,m_iCapacity,0);
                m_vVxr = new double[m_iCapacity];
                fill_n(m_vVxr,m_iCapacity,0);
                m_vPyl = new double[m_iCapacity];
                fill_n(m_vPyl,m_iCapacity,0);
                m_vPyr = new double[m_iCapacity];
                fill_n(m_vPyr,m_iCapacity,0);
                m_vVyl = new double[m_iCapacity];
                fill_n(m_vVyl,m_iCapacity,0);
                m_vVyr = new double[m_iCapacity];
                fill_n(m_vVyr,m_iCapacity,0);
                m_vVtx = new double[m_iCapacity];
                fill_n(m_vVtx,m_iCapacity,0);
                m_vPtx = new double[m_iCapacity];
                fill_n(m_vPtx,m_iCapacity,0);
                m_vVty = new double[m_iCapacity];
                fill_n(m_vVty,m_iCapacity,0);
                m_vPty = new double[m_iCapacity];
                fill_n(m_vPty,m_iCapacity,0);
                m_vVolumetx = new double[m_iCapacity];
                fill_n(m_vVolumetx,m_iCapacity,0);
                m_vVolumety = new double[m_iCapacity];
                fill_n(m_vVolumety,m_iCapacity,0);
                m_vNeighSize = new int[4*m_iCapacity];
                fill_n(m_vNeighSize,4*m_iCapacity,0);
                m_vNeighList = new int[40*m_iCapacity];
                fill_n(m_vNeighList,40*m_iCapacity,0);
                m_vNeighOfParticle = new int[m_iCapacity];
                fill_n(m_vNeighOfParticle,m_iCapacity,0);
#endif

                m_vVolume_x = new double[m_iCapacity];
                m_vVolume_y = new double[m_iCapacity];
                m_vVolume_z = new double[m_iCapacity];
		m_vDensity = new double[m_iCapacity];
                fill_n(m_vVolume_x,m_iCapacity,0);
                fill_n(m_vVolume_y,m_iCapacity,0);
                fill_n(m_vVolume_z,m_iCapacity,0);
		fill_n(m_vDensity,m_iCapacity,0);
		m_vTemp1VelocityU = new double[m_iCapacity];
		if(m_iDimension==2 || m_iDimension==3) m_vTemp1VelocityV = new double[m_iCapacity];
		if(m_iDimension==3) m_vTemp1VelocityW = new double[m_iCapacity];
		m_vTemp1Volume = new double[m_iCapacity];
		m_vTemp1Pressure = new double[m_iCapacity];
		m_vTemp1SoundSpeed = new double[m_iCapacity];
		m_vTemp1PositionX = new double[m_iCapacity];
		
		fill_n(m_vTemp1VelocityU,m_iCapacity,0);
		if(m_iDimension==2 || m_iDimension==3) fill_n(m_vTemp1VelocityV,m_iCapacity,0);	
		if(m_iDimension==3) fill_n(m_vTemp1VelocityW,m_iCapacity,0);
		fill_n(m_vTemp1Volume,m_iCapacity,0);
		fill_n(m_vTemp1Pressure,m_iCapacity,0);
		fill_n(m_vTemp1SoundSpeed,m_iCapacity,0);
		fill_n(m_vTemp1PositionX,m_iCapacity,0);

		m_vTemp2VelocityU = new double[m_iCapacity];
		if(m_iDimension==2 || m_iDimension==3) m_vTemp2VelocityV = new double[m_iCapacity];
		if(m_iDimension==3) m_vTemp2VelocityW = new double[m_iCapacity];
		m_vTemp2Volume = new double[m_iCapacity];
		m_vTemp2Pressure = new double[m_iCapacity];
		m_vTemp2SoundSpeed = new double[m_iCapacity];

		fill_n(m_vTemp2VelocityU,m_iCapacity,0);
		if(m_iDimension==2 || m_iDimension==3) fill_n(m_vTemp2VelocityV,m_iCapacity,0);	
		if(m_iDimension==3) fill_n(m_vTemp2VelocityW,m_iCapacity,0);
		fill_n(m_vTemp2Volume,m_iCapacity,0);
		fill_n(m_vTemp2Pressure,m_iCapacity,0);
		fill_n(m_vTemp2SoundSpeed,m_iCapacity,0);
		
		if(m_iDimension==2 || m_iDimension==3) {
			
			m_vLPFOrderRight = new int[m_iCapacity];
			m_vLPFOrderLeft = new int[m_iCapacity];
			m_vLPFOrderNorth = new int[m_iCapacity];
			m_vLPFOrderSouth = new int[m_iCapacity]; 
			if(m_iDimension==3) m_vLPFOrderUp = new int[m_iCapacity];
			if(m_iDimension==3) m_vLPFOrderDown = new int[m_iCapacity];	
			
			fill_n(m_vLPFOrderRight,m_iCapacity,0);
			fill_n(m_vLPFOrderLeft,m_iCapacity,0);	
			fill_n(m_vLPFOrderNorth,m_iCapacity,0);
			fill_n(m_vLPFOrderSouth,m_iCapacity,0);
			if(m_iDimension==3) fill_n(m_vLPFOrderUp,m_iCapacity,0);
			if(m_iDimension==3) fill_n(m_vLPFOrderDown,m_iCapacity,0);
			
			if(init.getUseLimiter()) {
				m_vLPFFirstOrderRight = new int[m_iCapacity];
				m_vLPFFirstOrderLeft = new int[m_iCapacity];
				m_vLPFFirstOrderNorth = new int[m_iCapacity];
				m_vLPFFirstOrderSouth = new int[m_iCapacity]; 
				if(m_iDimension==3) m_vLPFFirstOrderUp = new int[m_iCapacity];
				if(m_iDimension==3) m_vLPFFirstOrderDown = new int[m_iCapacity];	
				
				fill_n(m_vLPFFirstOrderRight,m_iCapacity,0);
				fill_n(m_vLPFFirstOrderLeft,m_iCapacity,0);	
				fill_n(m_vLPFFirstOrderNorth,m_iCapacity,0);
				fill_n(m_vLPFFirstOrderSouth,m_iCapacity,0);
				if(m_iDimension==3) fill_n(m_vLPFFirstOrderUp,m_iCapacity,0);
				if(m_iDimension==3) fill_n(m_vLPFFirstOrderDown,m_iCapacity,0);
			}


			m_vNeighbourList = new int[m_iCapacity*m_iMaxNeighbourNum]; // whole
			m_vNeighbourListRight = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // right
			m_vNeighbourListLeft = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // left
			m_vNeighbourListNorth = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // north
			m_vNeighbourListSouth = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // south
			if(m_iDimension==3) m_vNeighbourListUp = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // up
			if(m_iDimension==3)	m_vNeighbourListDown = new int[m_iCapacity*m_iMaxNeighbourNumInOneDir]; // down		
			
			fill_n(m_vNeighbourList,m_iCapacity*m_iMaxNeighbourNum,0);
			fill_n(m_vNeighbourListRight,m_iCapacity*m_iMaxNeighbourNumInOneDir,0);
			fill_n(m_vNeighbourListLeft,m_iCapacity*m_iMaxNeighbourNumInOneDir,0);
			fill_n(m_vNeighbourListNorth,m_iCapacity*m_iMaxNeighbourNumInOneDir,0);
			fill_n(m_vNeighbourListSouth,m_iCapacity*m_iMaxNeighbourNumInOneDir,0);
			if(m_iDimension==3) fill_n(m_vNeighbourListUp,m_iCapacity*m_iMaxNeighbourNumInOneDir,0);
			if(m_iDimension==3) fill_n(m_vNeighbourListDown,m_iCapacity*m_iMaxNeighbourNumInOneDir,0);	
			
			m_vNeighbourListSize = new int [m_iCapacity];
			m_vNeighbourListRightSize = new int [m_iCapacity];
			m_vNeighbourListLeftSize = new int [m_iCapacity];
			m_vNeighbourListNorthSize = new int [m_iCapacity];
			m_vNeighbourListSouthSize = new int [m_iCapacity];
			if(m_iDimension==3) m_vNeighbourListUpSize = new int [m_iCapacity];
			if(m_iDimension==3) m_vNeighbourListDownSize = new int [m_iCapacity];
			
			fill_n(m_vNeighbourListSize,m_iCapacity,0);
			fill_n(m_vNeighbourListRightSize,m_iCapacity,0);
			fill_n(m_vNeighbourListLeftSize,m_iCapacity,0);
			fill_n(m_vNeighbourListNorthSize,m_iCapacity,0);
			fill_n(m_vNeighbourListSouthSize,m_iCapacity,0);
			if(m_iDimension==3) fill_n(m_vNeighbourListUpSize,m_iCapacity,0);
			if(m_iDimension==3) fill_n(m_vNeighbourListDownSize,m_iCapacity,0);
		}

		if(m_iDimension == 1) { // for 1D limiter
			m_vDD1 = new double[m_iCapacity];
			m_vDD2Left = new double[m_iCapacity];
			m_vDD2Right = new double[m_iCapacity];
			m_vDD3Left = new double[m_iCapacity];
			m_vDD3Right = new double[m_iCapacity];
			m_vCumP = new double[m_iCapacity];
			m_vPositionXm = new double[m_iCapacity];
			fill_n(m_vDD1,m_iCapacity,0);
			fill_n(m_vDD2Left,m_iCapacity,0);
			fill_n(m_vDD2Right,m_iCapacity,0);
			fill_n(m_vDD3Left,m_iCapacity,0);
			fill_n(m_vDD3Right,m_iCapacity,0);
			fill_n(m_vCumP,m_iCapacity,0);
			fill_n(m_vPositionXm,m_iCapacity,0);
		}

	}
	catch(bad_alloc& ba) {
		cerr<<"std::bad_alloc caught during initialization: "<<ba.what()<<endl;
		assert(false);
	}
}

void ParticleData::swap(size_t i, size_t j){
        double *x = m_vPositionX;
        double *y = m_vPositionY;
        double *z = m_vPositionZ;
        double *pressure = m_vPressure;
        double *velocityU = m_vVelocityU;
        double *velocityV = m_vVelocityV;
        double *velocityW = m_vVelocityW;
        double *volume = m_vVolume;
        double *volumeold = m_vVolumeOld;
        double *localParSpacing = m_vLocalParSpacing;
        double *mass = m_vMass;
        double *sound = m_vSoundSpeed;
	    double *leftintegral = m_vLeftIntegral;
	    double *rightintegral = m_vRightIntegral;
	    double *Deltaq = m_vDeltaq;
	    double *Qplusminus = m_vQplusminus;
        double *temperature = m_vTemperature;
        bool *LeftInflow = m_bLeftInflow;
	    int* pelletid = m_vPelletID;
        int* timetrack = m_vTimeTrack;
    
    std::swap(x[i],x[j]);
	std::swap(y[i],y[j]);
    std::swap(temperature[i],temperature[j]);
	if(m_iDimension == 3) std::swap(z[i],z[j]);
	std::swap(pressure[i],pressure[j]);
	std::swap(velocityU[i],velocityU[j]);
        std::swap(velocityV[i],velocityV[j]);
        if(m_iDimension == 3) std::swap(velocityW[i],velocityW[j]);
	std::swap(volume[i],volume[j]);
	std::swap(volumeold[i],volumeold[j]);
	std::swap(localParSpacing[i],localParSpacing[j]);
	std::swap(mass[i],mass[j]);
	std::swap(sound[i],sound[j]);
	std::swap(LeftInflow[i],LeftInflow[j]);
	std::swap(leftintegral[i],leftintegral[j]);
	std::swap(rightintegral[i],rightintegral[j]);
	std::swap(Deltaq[i],Deltaq[j]);
	std::swap(Qplusminus[i],Qplusminus[j]);
    std::swap(pelletid[i],pelletid[j]);
    std::swap(timetrack[i],timetrack[j]);
}

void ParticleData::makezero(size_t i){
	m_vPositionX[i]=0;
	m_vPositionY[i]=0;
	if(m_iDimension==3) m_vPositionZ[i]=0;
	m_vVelocityU[i]=0;
	m_vVelocityV[i]=0;
	if(m_iDimension==3) m_vVelocityW[i]=0;
	m_vPressure[i]=0;
	m_vVolume[i]=0;
	m_vVolumeOld[i]=0;
	m_vVolumeVoronoi[i]=0;
	m_vSoundSpeed[i]=0;
	m_vLocalParSpacing[i]=0;
	m_vMass[i]=0;
	m_vLeftIntegral[i]=0;
	m_vRightIntegral[i]=0;
	m_vDeltaq[i]=0;
	m_vQplusminus[i]=0;
	m_vVolume_x[i]=0;
	m_vVolume_y[i]=0;
	m_vVolume_z[i]=0;
	m_vDensity[i]=0;
	m_bLeftInflow[i]=false;
	m_vTemp1VelocityU[i]=0;
	m_vTemp1VelocityV[i]=0;
	if(m_iDimension==3)
		m_vTemp1VelocityW[i]=0;
	m_vTemp1Pressure[i]=0;
	m_vTemp1Volume[i]=0;
	m_vTemp1SoundSpeed[i]=0;
	m_vTemp1PositionX[i]=0;	
        m_vTemp2VelocityU[i]=0;
        m_vTemp2VelocityV[i]=0;
        if(m_iDimension==3)
                m_vTemp2VelocityW[i]=0;
        m_vTemp2Pressure[i]=0;
        m_vTemp2Volume[i]=0;
        m_vTemp2SoundSpeed[i]=0;
	m_vLPFOrderRight[i]=1;
	m_vLPFOrderLeft[i]=1;
	m_vLPFOrderNorth[i]=1;
	m_vLPFOrderSouth[i]=1;
	if(m_iDimension==3)
	{
		m_vLPFOrderUp[i]=1;
		m_vLPFOrderDown[i]=1;
	}
	m_vNeighbourListSize[i]=0;
	m_vNeighbourListRightSize[i]=0;
        m_vNeighbourListLeftSize[i]=0;
        m_vNeighbourListNorthSize[i]=0;
        m_vNeighbourListSouthSize[i]=0;
	if(m_iDimension==3)
	{
	        m_vNeighbourListUpSize[i]=0;
	        m_vNeighbourListUpSize[i]=0;
	}
	m_vPelletID[i]=-1;

}
void ParticleData::augmentAllDataArrays(size_t newCapacity) {
	
	//size_t oldCapacity = m_iCapacity;
	//cout<<"Capacity change from "<<m_iCapacity;
	m_iCapacity = newCapacity;
	//cout<<" to "<<m_iCapacity<<endl;
	
	// keep the old data
	
	augment<double>(m_vPositionX, newCapacity, true);
	augment<double>(m_vPositionY, newCapacity, true);
	augment<double>(m_vPositionZ, newCapacity, true);	

	augment<double>(m_vVelocityU, newCapacity, true);
	if(m_iDimension==2 || m_iDimension==3) 
		augment<double>(m_vVelocityV, newCapacity, true);
	if(m_iDimension==3) 
		augment<double>(m_vVelocityW, newCapacity, true);

	augment<double>(m_vPressure, newCapacity, true);
	augment<double>(m_vVolume, newCapacity, true);
	augment<double>(m_vVolumeOld, newCapacity, true);
	augment<double>(m_vVolumeVoronoi, newCapacity, true);
	augment<double>(m_vSoundSpeed, newCapacity, true);
	augment<double>(m_vLocalParSpacing, newCapacity, true);
        augment<double>(m_vMass, newCapacity, true);
        augment<double>(m_vLeftIntegral, newCapacity, true);
        augment<double>(m_vRightIntegral, newCapacity, true);
	augment<double>(m_vDeltaq, newCapacity, true);
	augment<double>(m_vQplusminus, newCapacity, true);
        augment<double>(m_vVolume_x, newCapacity, true);
        augment<double>(m_vVolume_y, newCapacity, true);
        augment<double>(m_vVolume_z, newCapacity, true);
	augment<double>(m_vDensity, newCapacity, true);
        augment<bool>(m_bLeftInflow, newCapacity, true);
	augment<int>(m_vPelletID, newCapacity, true);	
	if(m_iDimension==2 || m_iDimension==3)
		augment<int>(m_vObjectTag, newCapacity, true);
	
	
	// do not need to keep old data
	if(m_iDimension==2 || m_iDimension==3) {
		// Temp1
#ifdef LW_DEBUG
                augment<double>(m_vPhi, newCapacity, false);
                augment<double>(m_vPError0, newCapacity, false);
                augment<double>(m_vPError1, newCapacity, false);
                augment<double>(m_vVelError0, newCapacity, false);
                augment<double>(m_vVelError1, newCapacity, false);
                augment<double>(m_vPxl, newCapacity, false);
                augment<double>(m_vPxr, newCapacity, false);
                augment<double>(m_vVxl, newCapacity, false);
                augment<double>(m_vVxr, newCapacity, false);
                augment<double>(m_vPyl, newCapacity, false);
                augment<double>(m_vPyr, newCapacity, false);
                augment<double>(m_vVyl, newCapacity, false);
                augment<double>(m_vVyr, newCapacity, false);
                augment<double>(m_vVtx, newCapacity, false);
                augment<double>(m_vVty, newCapacity, false);
                augment<double>(m_vPtx, newCapacity, false);
                augment<double>(m_vPty, newCapacity, false);
                augment<double>(m_vVolumetx, newCapacity, false);
                augment<double>(m_vVolumety, newCapacity, false);
		augment<int>(m_vNeighSize, 4*newCapacity, false);
		augment<int>(m_vNeighList, 40*newCapacity, false);
                augment<int>(m_vNeighOfParticle, newCapacity, false);
		augment<int>(m_vIfSPHDensity, newCapacity, true);
#endif
		augment<double>(m_vTemp1VelocityU, newCapacity, false);
		augment<double>(m_vTemp1VelocityV, newCapacity, false);
		if(m_iDimension==3) 
			augment<double>(m_vTemp1VelocityW, newCapacity, false);
		augment<double>(m_vTemp1Pressure, newCapacity, false);
		augment<double>(m_vTemp1Volume, newCapacity, false);
		augment<double>(m_vTemp1SoundSpeed, newCapacity, false);
		augment<double>(m_vTemp1PositionX, newCapacity, false);
		
		// Temp2
		augment<double>(m_vTemp2VelocityU, newCapacity, false);
		augment<double>(m_vTemp2VelocityV, newCapacity, false);
		if(m_iDimension==3) 
			augment<double>(m_vTemp2VelocityW, newCapacity, false);
		augment<double>(m_vTemp2Pressure, newCapacity, false);
		augment<double>(m_vTemp2Volume, newCapacity, false);
		augment<double>(m_vTemp2SoundSpeed, newCapacity, false);	
		
		augment<int>(m_vLPFOrderRight, newCapacity, false);
		augment<int>(m_vLPFOrderLeft, newCapacity, false);
		augment<int>(m_vLPFOrderNorth, newCapacity, false);
		augment<int>(m_vLPFOrderSouth, newCapacity, false);
		if(m_iDimension==3) {
			augment<int>(m_vLPFOrderUp, newCapacity, false);
			augment<int>(m_vLPFOrderDown, newCapacity, false);
		}
		
		augment<int>(m_vNeighbourListSize, newCapacity, false);
		augment<int>(m_vNeighbourListRightSize, newCapacity, false);
		augment<int>(m_vNeighbourListLeftSize, newCapacity, false);
		augment<int>(m_vNeighbourListNorthSize, newCapacity, false);
		augment<int>(m_vNeighbourListSouthSize, newCapacity, false);
		if(m_iDimension==3) {
			augment<int>(m_vNeighbourListUpSize, newCapacity, false);
			augment<int>(m_vNeighbourListDownSize, newCapacity, false);
		}
		
		augment<int>(m_vNeighbourList, newCapacity*m_iMaxNeighbourNum, false);
	
		augment<int>(m_vNeighbourListRight, newCapacity*m_iMaxNeighbourNumInOneDir, false);
	
		augment<int>(m_vNeighbourListLeft, newCapacity*m_iMaxNeighbourNumInOneDir, false);

		augment<int>(m_vNeighbourListNorth, newCapacity*m_iMaxNeighbourNumInOneDir, false);
		
		augment<int>(m_vNeighbourListSouth, newCapacity*m_iMaxNeighbourNumInOneDir, false);
		
		if(m_iDimension==3) {
			augment<int>(m_vNeighbourListUp, newCapacity*m_iMaxNeighbourNumInOneDir, false);
			augment<int>(m_vNeighbourListDown, newCapacity*m_iMaxNeighbourNumInOneDir, false);
		}


	}	


}

ParticleData::~ParticleData() {
  
	delete[] m_vPositionX;
	delete[] m_vPositionY;
	delete[] m_vPositionZ;

	delete[] m_vVelocityU;
	if(m_iDimension==2 || m_iDimension==3) delete[] m_vVelocityV;
	if(m_iDimension==3) delete[] m_vVelocityW;

	delete[] m_vVolume;
	delete[] m_vVolumeOld;
	delete[] m_vVolumeVoronoi;
	delete[] m_vPressure;
	delete[] m_vTemperature;
    delete[] m_vTimeTrack;
    //delete[] m_vEnergy;
	delete[] m_vSoundSpeed;
	delete[] m_vLocalParSpacing;

	delete[] m_vObjectTag;

	delete[] m_vTemp1VelocityU;
	if(m_iDimension==2 || m_iDimension==3) delete[] m_vTemp1VelocityV;
	if(m_iDimension==3) delete[] m_vTemp1VelocityW;
	delete[] m_vTemp1Volume;
	delete[] m_vTemp1Pressure;
	delete[] m_vTemp1SoundSpeed;
	delete[] m_vTemp1PositionX;

	delete[] m_vTemp2VelocityU;
	if(m_iDimension==2 || m_iDimension==3) delete[] m_vTemp2VelocityV;
	if(m_iDimension==3) delete[] m_vTemp2VelocityW;
	delete[] m_vTemp2Volume;
	delete[] m_vTemp2Pressure;
	delete[] m_vTemp2SoundSpeed;
	delete[] m_vMass;
	
    if(m_iNumberofPellet){
    delete[] m_vLeftIntegral;
	delete[] m_vRightIntegral;
	delete[] m_vDeltaq;
	delete[] m_vQplusminus;
    delete[] m_vPelletID;
	delete[] m_vPelletPositionX;
	delete[] m_vPelletPositionY;
	delete[] m_vPelletPositionZ;
	delete[] m_vPelletRadius;
	delete[] m_vPelletInnerRadius;
	delete[] m_vPelletOuterRadius;
    delete[] m_vPelletEnergy;
	delete[] m_vPelletVelocity;
    delete[] m_vMassFlowRate;
    delete[] m_vPelletState;
    }
	delete[] m_vVolume_x;
    delete[] m_vVolume_y;
    delete[] m_vVolume_z;
	delete[] m_vDensity;
	delete[] m_bLeftInflow;
	
	delete[] m_vPhi;
    
	delete[] m_vIfSPHDensity;
#ifdef LW_DEBUG
	    delete[] m_vPError0;
        delete[] m_vPError1;
        delete[] m_vVelError0;
        delete[] m_vVelError1;
	    delete[] m_vPxl;
        delete[] m_vPxr;
        delete[] m_vVxl;
        delete[] m_vVxr;
        delete[] m_vPyl;
        delete[] m_vPyr;
        delete[] m_vVyl;
        delete[] m_vVyr;
        delete[] m_vVtx;
        delete[] m_vVty;
        delete[] m_vPtx;
        delete[] m_vPty;
        delete[] m_vVolumetx;
        delete[] m_vVolumety;
	delete[] m_vNeighList;
	delete[] m_vNeighSize;
	delete[] m_vNeighOfParticle;
#endif	
	if(m_iDimension==2 || m_iDimension==3) {

		delete[] m_vLPFOrderRight;
		delete[] m_vLPFOrderLeft;
		delete[] m_vLPFOrderNorth;
		delete[] m_vLPFOrderSouth;
		if(m_iDimension==3) delete[] m_vLPFOrderUp;
		if(m_iDimension==3) delete[] m_vLPFOrderDown;
		
		
		delete[] m_vNeighbourList;
		delete[] m_vNeighbourListRight;
		delete[] m_vNeighbourListLeft;
		delete[] m_vNeighbourListNorth;
		delete[] m_vNeighbourListSouth;	
		if(m_iDimension==3) delete[] m_vNeighbourListUp;	
		if(m_iDimension==3)	delete[] m_vNeighbourListDown;
		

		delete[] m_vNeighbourListSize;
		delete[] m_vNeighbourListRightSize;
		delete[] m_vNeighbourListLeftSize;
		delete[] m_vNeighbourListNorthSize;
		delete[] m_vNeighbourListSouthSize;
		delete[] m_vNeighbourListUpSize;
		delete[] m_vNeighbourListDownSize;
	
	}	

	if(m_iDimension == 1) {
		delete[] m_vDD1;
		delete[] m_vDD2Left;
		delete[] m_vDD2Right;
		delete[] m_vDD3Left;
		delete[] m_vDD3Right;
		delete[] m_vCumP;
		delete[] m_vPositionXm;
	}
}


////////////////////////////////////////////////////////////////////////////////
// End : ParticleData
////////////////////////////////////////////////////////////////////////////////
