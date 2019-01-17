#ifndef __PELLET_SOLVER__
#define __PELLET_SOLVER__

#include "lp_solver.h"
#include "particle_data.h"
#include "initializer.h"
#include <algorithm>
#include <math.h>
class PelletSolver {

    public: 
       PelletSolver(const Initializer& init ,ParticleData* pdata);
        ~PelletSolver(){};
        
	void calculateHeatDeposition( double dt);
	void computeIntegralSpherical();
	void updateStatesByLorentzForce(double dt);

 	double neon_radiation_power_density(double,double);

	double neon_radiation_data(int,int);
       
    private:

        
        EOS* m_pEOS;
         
        ParticleData* m_pPelletData;
};























#endif
