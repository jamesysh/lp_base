/**
 * \file   eos.h
 *
 * \brief  This header file contains classes for the calculation of the Equation of State (EOS)
 *
 */


#ifndef __EOS_H__
#define __EOS_H__

#include <vector>

class EOS {
protected:
	int m_iEOSChoice; ///< The eos choice: 1=Polytropic gas; 2=Stiffened Polytropic gas; 3=Saha Eos
public:
	/// Destructor
	virtual ~EOS() {};

	/// Getter function of the protected data member m_iEOSChoice
	int getEOSChoice() {return m_iEOSChoice;}
	
	virtual void getParameters(std::vector<double>& params) = 0;
	virtual double getEnergy(double pressure, double density) = 0;
 	virtual double getTemperature(double pressure, double density) = 0;
      	virtual double getSoundSpeed(double pressure, double density) = 0;
	virtual double getElectricConductivity(double pressure, double density) = 0;
    virtual void diagnosis(double rho0, double rho1, double t0, double t1) = 0;
};



class PolytropicGasEOS : public EOS {
protected:
	double m_fGamma; ///< The parameter \e gamma
    
    int m_iPelletMaterial;
public:
	PolytropicGasEOS(double gamma,int pelletmaterial) : m_fGamma(gamma),m_iPelletMaterial(pelletmaterial) {m_iEOSChoice=1;}
	
	// Destructor
	virtual ~PolytropicGasEOS() {}	
	virtual void getParameters(std::vector<double>& params) {params.push_back(m_fGamma);}
	virtual double getEnergy(double pressure, double density);
	virtual double getTemperature(double pressure, double density);
	virtual double getSoundSpeed(double pressure, double density);
	virtual double getElectricConductivity(double pressure, double density);
    virtual void diagnosis(double rho0, double rho1, double t0, double t1){};

};


class StiffPolytropicGasEOS : public EOS {
protected:
	double m_fGamma; ///< The parameter \e gamma
	double m_fPinf; ///< The parameter pressure infinity 
	double m_fEinf; ///< The parameter energy infinity

public:
	/// Constructor
        StiffPolytropicGasEOS(double gamma, double pinf, double einf):
		m_fGamma(gamma), m_fPinf(pinf), m_fEinf(einf) {m_iEOSChoice=2;}
	
	/// Destructor
	virtual ~StiffPolytropicGasEOS() {}
	
	virtual void getParameters(std::vector<double>& params) { 
		params.push_back(m_fGamma);
		params.push_back(m_fPinf);
		params.push_back(m_fEinf);
	}

	virtual double getEnergy(double pressure, double density);
	virtual double getTemperature(double pressure, double density);
	virtual double getSoundSpeed(double pressure, double density);
	virtual double getElectricConductivity(double pressure, double density);
     virtual void diagnosis(double rho0, double rho1, double t0, double t1){};

};



class SahaEOS : public EOS {
protected:
	double m_fGamma; ///< The parameter \e gamma

public:
	/// Constructor
	SahaEOS(double gamma) : m_fGamma(gamma) {m_iEOSChoice=3;}
	
	// Destructor
	virtual ~SahaEOS() {}	
	virtual void getParameters(std::vector<double>& params);

	virtual double getEnergy(double pressure, double density);
	virtual double getTemperature(double pressure, double density);
	virtual double getSoundSpeed(double pressure, double density);
	virtual double getElectricConductivity(double pressure, double density);
    virtual void diagnosis(double rho0, double rho1, double p0, double p1);

};



#endif // __EOS_H__
