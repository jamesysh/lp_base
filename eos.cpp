#include "eos.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstring>
#include <unistd.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
using namespace std;
#define INTERP_TYPE gsl_interp2d_bicubic // bicubic or bilinear

static const gsl_interp2d_type *T_sound_speed = INTERP_TYPE;
static gsl_spline2d *spline_sound_speed;
static const gsl_interp2d_type *T_temperature = INTERP_TYPE;
static gsl_spline2d *spline_temperature;
static const gsl_interp2d_type *T_conductivity = INTERP_TYPE;
static gsl_spline2d *spline_conductivity;


////////////////////////////////////////////////////////////////////////////////
// Start of PolytropicGasEOS
////////////////////////////////////////////////////////////////////////////////

double PolytropicGasEOS::getEnergy(double pressure, double density) {
	if(((m_fGamma - 1.) * density) != 0) 
		return ( pressure / ((m_fGamma - 1.) * density) );
	else { // divide by zero
		std::cout<<"Error (Divide by zero)! Computing energy by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", density = "<<density<<std::endl;
		assert(false);
	}
}

double PolytropicGasEOS::getSoundSpeed(double pressure, double density) {
	double cs;
	if(density != 0)
		cs = m_fGamma * pressure / density;
	else {
		std::cout<<"Error (Divide by zero density)! Computing sound speed by EOS: "<<std::endl;
		//std::cout<<"density = "<<density<<std::endl;
		assert(false);
	}
	if(cs > 0) 
		return sqrt(cs);
	else if(cs==0)
		return cs;
	else { // taking square root of a negative number
		std::cout<<"Error (Taking suqre root of a negative number)! Computing sound speed by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", pressure = "<<pressure<<", density = "<<density<<std::endl;
		assert(false);
	}
}

double PolytropicGasEOS::getElectricConductivity(double pressure, double density) {
if(m_iPelletMaterial == 0)
    return 0.0;
else if(m_iPelletMaterial == 1)
{   
    double cond;
    double T = getTemperature(pressure,density);
    double m = 1 ;
    double nt,fi,lnL,f_one;   
    nt =  density/3.34e-24;
    fi = m;
    if(fi < 1.e-20)
        cond = 0.;
    else{
        lnL = 3.6e9*pow(T,1.5)/sqrt(fi*nt);
        cond = 1e11/1.15/(lnL/pow(T,1.5) + 0.54*pow(T,0.059)*(1/fi-1));

        
        }
      
    return cond;
    }
}

double PolytropicGasEOS::getTemperature(double pressure, double density) {
  double R,mu;
  
  if(m_iPelletMaterial == 0 ){
    mu = 20.18;
    R = 83.14;
  }
  else if(m_iPelletMaterial == 1){
       R = 83.14;
       mu = 2.014;
      }
  return mu*pressure/(R*density)/11604.525;
}

////////////////////////////////////////////////////////////////////////////////
// End of PolytropicGasEOS
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
// Start : StiffPolytropicEOS
////////////////////////////////////////////////////////////////////////////////

double StiffPolytropicGasEOS::getEnergy(double pressure, double density) {
	if(((m_fGamma - 1.) * density) != 0) 
		return ( (pressure + m_fGamma * m_fPinf) / ((m_fGamma - 1.) * density) - m_fEinf);
	else {
		std::cout<<"Error (Divide by zero)! Computing energy by EOS: "<<std::endl;
		std::cout<<"gamma = "<<m_fGamma<<", density = "<<density<<std::endl;
		assert(false);
	}
  
}


double StiffPolytropicGasEOS::getSoundSpeed(double pressure, double density) {
	double cs;
	if(density != 0)
		cs = m_fGamma * (pressure + m_fPinf) / density;
	else {
		std::cout<<"Error (Divide by zero density)! Computing sound speed by EOS: "<<std::endl;
		std::cout<<"density = "<<density<<std::endl;
		assert(false);
	}
	if(cs > 0)
		return sqrt(cs);
	else if(cs==0)
		return cs;
	else { 
		std::cout<<"Error (Taking suqre root of a negative number)! Computing sound speed by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", pinf = "<<m_fPinf<<", pressure = "<<pressure<<", density = "<<density<<std::endl;
		assert(false);
	}  
}

double  StiffPolytropicGasEOS::getElectricConductivity(double pressure, double density) {
  return 0.0;
}

double  StiffPolytropicGasEOS::getTemperature(double pressure, double density) {
  double mu = 20.18;
  double R = 83.14;
  
  return mu*(pressure + m_fPinf)/(R*density);
}

////////////////////////////////////////////////////////////////////////////////
// End of StiffPolytropicEOS
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
// Start : SahaEOS
///////////////////////////////////////////////////////////////////////////////


void  SahaEOS::getParameters(std::vector<double>& params)
{
  FILE *fp_rho, *fp_pres, *fp_sound_speed, *fp_temperature, *fp_conductivity;
  int dim = 1000;
  double* rho = new double[dim];
  double* pres = new double[dim];
  double* sound_speed = new double[dim*dim];
  double* temperature = new double[dim*dim];
  double* conductivity = new double[dim*dim];
  double rho_entry, pres_entry, sc_entry, t_entry, cond_entry;
  int i,j;

  printf("ENTERED  SahaEOS::getParameters \n");
  char* buffer;
  buffer = getcwd(NULL,0);
  string sub_dir = "/tables_1000/";
  string main_dir(buffer);
  string dir = main_dir+sub_dir;
  
  if (!(fp_rho = fopen((dir+"rho.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS rho.dat file \n");

  if (!(fp_pres = fopen((dir+"pres.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS pres.dat file \n");

  if (!(fp_sound_speed = fopen((dir+"sound_speed.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS sound_speed.dat file \n");

  if (!(fp_temperature = fopen((dir+"temp.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS temp.dat file \n");

  if (!(fp_conductivity = fopen((dir+"conductivity.dat").c_str(),"r")))
    printf("CANNOT OPEN Saha EOS conductivity.dat file \n");

  printf("Opend files in  SahaEOS::getParameters \n");
 
  i = 0;
  while (fscanf(fp_rho,"%lf",&rho_entry)!=EOF)
    rho[i++] = rho_entry;
  i = 0;
  while (fscanf(fp_pres,"%lf",&pres_entry)!=EOF)
    pres[i++] = pres_entry;

  spline_sound_speed = gsl_spline2d_alloc(T_sound_speed, dim, dim);
  spline_temperature = gsl_spline2d_alloc(T_temperature, dim, dim);
  spline_conductivity = gsl_spline2d_alloc(T_conductivity, dim, dim);

  printf("BUILT 2D splies in  SahaEOS::getParameters \n");

  for (i=0; i<dim; ++i)
    for (j=0; j<dim; ++j)
    {
      fscanf (fp_sound_speed,"%lf",&sc_entry);
      gsl_spline2d_set(spline_sound_speed, sound_speed, i,j, sc_entry);

      fscanf (fp_temperature,"%lf",&t_entry);
      gsl_spline2d_set(spline_temperature, temperature, i,j, t_entry);
      
      fscanf (fp_conductivity,"%lf",&cond_entry);
      gsl_spline2d_set(spline_conductivity, conductivity, i,j, cond_entry);           
    }
  
  fclose(fp_rho);
  fclose(fp_pres);
  fclose(fp_sound_speed);
  fclose(fp_temperature);
  fclose(fp_conductivity);
  
  gsl_spline2d_init(spline_sound_speed, rho, pres, sound_speed, dim, dim);
  gsl_spline2d_init(spline_temperature, rho, pres, temperature, dim, dim);
  gsl_spline2d_init(spline_conductivity, rho, pres, conductivity, dim, dim);
  delete[] rho;
  delete[] pres;
  delete[] sound_speed;
  delete[] temperature;
  delete[] conductivity;
}



double SahaEOS::getEnergy(double pressure, double density) {
  double energy = 0.0;

  return energy;
}


double SahaEOS::getSoundSpeed(double pressure, double density) {
  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double sc;

  if (pressure < 1.e-10 || pressure > 200.0 || density < 1.e-12 || density > 2)
    sc = sqrt(1.67*pressure/density);
  else
    sc =  gsl_spline2d_eval(spline_sound_speed, density, pressure, densacc, presacc);

  //  printf("rho = %e  P = %e  sc = %e\n",density,pressure,sc);

  return sc;
}

double SahaEOS::getElectricConductivity(double pressure, double density) {
  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double cond;

   if (pressure < 1.e-10 || pressure > 200.0 || density < 1.e-12 || density > 2)
      cond = 0.0;
    else
    cond =  gsl_spline2d_eval(spline_conductivity, density, pressure, densacc, presacc);
 
  return cond;
}

double SahaEOS::getTemperature(double pressure, double density) {
  static gsl_interp_accel *densacc = gsl_interp_accel_alloc();
  static gsl_interp_accel *presacc = gsl_interp_accel_alloc();
  double T;
  double mu = 20.18;
  double R = 83.14;
  
   if (pressure < 1.e-10 || pressure > 200.0 || density < 1.e-12 || density > 2)
      T = mu*pressure/(R*density)/11604.525;
    else
    T =  gsl_spline2d_eval(spline_temperature, density, pressure, densacc, presacc);
 
  return T;
}
void SahaEOS::diagnosis(double rho0, double rho1, double p0, double p1){
    
    int rho_n = 5; 
    double rho_offset = (rho1-rho0)/(rho_n-1);
    int p_n = 20000;
    double p_offset = (p1-p0)/(p_n-1);
    for(int rho_i = 0; rho_i<rho_n; rho_i++){
          string PT = "pres_temp_"+to_string(rho0+rho_i*rho_offset)+".txt";
           FILE * file1 = fopen(PT.c_str(),"w");
  
          string CT = "cond_temp_"+to_string(rho0+rho_i*rho_offset)+".txt";
          FILE* file2 = fopen(CT.c_str(),"w");

        for(int p_i = 0; p_i<p_n; p_i++){

            double temp = getTemperature(p0+p_i*p_offset,rho0+rho_i*rho_offset); 
            fprintf(file1, "%.16g %.16g\n", p0+p_i*p_offset,temp);
            double conductivity = getElectricConductivity(p0+p_i*p_offset,rho0+rho_i*rho_offset);
            fprintf(file2, "%.16g %.16g\n",conductivity,temp);
            }
        fclose(file1);
        fclose(file2);
        
        }

    
    }


////////////////////////////////////////////////////////////////////////////////
// End of SahaEOS
////////////////////////////////////////////////////////////////////////////////
