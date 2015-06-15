#ifndef _ENERGY
#define	_ENERGY

#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>

#include "mocha_numbers.h"
#include "mocha_types.h"


namespace mocha_energy_nmspc
{
  
  
  using namespace std;
  using namespace votca::tools;
  using namespace mocha_types_nmspc;
  using namespace mocha_numbers_nmspc;
  
  
  
  
  void energy_lj(mocha_distance d, mocha_lj_params lj_params, double *energy, double *virial)
  {
    
    mocha_distance s;
    double s_rep, s_att, erg, vir;
    
    erg = 0.0;
    vir = 0.0;

    if(lj_params.epsilon > 0.0)
    {
        if((d.r <= lj_params.cutoff))
        {
            s.r = d.r / lj_params.sigma;
            s_rep = s.invpow(lj_params.Erep);
            s_att = s.invpow(lj_params.Eatt);
            
            erg = 4.0 * lj_params.epsilon * (s_rep - s_att) - lj_params.shift;
            if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
                vir = - 4.0 * lj_params.epsilon * (lj_params.Erep * s_rep - lj_params.Eatt * s_att);
            else
                vir = 0.0;
        }
    }

    *energy = erg;
    *virial = vir;
    
  }



  
  void energy_lj_CG(mocha_distance d, mocha_lj_params lj_params, double *energy, double *virial)
  {

    mocha_distance s;
    double rshifted;
    //double s_rep, s_att, erg, vir;
    double s_rep, erg, vir;

    erg = 0.0;
    vir = 0.0;

    if(lj_params.epsilon > 0.0)
    {
        rshifted = d.r - 0.005;
        if((rshifted <= lj_params.cutoff))
        {
            s.r = rshifted / lj_params.sigma;
            s_rep = s.invpow(lj_params.Erep);
            //_att = s.invpow(lj_params.Eatt);

            //erg = 4.0 * lj_params.epsilon * (s_rep - s_att) - lj_params.shift;
            erg = 4.0 * lj_params.epsilon * ((s_rep) - 0.25);
            if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
            
                vir = - 4.0 * lj_params.epsilon * (lj_params.Erep * s_rep) * d.r / rshifted;
            else
                vir = 0.0;
        }
    }

    *energy = erg;
    *virial = vir;

  }

  

  //SHIFTED NOW!!! SHIFT HARDCODED!!!
  void energy_WCA(mocha_distance d, mocha_lj_params lj_params, double *energy, double *virial)
  {
    //cout << "WCA CALLED!\n";
    mocha_distance s;
    double s_rep, s_att, erg, vir, shift;
   
    erg = 0.0;
    vir = 0.0;
    shift = 0.147;  //SHIFTED NOW!!! SHIFT HARDCODED!!!

    if(lj_params.epsilon > 0.0)
    {
        if((d.r <= lj_params.cutoff + shift))
        {
            s.r = (d.r - shift) / lj_params.sigma;
            s_rep = s.invpow(lj_params.Erep);
            s_att = s.invpow(lj_params.Eatt);
            
            erg = 4.0 * lj_params.epsilon * (s_rep - s_att) + lj_params.epsilon;
            if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
                vir = - 4.0 * lj_params.epsilon * (d.r/(d.r - shift)) * (lj_params.Erep * s_rep - lj_params.Eatt * s_att);
            else
                vir = 0.0;
        }
    }

    *energy = erg;
    *virial = vir;
    
  }
  
  //IBI FITTED LJ POTENTIAL!!!
  void energy_IBI_LJ(mocha_distance d, mocha_lj_params lj_params, double beta, double *energy, double *virial)
  {
  // HARDCODED...      
      
    mocha_distance s;
    double s_rep, s_att, erg, vir, eps, sig, cuto, shifted;
    
    erg = 0.0;
    vir = 0.0;
    
    eps = 0.17;
    sig = 0.3298;
    cuto = 0.9;//pow(2.0, 1.0/6.0) * sigma;
    shifted = -0.0016425;//0.147;

    if(lj_params.epsilon > 0.0)
    {
        if((d.r <= cuto))
        {
            s.r = d.r / sig;
            s_rep = s.invpow(lj_params.Erep);
            s_att = s.invpow(lj_params.Eatt);
            
            erg = 4.0 * eps * beta * (s_rep - s_att) - shifted * beta;
            if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
                vir = - 4.0 * eps * beta * (lj_params.Erep * s_rep - lj_params.Eatt * s_att);
            else
                vir = 0.0;
        }
    }

    *energy = erg;
    *virial = vir;
    
  }
  
  
  
  
  void energy_SG(mocha_distance d, mocha_lj_params lj_params, double beta, double *energy, double *virial)
  {
    
    mocha_distance s;
    double erg, vir, C6, C8, C9, C10, alp, bet, gam, rm, e_const, l_const, e_shift;
    
    C6 = 12.14;
    C8 = 215.2;
    C9 = 143.1;
    C10 = 4813.9;
    alp = 1.713;
    bet = 1.5671;
    gam = 0.00993;
    rm = 8.32;
    e_const = 2625.497656;
    l_const = 0.052917721092;
    //e_shift = 0.000735598179617464; (cutoff at 1nm)
    e_shift = 0.001401719157035950; //(cutoff at 0.9nm)
    
    erg = 0.0;
    vir = 0.0;

    s.r = d.r / l_const;
    
    if(lj_params.epsilon > 0.0)
    {
        if(d.r <= 0.9)
        {          

            if(d.r <= 0.44027543948544)
            {
                erg = e_const*beta*( exp(alp - bet*(s.r) - gam*(s.pow(2))) - (C6*(s.invpow(6)) + C8*(s.invpow(8)) - C9*(s.invpow(9)) + C10*(s.invpow(10))) * exp(-(rm/(s.r)-1.0)*(rm/(s.r)-1.0)) ) + e_shift*beta;

                if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
                    vir = e_const*beta*( (-bet*s.r-2.0*gam*s.pow(2))*exp(alp - bet*(s.r) - gam*(s.pow(2))) + (6.0*C6*(s.invpow(6)) + 8.0*C8*(s.invpow(8)) - 9.0*C9*(s.invpow(9)) + 10.0*C10*(s.invpow(10))) * exp(-(rm/(s.r)-1)*(rm/(s.r)-1)) - (C6*(s.invpow(6)) + C8*(s.invpow(8)) - C9*(s.invpow(9)) + C10*(s.invpow(10))) * exp(-(rm/(s.r)-1.0)*(rm/(s.r)-1.0)) * 2.0 * (rm/s.r -1.0) *rm/s.r );          
            }
            else
            {
                erg = e_const*beta*( exp(alp - bet*(s.r) - gam*(s.pow(2))) - (C6*(s.invpow(6)) + C8*(s.invpow(8)) - C9*(s.invpow(9)) + C10*(s.invpow(10))) ) + e_shift*beta;

                if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
                    vir = e_const*beta*( (-bet*s.r-2.0*gam*s.pow(2))*exp(alp - bet*(s.r) - gam*(s.pow(2))) + (6.0*C6*(s.invpow(6)) + 8.0*C8*(s.invpow(8)) - 9.0*C9*(s.invpow(9)) + 10.0*C10*(s.invpow(10))) );              
            }

        }    
    }
    *energy = erg;
    *virial = vir;
    
  }
  
  
  

  void energy_fene(mocha_distance d, double rmax, double kappa, double *energy, double *virial)
  {
    
    double erg, vir;
    mocha_distance temp;
    
    erg = 0.0;
    vir = 0.0;

    if(kappa > 0.0)
    {
      if(d.r >= 0.999999*rmax)
          temp.r = 0.999999;
      else
          temp.r = d.r / rmax;

      erg = - 0.5 * kappa * rmax * rmax * log( 1.0 -  temp.pow(2) );
      if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
	vir = kappa * temp.pow(2) / ( 1.0 -  temp.pow(2) );
      else
          vir = 0.0;

    }
    
    *energy = erg;
    *virial = vir;
    
  }
  
  
  
  
  
  






  void energy_quartic(mocha_distance d, double rmax, double kappa, double *energy, double *virial)
  {

    double erg, vir;

    erg = 0.0;
    vir = 0.0;

    if(kappa > 0.0)
    {

        erg = 0.25 * kappa * (d.pow(2) - rmax*rmax) * (d.pow(2) - rmax*rmax);

      if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
          vir = kappa * (d.pow(2) - rmax*rmax) * d.pow(2);
      else
          vir = 0.0;

    }

    *energy = erg;
    *virial = vir;

  }
  
  
  void energy_quadratic(mocha_distance d, double rmax, double kappa, double *energy, double *virial)
  {

    double erg, vir;

    erg = 0.0;
    vir = 0.0;

    if(kappa > 0.0)
    {

        erg = 0.5 * kappa * (d.pow(2));

      if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
          vir = kappa * (d.pow(2));
      else
          vir = 0.0;

    }

    *energy = erg;
    *virial = vir;

  }
  
  
  void energy_electric(mocha_distance d, double q1, double q2, double cutoff, double inv_cutoff, double beta, double *energy, double *virial)
  {

    double erg, vir;

    erg = 0.0;
    vir = 0.0;

    if(d.r <= cutoff)
    {

        erg = beta*q1*q2*DIELECTRIC*(d.invpow(1) - inv_cutoff);

      if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
          vir = -erg;
      else
          vir = 0.0;

    }

    *energy = erg;
    *virial = vir;

  }
  
}



#endif






