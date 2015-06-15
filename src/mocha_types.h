#ifndef _TYPES
#define	_TYPES

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include <omp.h>

#include "vec.h"


namespace mocha_types_nmspc
{
  
  using namespace std;
  using namespace votca::tools;
  
  
  class mocha_atom
  {
    
    public:
      double mass;
      double charge;
      
      double lambda;
      //double new_lambda;
      
      vec position;
      vec position_unfolded;      
      vec dist_trav;      
      string type;
      int interaction_number;
      void copy(mocha_atom Atom)
      {
	position = Atom.position;
	position_unfolded = Atom.position_unfolded;
	type = Atom.type;
	mass = Atom.mass;
        interaction_number = Atom.interaction_number;
        charge = Atom.charge;
      }
      void initialize()
      {
	mass = 0.0;
	position = vec(0.0, 0.0, 0.0);
	position_unfolded = vec(0.0, 0.0, 0.0);
        dist_trav = vec(0.0, 0.0, 0.0);
	type = "";
        interaction_number = -1;
        charge = 0.0;
        
        lambda = 1.0;
        //new_lambda = 1.0;
      }
      
  };
  
  
  
  struct pairlist
  {
    int a, b;
    double k;
    double rzero;
    double rmax;
  };
  
  
  struct bonded_molecule_pair
  {
    int a, b;
    int Nbonds;
    struct pairlist *pairs;
  };
  
  
  class mocha_lj_params
  {
    
    public:
      string type;
      int Erep;
      int Eatt;
      int label;
      double sigma;
      double epsilon;
      double cutoff;
      double shift;
      double mass;
      double charge;
      void copy(mocha_lj_params ljB)
      {
	type = ljB.type;
	label = ljB.label;
	Erep = ljB.Erep;
	Eatt = ljB.Eatt;
	sigma = ljB.sigma;
	epsilon = ljB.epsilon;
	cutoff = ljB.cutoff;
	shift = ljB.shift;
	mass = ljB.mass;
        charge = ljB.charge;
      }
      
  };
  
  
  class mocha_box
  {
    
    public:
      vec side;
      vec hside;
      double min_side;
      double min_hside;
      double rmin;
      double rmax;
      double hregsize;
      double invhregsize;
      double volume;
      void copy(mocha_box in_box)
      {
	side = in_box.side;
	hside = in_box.hside;
	min_side = in_box.min_side;
	min_hside = in_box.min_hside;
        rmin = in_box.rmin;
        rmax = in_box.rmax;
        hregsize = in_box.hregsize;
        invhregsize = in_box.invhregsize;
        volume = in_box.volume;
      }
      void initialize()
      {
	min_side = 0.0;
	min_hside = 0.0;
	side = vec(0.0, 0.0, 0.0);
	hside = vec(0.0, 0.0, 0.0);
        rmin = 0.0;
        rmax = 0.0;
        hregsize = 0.0;
        invhregsize = 0.0;
        volume = 0.0;
      }
      void get_mins()
      {
	min_side = side.x();
	if(side.y() < side.x()) min_side = side.y();
	if(side.z() < side.x()) min_side = side.z();
	min_hside = hside.x();
	if(hside.y() < hside.x()) min_hside = hside.y();
	if(hside.z() < hside.x()) min_hside = hside.z();
        volume = side.x() * side.y() * side.z();
      }
      void get_hadress_geometry(double inRMIN, double inRMAX)
      {
          
          rmin = 1.0;
          rmax = 1.0;
          
          if(SPHERE_SGNL == true)
          {
              rmin = inRMIN * min_hside;
              rmax = inRMAX * min_hside;
          }
          
          if(XSLAB_SGNL == true)
          {
              rmin = inRMIN * hside.x();
              rmax = inRMAX * hside.x();
          }

          if(YSLAB_SGNL == true)
          {
              rmin = inRMIN * hside.y();
              rmax = inRMAX * hside.y();
          }

          if(ZSLAB_SGNL == true)
          {
              rmin = inRMIN * hside.z();
              rmax = inRMAX * hside.z();
          }

          hregsize = rmax - rmin;

          if(HADRESS_SGNL == true)
          {
              if(hregsize == 0.0)
              {
                  cout << "The hybrid region has size zero! Exiting...";
                  exit (1);
              }
              invhregsize = 1.0/hregsize;
          }

        }
      
  };
  
  
  vec pbcdist(vec u, vec v, mocha_box box)
  {
    
    double x, y, z;

    x = u.x() - v.x();
    y = u.y() - v.y();
    z = u.z() - v.z();

    if(x >= box.hside.x()) x = x - box.side.x();
    if(x < - box.hside.x()) x = x + box.side.x();
    if(y >= box.hside.y()) y = y - box.side.y();
    if(y < - box.hside.y()) y = y + box.side.y();
    if(z >= box.hside.z()) z = z - box.side.z();
    if(z < - box.hside.z()) z = z + box.side.z();

    return vec(x, y, z);
    
  }

  class data_list_struct
  {
  public:
      double x_0;
      double x_1;
      double x_2;
      double x_3;
      double x_4;
      void clear()
      {
          x_0 = 0;
          x_1 = 0;
          x_2 = 0;
          x_3 = 0;
          x_4 = 0;
      }
  };
  
  
  
  
  
}

#endif
