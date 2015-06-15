#ifndef _INTERACT
#define	_INTERACT

#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>


namespace mocha_interaction_nmspc
{
  
  using namespace std;
  using namespace mocha_types_nmspc;
  using namespace mocha_numbers_nmspc;
  using namespace mocha_energy_nmspc;
  
  
  
  class mocha_internal_interaction
  {
    
    public:
      int Natoms;
      int Nbonds;
      string type;
      struct pairlist *pairs;
      
      void allocate_memory();
      void copy(mocha_internal_interaction intA);
      int isbonded(int a, int b);
      
      
  };
  
  
  class mocha_ext_bonded_interaction
  {
    
    public:
      int Nmolecular_bonds;/** number of molecule pairs taking place in an intermolecular bonded interaction **/
      struct bonded_molecule_pair *pair_list; /** list of lists of bonds **/
      void allocate_memory(int n)
      {
	Nmolecular_bonds = n;
	pair_list = new struct bonded_molecule_pair [n];
      }
      
  };
  
  
  class mocha_nonbonded_interaction
  {
    
    public:
      int Ninteractions;
      mocha_lj_params *lj;
      int *label;
      void allocate_memory();
      void initialize_shift();
      int get_label(string a, string b);
      void copy(mocha_nonbonded_interaction intA);
      
      
  };
  
  
  
  
  
  class mocha_interaction
  {
    
    public:
      mocha_internal_interaction *internal;
      mocha_ext_bonded_interaction ext_bonded;
      mocha_nonbonded_interaction ext_nonbonded;
      
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  /** **************************************************** **/
  
  
  
  void mocha_internal_interaction::allocate_memory()
  {
    pairs = new struct pairlist [Nbonds];
    
  }
  
  
  void mocha_internal_interaction::copy(mocha_internal_interaction intA)
  {
    
    int i;
    
    Natoms = intA.Natoms;
    Nbonds = intA.Nbonds;
    type = intA.type;
    allocate_memory();
    
    for(i = 0; i < Nbonds; i++)
    {
      pairs[i].a = intA.pairs[i].a;
      pairs[i].b = intA.pairs[i].b;
      pairs[i].k = intA.pairs[i].k;
      pairs[i].rzero = intA.pairs[i].rzero;
      pairs[i].rmax = intA.pairs[i].rmax;
    }
    
    
  }
  
  
  
  int mocha_internal_interaction::isbonded(int a, int b)
  {
    
    int i, aa, bb, z;
    
    z = 0;
    
    for(i = 0; i < Nbonds; i++)
    {
      
      aa = pairs[i].a;
      bb = pairs[i].b;
      
      if( ( (a == aa) && (b == bb) ) || ( (a == bb) && (b == aa) ) ) z = 1;
      
    }
    
    return z;
    
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  void mocha_nonbonded_interaction::allocate_memory()
  {
      lj = new mocha_lj_params [Ninteractions];
  }
  
  
  
  
  
  void mocha_nonbonded_interaction::copy(mocha_nonbonded_interaction intA)
  {
    
    int i;
    
    Ninteractions = intA.Ninteractions;
    lj = new mocha_lj_params [intA.Ninteractions];

    for(i = 0; i < Ninteractions; i++) lj[i].copy(intA.lj[i]);
    
    
  }
  
  
  
  void mocha_nonbonded_interaction::initialize_shift()
  {
    
    int i;
    mocha_distance s;
    double s_rep, s_att;
    
    for(i = 0; i < Ninteractions; i++)
    {
      
      s.r = lj[i].cutoff / lj[i].sigma;
      s_rep = s.invpow(lj[i].Erep);
      s_att = s.invpow(lj[i].Eatt);
      lj[i].shift = 4.0 * lj[i].epsilon * ( s_rep - s_att );
      
    }
    
    
  }
  
  
  
  
  int mocha_nonbonded_interaction::get_label(string a, string b)
  {
    
    int i, lab;
    string typeA, typeB, typeC;

    lab = -1;
    
    for(i = 0; i < Ninteractions; i++)
    {
      typeA = a;
      typeA.append("_");
      typeA.append(b);
      
      typeB = b;
      typeB.append("_");
      typeB.append(a);
      
      typeC = lj[i].type;

      if( (typeA.compare(typeC)==0) || (typeB.compare(typeC)==0) )
      {
	lab = i;
	break;
      }
      
    }
    
    if(lab == -1)
    {
        cout << "Nonbonded interaction between " << a << " and " << b << " is not defined in lj_parameters.dat\nExiting..." << endl;
        exit (1);
    }
    
    return lab;
    
  }
 
  
  
  
}

#endif





