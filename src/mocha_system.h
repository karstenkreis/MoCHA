#ifndef _SYS
#define _SYS

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include <random>
//#include <omp.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/times.h>
#include <boost/algorithm/string.hpp>
#include <chrono>

#include "mocha_molecules.h"
#include "mocha_interaction.h"
#include "mocha_matrix_inversion.h"
//#include "mocha_types.h"
//#include "my_random.h"

#define MAXINT (70)
#define PARALLEL___
#define PARALLEL_SWEEP___
#define PARALLEL_NB___
#define PARALLEL_DENS___
#define DORHO

namespace mocha_system_nmspc
{

  using namespace std;
  using namespace votca::tools;
  using namespace mocha_molecules_nmspc;
  using namespace mocha_interaction_nmspc;
  using namespace mocha_types_nmspc;
  using namespace boost;
  using namespace mocha_matrixInv_nmspc;
 // using namespace mocha_myrandom_nmspc;



  class mocha_system
  {

    public:

      string input_folder;
      int Natoms;
      int Nmolecules;
      int Nmol_types;
      int Nthreads = 1;
      struct molecules_per_type_struct
      {
          int *number;
          int *atoms;
          string *type;
      } molecules_per_type;
      struct free_energy_comp_struct
      {
          double *f;
          double *p;
      } *free_energy_comp;
      struct density_vector_struct
      {
          int *Xdens;
          int *Ydens;
          int *Zdens;
          int *Rdens;
      };
      struct density_struct
      {
          int Xbins;
          int Ybins;
          int Zbins;
          int Rbins;
          double invXdelta;
          double invYdelta;
          double invZdelta;
          double invRdelta;
          density_vector_struct *dens_vector;
      } *density, *FEC_density;
      int fec_bins;

      struct gyration_struct
      {
         int bins;
         double invDelta;
         int* counts;
         double* gyr;
      } *gyration;

      int partQM;
      double tot_energy;
      double rad_of_gyr;
      double lambda_AA;
      double lambda_CG;
      struct lj_table_struct
      {
          int *row;
      } *lj_table;
      struct fec_poly_params_struct
      {
          double *a;
      } *fec_poly_params, *b, *xmatrix, *fec_poly_params_rho, *b_rho, *xmatrix_rho;
      int Num_fec_poly_params = 6;
      double epsilon_param_rho = 0.5;
      bool UpdateFlag = false;




      double NB_CUTOFF;
      double NB_SKIN_HALF;
      double KTI_BEGIN;
      double KTI_END;
      double LAMBDA_VALUE;
      double INV_NB_CUTOFF;
      int INNER_SWEEPS;
      int OUTER_SWEEPS;
      int EQUILIBRATE;
      double BETA;
      int KIRKW_STEPS;
      double RMIN;
      double RMAX;
      double SIGMA_DISPL;
      double SIGMA_ROTAT;
      double SIGMA_BREATH;
      double SIGMA_COLLAPSE;
      bool CGENERGY;
      bool AAENERGY;
      bool USEFEC;
      bool USEPOLYFEC;
      int XDENS;
      int YDENS;
      int ZDENS;
      int RDENS;
      int GYR_BINS;
      int PRINT_TRAJ;

      int *ATTEMPTS;
      int *SUCCESSES;
      int VerletListUpdates;

      double EnergyTime;

      mocha_box box;
      mocha_molecule *molecule;
      mocha_interaction interaction;

      /************************************************************************/
      /** I/O ROUTINES **/
      void get_threads();
      void read_params();
      void dump_info();
      void get_data_folder(int in_argc, char *in_argv[]);
      void append_GRO_conf_unfolded(const char * filename, int t);
      void append_XYZ_conf_unfolded(const char * filename, int t);
      void dump_GRO_conf_unfolded(const char * filename, int t);
      void dump_GRO_conf(const char * filename, int t);
      void dump_GRO_CG_conf(const char * filename, int t);
      void dump_MOC_conf(const char * filename, int t);

      /************************************************************************/
      /** SYSTEM ROUTINES **/
      void initialize_molecules();
      void fill_molecules();
      void fill_GRO_molecules();
      void fill_MOC_molecules();
      void read_GRO_input_file();
      void read_MOC_input_file();
      void read_input_file();
      void initialize_interaction();
      void initialize_attempt_success();
      void read_internal_interaction_file();
      void read_bonded_intramol_interaction_file();
      void read_lj_parameters();
      void total_energy();
      void update_neighbors();
      void update_neighbors_PI();
      void com_shift_center();
      void com_shift_corner();
      void get_molecule_masses();
      void update_molecule_com();
      void sweep();
      void parallel_sweep();
      void build_new_atoms();
      void nonbonded_PI_molecular_energy(int ii, int new_coord, double *energy, double *virial, double *unw_energy, double *unw_virial);
      void nonbonded_AA_molecular_energy(int ii, int new_coord, double *energy, double *virial, double *unw_energy, double *unw_virial);
      void nonbonded_CG_molecular_energy(int ii, int new_coord, double *energy, double *virial, double *unw_energy, double *unw_virial);
      void bonded_AA_molecular_energy(int mol_index, int new_coord, double *energy, double *virial);
      void internal_energy(int mol_index, int new_coord, double *energy, double *virial, double *e_diff);
      //void delta_energy(int i, double *delta_erg, double *int_energy, double *ext_energy, double *kirkw_diff, double *virial);
      data_list_struct delta_energy_function(int i);
      void kirkwood_energy_difference(double *kirkw, double *virial);
      void number_of_molecules_per_type(void);
      void read_and_allocate_fec(void);
      void initialize_density_structure(void);
      void initialize_gyration_structure(void);
      void calc_density(int direction);
      void calc_gyration();
      void calc_partQM();
      void mean_gyration();
      void print_density(int ticket);
      void print_gyration(int ticket);
      void test_func();
      void initialize_fec_poly_params(void);
      double poly_fec(int type, double lambda);
      void poly_fec_coefficient_update(void);
      void poly_fec_coefficient_calc(void);
      //void init_mersenne_twister();


  };


  /** **************************************************** **/





  void mocha_system::poly_fec_coefficient_calc(void)
  {

      int i, j, n;
      double **m, **invm;
      int k, l, nbins, type;
      double lambda, deltarho, x;

      m = new double*[Num_fec_poly_params];
      invm = new double*[Num_fec_poly_params];
      for(i = 0; i < Num_fec_poly_params; i++)
      {
          m[i] = new double[Num_fec_poly_params];
          invm[i] = new double[Num_fec_poly_params];
      }

      for(n = 0; n < Nmol_types; n++)
      {

          for(i = 0; i < Num_fec_poly_params; i++)
          {
              for(j = 0; j < Num_fec_poly_params; j++) m[i][j] = inverse_frame*xmatrix[n].a[i + Num_fec_poly_params*j];
          }

          my_invert_symmetric_matrix(m, invm, Num_fec_poly_params);

          for(i = 0; i < Num_fec_poly_params; i++)
          {
              fec_poly_params[n].a[i] = 0.0;
              for(j = 0; j < Num_fec_poly_params; j++) fec_poly_params[n].a[i] += invm[i][j]*inverse_frame*b[n].a[j];
          }

      }


#ifdef DORHO
      if(inverse_frame < 1.0)
      {
          nbins = FEC_density[0].Xbins;
          for(type = 0; type < Nmol_types; type++)
          {
              for(i = 0; i < Num_fec_poly_params; i++) b_rho[type].a[i] = 0.0;nbins = FEC_density[0].Xbins;
              for(k = 0; k < nbins; k++)
              {

                  x = k / FEC_density[0].invXdelta - box.hside.getX();
                  x = fabs(x);

                  if((x > box.rmin) && (x < box.rmax))
                  {
                      l = k;
                  }
                  else
                      continue;

                  deltarho = FEC_density[0].dens_vector[type].Xdens[l] / ( (3*OUTPUTFREQ) * (double)molecules_per_type.number[type] / (double)FEC_density[0].Xbins) - 1.0;
                  lambda = molecule[0].lambda_function(x+0.5/FEC_density[0].invXdelta);

                  for(i = 0; i < Num_fec_poly_params; i++) b_rho[type].a[i] += -epsilon_param_rho * deltarho * pow(lambda, i + 1) / ((double)i + 1.0);

              }
          }
          for(type = 0; type < Nmol_types; type++) for(k = 0; k < nbins; k++) FEC_density[0].dens_vector[type].Xdens[k] = 0;
      }

      for(n = 0; n < Nmol_types; n++)
      {
          for(i = 0; i < Num_fec_poly_params; i++)
          {
              for(j = 0; j < Num_fec_poly_params; j++) fec_poly_params_rho[n].a[i] += xmatrix_rho[n].a[i + Num_fec_poly_params*j]*b_rho[n].a[j];
          }
      }
#endif



      for(i = 0; i < Num_fec_poly_params; i++)
      {
          free(m[i]);
          free(invm[i]);
      }

      free(m);
      free(invm);

  }















  void mocha_system::poly_fec_coefficient_update(void)
  {

      int i, j, n, type;
      double lambda, h;

      for(n = 0; n < Nmolecules; n++)
      {

          if( molecule[n].is_in_HR == false) continue;
          lambda = molecule[n].lambda;
          type = molecule[n].type_number;
          h = molecule[n].unweighted_energy_diff;

          for(i = 0; i < Num_fec_poly_params; i++)
          {
              b[type].a[i] += h*pow(lambda, i);
              for(j = 0; j < Num_fec_poly_params; j++) xmatrix[type].a[i + Num_fec_poly_params*j] += pow(lambda, i+j);
          }

      }

  }


  void mocha_system::initialize_fec_poly_params(void)
  {

      int i, j, k, l, n, nbins, type;
      double lambda, x;
      double **m, **invm;

      m = new double*[Num_fec_poly_params];
      invm = new double*[Num_fec_poly_params];
      for(i = 0; i < Num_fec_poly_params; i++)
      {
          m[i] = new double[Num_fec_poly_params];
          invm[i] = new double[Num_fec_poly_params];
      }




      fec_poly_params = new fec_poly_params_struct[Nmol_types];
      b = new fec_poly_params_struct[Nmol_types];
      xmatrix = new fec_poly_params_struct[Nmol_types];
      for(i = 0; i < Nmol_types; i++)
      {
          fec_poly_params[i].a = new double[Num_fec_poly_params];
          b[i].a = new double[Num_fec_poly_params];
          xmatrix[i].a = new double[Num_fec_poly_params*Num_fec_poly_params];
          for(j = 0; j < Num_fec_poly_params; j++)
          {
              fec_poly_params[i].a[j] = 0.0;
              b[i].a[j] = 0.0;
              for(k = 0; k < Num_fec_poly_params; k++) xmatrix[i].a[j + Num_fec_poly_params*k] = 0.0;
          }
      }

      fec_poly_params_rho = new fec_poly_params_struct[Nmol_types];
      b_rho = new fec_poly_params_struct[Nmol_types];
      xmatrix_rho = new fec_poly_params_struct[Nmol_types];
      for(i = 0; i < Nmol_types; i++)
      {
          fec_poly_params_rho[i].a = new double[Num_fec_poly_params];
          b_rho[i].a = new double[Num_fec_poly_params];
          xmatrix_rho[i].a = new double[Num_fec_poly_params*Num_fec_poly_params];
          for(j = 0; j < Num_fec_poly_params; j++)
          {
              fec_poly_params_rho[i].a[j] = 0.0;
              b_rho[i].a[j] = 0.0;
              for(k = 0; k < Num_fec_poly_params; k++) xmatrix_rho[i].a[j + Num_fec_poly_params*k] = 0.0;
          }
      }

      nbins = FEC_density[0].Xbins;
      for(type = 0; type < Nmol_types; type++)
          {

              for(k = 0; k < nbins; k++)
              {

                  x = k / FEC_density[0].invXdelta - box.hside.getX();
                  x = fabs(x);

                  if((x > box.rmin) && (x < box.rmax))
                  {
                      l = k;
                  }
                  else
                      continue;

                  lambda = molecule[0].lambda_function(x+0.5/FEC_density[0].invXdelta);

                  for(i = 0; i < Num_fec_poly_params; i++)
                  {
                      for (j = 0; j < Num_fec_poly_params; j++) xmatrix_rho[type].a[i + Num_fec_poly_params * j] += pow(lambda, i + j + 2) / (((double)i + 1.0)*((double)j + 1.0));
                  }
              }
          }


      for(i = 0; i < Num_fec_poly_params; i++)
      {
          for(j = 0; j < Num_fec_poly_params; j++) m[i][j] = 0.0;
          for(j = 0; j < Num_fec_poly_params; j++) invm[i][j] = 0.0;
      }

      for(n = 0; n < Nmol_types; n++)
      {

          for(i = 0; i < Num_fec_poly_params; i++)
          {
              for(j = 0; j < Num_fec_poly_params; j++) m[i][j] = xmatrix_rho[n].a[i + Num_fec_poly_params*j];
          }

          my_invert_symmetric_matrix(m, invm, Num_fec_poly_params);

          for(i = 0; i < Num_fec_poly_params; i++)
          {
              for(j = 0; j < Num_fec_poly_params; j++) xmatrix_rho[n].a[i + Num_fec_poly_params*j] = invm[i][j];
          }

      }

      for(i = 0; i < Num_fec_poly_params; i++)
      {
          free(m[i]);
          free(invm[i]);
      }

      free(m);
      free(invm);



  }


  double mocha_system::poly_fec(int type, double lambda)
  {
      int i;
      double x, y;
      x = 0.0;
      y = 1.0;

      for(i = 0; i < Num_fec_poly_params; i++)
      {
          //y = 1.0;
          //for(j = 1; j <= i+1; j++) y = y*lambda;
          y = pow(lambda, i+1)/((double)i + 1.0);
          x += y*(fec_poly_params[type].a[i] + fec_poly_params_rho[type].a[i]);
      }
      //cout << "aaa " << lambda << " " << x << endl;
      return x;

  }




  void mocha_system::read_and_allocate_fec(void)
  {

      stringstream filename;
      filename << input_folder << "/fec.dat";

      int i, j, k, l, n, t;
      double f, p;
      string type;

      ifstream inFile(filename.str().c_str(), ios::in);

      inFile >> n >> t;
      fec_bins = n;

      free_energy_comp = new free_energy_comp_struct[n];
      for(i = 0; i < n; i++)
      {
          free_energy_comp[i].f = new double[t];
          free_energy_comp[i].p = new double[t];
      }

      for(k = 0; k < n*t; k++)
      {

          inFile >> i >> f >> p >> type;
          trim(type);

          l = -1;
          for(j = 0; j < Nmol_types; j++)
          {
              if(type.compare(molecules_per_type.type[j])==0)
              {
                  l = j;
                  break;
              }
          }

          if(l<0)
          {
              cout << "Warning! In free energy compensation file the molecule type label is missing! Exiting......\n";
              exit(1);
          }

          if(ATOMLAM == true)
          {
              free_energy_comp[i].f[l] = (double)f*Nmolecules/Natoms;
              free_energy_comp[i].p[l] = (double)p*Nmolecules/Natoms;
          }
          else
          {
              free_energy_comp[i].f[l] = f;
              free_energy_comp[i].p[l] = p;
          }

          //cout << "free_energy_comp[i].f[l] = " << free_energy_comp[i].f[l] << endl;
          //cout << "free_energy_comp[i].p[l] = " << free_energy_comp[i].p[l] << endl;



      }


  }


  void mocha_system::initialize_molecules()
  {

    int i;

    molecule = new mocha_molecule[Nmolecules];

    for(i = 0; i < Nmolecules; i++) molecule[i].Natoms = 0;

  }


  void mocha_system::get_molecule_masses()
  {

    int i;

    for(i = 0; i < Nmolecules; i++)
    {
      molecule[i].calc_mass();
    }

  }



  void mocha_system::update_molecule_com()
  {

    int i;
    vec temp;

    for(i = 0; i < Nmolecules; i++)
    {
        temp = molecule[i].update_com();
        molecule[i].get_lambda(0);
    }

  }




  void mocha_system::initialize_interaction()
  {
    interaction.internal = new mocha_internal_interaction [Nmol_types];
  }



  void mocha_system::initialize_attempt_success()
  {
    int i;
    ATTEMPTS = new int [4];
    SUCCESSES = new int [4];
    for (i = 0; i < 4; i++)
    {
        ATTEMPTS[i] = 0;
        SUCCESSES[i] = 0;
    }
  }



  void mocha_system::read_input_file()
  {

    int i, foo, local_Natoms, local_Nmolecules, local_Nmol_types, mol_index;
    double x, y, z;
    string local_label, mol_type;

    local_Natoms = 0;
    local_Nmolecules = 0;

    stringstream filename1, filename2;
    filename1 << input_folder << "/input.xyz";

    ifstream inFile(filename1.str().c_str(), ios::in);
    inFile >> local_Natoms;
    inFile >> x >> y >> z;

    box.side.x() = x;
    box.side.y() = y;
    box.side.z() = z;

    box.hside = box.side / 2.0;
    box.get_mins();

    box.get_hadress_geometry(RMIN, RMAX);

    foo = 0;

    for(i = 0; i < local_Natoms; i++)
    {

      inFile >> local_label >> x >> y >> z >> mol_index >> mol_type;

      if(i == 0)
      {
	foo = mol_index;
	local_Nmolecules++;
      }
      else
      {
	if(mol_index != foo)
	{
	  local_Nmolecules++;
	  foo = mol_index;
	}
      }

    }

    inFile.close();

    /** read the number of different molecule types **/
    filename2 << input_folder << "/internal_topology.dat";
    ifstream interactFile(filename2.str().c_str(), ios::in);
    interactFile >> local_Nmol_types;
    interactFile.close();

    Natoms = local_Natoms;
    Nmolecules = local_Nmolecules;
    Nmol_types = local_Nmol_types;


  }











  void mocha_system::read_GRO_input_file()
  {

    int i, foo, local_Natoms, local_Nmolecules, local_Nmol_types, mol_index;
    double x, y, z, v, temp;
    string local_label, mol_type, foo_string;
    char buffer[100];

    istringstream value;

    local_Natoms = 0;
    local_Nmolecules = 0;

    stringstream filename1, filename2, stream;
    filename1 << input_folder << "/input.gro";

    ifstream inFile(filename1.str().c_str(), ios::in);

    getline(inFile, foo_string);
    //inFile >> local_Natoms;

    getline(inFile, foo_string);
    stream << foo_string;
    stream.str().copy(buffer,100,0);
    value.str(buffer);
    value >> local_Natoms;
    value.clear();
    stream.str("");
    memset(buffer, 0, 100);

    foo = 0;

    for(i = 0; i < local_Natoms; i++)
    {
/*
      inFile >> setw(5) >> mol_index\
                >> setw(5) >> mol_type\
                >> setw(5) >> local_label\
                >> setw(5) >> v\
                >> setw(8) >> setprecision(3) >> fixed >> x\
                >> setw(8) >> setprecision(3) >> fixed >> y\
                >> setw(8) >> setprecision(3) >> fixed >> z\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v;
 */


        getline(inFile, foo_string);

      stream << foo_string;
      stream.str().copy(buffer,5,0);
      value.str(buffer);
      value >> mol_index;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,5,5);
      mol_type = buffer;
      trim(mol_type);
      memset(buffer, 0, 100);
      stream.str("");

      stream << foo_string;
      stream.str().copy(buffer,5,10);
      local_label = buffer;
      trim(local_label);
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,5,15);
      value.str(buffer);
      value >> temp;
      v = (int)temp;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,8,20);
      value.str(buffer);
      value >> x;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,8,28);
      value.str(buffer);
      value >> y;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,8,36);
      value.str(buffer);
      value >> z;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);



      if(i == 0)
      {
	foo = mol_index;
	local_Nmolecules++;
      }
      else
      {
	if(mol_index != foo)
	{
	  local_Nmolecules++;
	  foo = mol_index;
	}
      }

    }
    //cout << "local_Nmolecules   tttttt " << local_Nmolecules <<  endl;

    inFile >> x >> y >> z;

    box.side.x() = x;
    box.side.y() = y;
    box.side.z() = z;

    box.hside = box.side / 2.0;
    box.get_mins();

    box.get_hadress_geometry(RMIN, RMAX);

    inFile.close();

    /** read the number of different molecule types **/
    filename2 << input_folder << "/internal_topology.dat";
    ifstream interactFile(filename2.str().c_str(), ios::in);
    interactFile >> local_Nmol_types;
    interactFile.close();

    Natoms = local_Natoms;
    cout << "Natoms = " << Natoms << endl;
    Nmolecules = local_Nmolecules;
    cout << "Nmolecules = " << Nmolecules << endl;
    Nmol_types = local_Nmol_types;
    //cout << "Natoms " << Natoms << " Nmolecules " << Nmolecules << " Nmol_types " << Nmol_types <<  endl;


  }






  void mocha_system::read_MOC_input_file()
  {

    int i, foo, local_Natoms, local_Nmolecules, local_Nmol_types, mol_index;
    double x, y, z, v, temp;
    string local_label, mol_type, foo_string;
    char buffer[100];

    istringstream value;

    local_Natoms = 0;
    local_Nmolecules = 0;

    stringstream filename1, filename2, stream;
    filename1 << input_folder << "/input.moc";

    ifstream inFile(filename1.str().c_str(), ios::in);

    getline(inFile, foo_string);
    //inFile >> local_Natoms;

    getline(inFile, foo_string);
    stream << foo_string;
    stream.str().copy(buffer,100,0);
    value.str(buffer);
    value >> local_Natoms;
    value.clear();
    stream.str("");
    memset(buffer, 0, 100);

    foo = 0;

    for(i = 0; i < local_Natoms; i++)
    {
/*
      inFile >> setw(5) >> mol_index\
                >> setw(5) >> mol_type\
                >> setw(5) >> local_label\
                >> setw(5) >> v\
                >> setw(8) >> setprecision(3) >> fixed >> x\
                >> setw(8) >> setprecision(3) >> fixed >> y\
                >> setw(8) >> setprecision(3) >> fixed >> z\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v;
 */


        getline(inFile, foo_string);

      stream << foo_string;
      stream.str().copy(buffer,5,0);
      value.str(buffer);
      value >> mol_index;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,5,5);
      mol_type = buffer;
      trim(mol_type);
      memset(buffer, 0, 100);
      stream.str("");

      stream << foo_string;
      stream.str().copy(buffer,5,10);
      local_label = buffer;
      trim(local_label);
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,5,15);
      value.str(buffer);
      value >> temp;
      v = (int)temp;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,10,20);
      value.str(buffer);
      value >> x;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,10,30);
      value.str(buffer);
      value >> y;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo_string;
      stream.str().copy(buffer,10,40);
      value.str(buffer);
      value >> z;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);



      if(i == 0)
      {
	foo = mol_index;
	local_Nmolecules++;
      }
      else
      {
	if(mol_index != foo)
	{
	  local_Nmolecules++;
	  foo = mol_index;
	}
      }

    }
    //cout << "local_Nmolecules   tttttt " << local_Nmolecules <<  endl;

    inFile >> x >> y >> z;

    box.side.x() = x;
    box.side.y() = y;
    box.side.z() = z;

    box.hside = box.side / 2.0;
    box.get_mins();

    box.get_hadress_geometry(RMIN, RMAX);

    inFile.close();

    /** read the number of different molecule types **/
    filename2 << input_folder << "/internal_topology.dat";
    ifstream interactFile(filename2.str().c_str(), ios::in);
    interactFile >> local_Nmol_types;
    interactFile.close();

    Natoms = local_Natoms;
    Nmolecules = local_Nmolecules;
    Nmol_types = local_Nmol_types;
    //cout << "Natoms " << Natoms << " Nmolecules " << Nmolecules << " Nmol_types " << Nmol_types <<  endl;


  }






  void mocha_system::fill_molecules()
  {

    int i, mol_index, atom_index;
    double x, y, z;
    string foo, mol_type;

    stringstream filename;
    filename << input_folder << "/input.xyz";
    ifstream inFile(filename.str().c_str(), ios::in);
    inFile >> foo;
    inFile >> x >> y >> z;


    for(i = 0; i < Natoms; i++)
    {
      inFile >> foo >> x >> y >> z >> mol_index >> mol_type;
      molecule[mol_index].Natoms++;
    }

    int *n = new int[Nmolecules];
    for(i = 0; i < Nmolecules; i++)
    {
      molecule[i].allocate_memory(Nmolecules);
      molecule[i].sigma_displ = SIGMA_DISPL;
      molecule[i].sigma_rotat = SIGMA_ROTAT;
      molecule[i].sigma_breath = SIGMA_BREATH;
      molecule[i].sigma_collapse = SIGMA_COLLAPSE;
      n[i] = 0;
    }

    inFile.seekg (0, ios::beg);
    inFile >> foo;
    inFile >> x >> y >> z;

    for(i = 0; i < Natoms; i++)
    {

      inFile >> foo >> x >> y >> z >> mol_index >> mol_type;

      atom_index = n[mol_index];

      molecule[mol_index].atom[atom_index].position.x() = x;
      molecule[mol_index].atom[atom_index].position.y() = y;
      molecule[mol_index].atom[atom_index].position.z() = z;
      molecule[mol_index].atom[atom_index].type = foo;


      molecule[mol_index].index = mol_index;
      molecule[mol_index].type = mol_type;

      n[mol_index] += 1;

    }


    inFile.close();

    /* copy in each molecule's box structure the system's box data */
    for(i = 0; i < Nmolecules; i++)
    {
        molecule[i].box.copy(box);

    }


  }


  void mocha_system::fill_GRO_molecules()
  {

    int i, j, mol_index, atom_index, check, shift;
    double x, y, z, temp;
    string foo, mol_type, local_label;
    char buffer[100];

    stringstream filename, stream;
    istringstream value;
    filename << input_folder << "/input.gro";
    ifstream inFile(filename.str().c_str(), ios::in);

    getline(inFile, foo);
    getline(inFile, foo);

    check = 0;
    shift = 0;
    for(i = 0; i < Natoms; i++)
    {
      /*
       inFile >> setw(5) >> mol_index\
                >> setw(5) >> mol_type\
                >> setw(5) >> local_label\
                >> setw(5) >> j\
                >> setw(8) >> setprecision(3) >> fixed >> x\
                >> setw(8) >> setprecision(3) >> fixed >> y\
                >> setw(8) >> setprecision(3) >> fixed >> z\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v;
*/
        getline(inFile, foo);

            stream << foo;
            stream.str().copy(buffer, 5, 0);
            value.str(buffer);
            value >> mol_index;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 5, 5);
            mol_type = buffer;
	    trim(mol_type);
            memset(buffer, 0, 100);
            stream.str("");

            stream << foo;
            stream.str().copy(buffer, 5, 10);
            local_label = buffer;
	    trim(local_label);
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 5, 15);
            value.str(buffer);
            value >> temp;
            j = (int) temp;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 8, 20);
            value.str(buffer);
            value >> x;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 8, 28);
            value.str(buffer);
            value >> y;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 8, 36);
            value.str(buffer);
            value >> z;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            if(check==0)
            {
                if(mol_index == 1) shift = 1;

                check = 1;

                if((mol_index!=0)&&(mol_index!=1))
                {
                    cout << "The molecules in the .gro file are not properly indexed! Exiting...\n";
                    exit(1);
                }
            }

            if(shift==1) mol_index = mol_index-1;

            molecule[mol_index].Natoms++;
    }

    int *n = new int[Nmolecules];
    for(i = 0; i < Nmolecules; i++)
    {
      molecule[i].allocate_memory(Nmolecules);
      molecule[i].sigma_displ = SIGMA_DISPL;
      molecule[i].sigma_rotat = SIGMA_ROTAT;
      molecule[i].sigma_breath = SIGMA_BREATH;
      molecule[i].sigma_collapse = SIGMA_COLLAPSE;
      n[i] = 0;
    }

    inFile.seekg (0, ios::beg);

    getline(inFile, foo);
    getline(inFile, foo);


    check = 0;
    shift = 0;
    for(i = 0; i < Natoms; i++)
    {

      /*inFile >> setw(5) >> mol_index\
                >> setw(5) >> mol_type\
                >> setw(5) >> local_label\
                >> setw(5) >> j\
                >> setw(8) >> setprecision(3) >> fixed >> x\
                >> setw(8) >> setprecision(3) >> fixed >> y\
                >> setw(8) >> setprecision(3) >> fixed >> z\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v;
       */

      getline(inFile, foo);

      stream << foo;
      stream.str().copy(buffer,5,0);
      value.str(buffer);
      value >> mol_index;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,5,5);
      mol_type = buffer;
      trim(mol_type);
      memset(buffer, 0, 100);
      stream.str("");

      stream << foo;
      stream.str().copy(buffer,5,10);
      local_label = buffer;
      trim(local_label);
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,5,15);
      value.str(buffer);
      value >> temp;
      j = (int)temp;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,8,20);
      value.str(buffer);
      value >> x;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,8,28);
      value.str(buffer);
      value >> y;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,8,36);
      value.str(buffer);
      value >> z;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);



      if(check==0)
            {
                if(mol_index == 1) shift = 1;

                check = 1;

                if((mol_index!=0)&&(mol_index!=1))
                {
                    cout << "The molecules in the .gro file are not properly indexed! Exiting...\n";
                    exit(1);
                }
            }

            if(shift==1) mol_index = mol_index-1;

      atom_index = n[mol_index];

      molecule[mol_index].atom[atom_index].position.x() = x;
      molecule[mol_index].atom[atom_index].position.y() = y;
      molecule[mol_index].atom[atom_index].position.z() = z;
      molecule[mol_index].atom[atom_index].type = local_label;
      //cout << "atom type " << molecule[mol_index].atom[atom_index].type << endl;


      molecule[mol_index].index = mol_index;
      molecule[mol_index].type = mol_type;
      //cout << "index " << mol_index << " type " << molecule[mol_index].type << endl;

      n[mol_index] += 1;

    }

    inFile >> x >> y >> z;
    //cout << "BOX " << x << " " << y << " " << z << endl;

    inFile.close();

    /* copy in each molecule's box structure the system's box data */
    for(i = 0; i < Nmolecules; i++)
    {
        molecule[i].box.copy(box);

    }


  }






  void mocha_system::fill_MOC_molecules()
  {

    int i, j, mol_index, atom_index, check, shift;
    double x, y, z, temp;
    string foo, mol_type, local_label;
    char buffer[100];

    stringstream filename, stream;
    istringstream value;
    filename << input_folder << "/input.moc";
    ifstream inFile(filename.str().c_str(), ios::in);

    getline(inFile, foo);
    getline(inFile, foo);

    check = 0;
    shift = 0;
    for(i = 0; i < Natoms; i++)
    {
      /*
       inFile >> setw(5) >> mol_index\
                >> setw(5) >> mol_type\
                >> setw(5) >> local_label\
                >> setw(5) >> j\
                >> setw(8) >> setprecision(3) >> fixed >> x\
                >> setw(8) >> setprecision(3) >> fixed >> y\
                >> setw(8) >> setprecision(3) >> fixed >> z\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v;
*/
        getline(inFile, foo);

            stream << foo;
            stream.str().copy(buffer, 5, 0);
            value.str(buffer);
            value >> mol_index;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 5, 5);
            mol_type = buffer;
	    trim(mol_type);
            memset(buffer, 0, 100);
            stream.str("");

            stream << foo;
            stream.str().copy(buffer, 5, 10);
            local_label = buffer;
	    trim(local_label);
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 5, 15);
            value.str(buffer);
            value >> temp;
            j = (int) temp;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 10, 20);
            value.str(buffer);
            value >> x;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 10, 30);
            value.str(buffer);
            value >> y;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            stream << foo;
            stream.str().copy(buffer, 10, 40);
            value.str(buffer);
            value >> z;
            value.clear();
            stream.str("");
            memset(buffer, 0, 100);

            if(check==0)
            {
                if(mol_index == 1) shift = 1;

                check = 1;

                if((mol_index!=0)&&(mol_index!=1))
                {
                    cout << "The molecules in the .gro file are not properly indexed! Exiting...\n";
                    exit(1);
                }
            }

            if(shift==1) mol_index = mol_index-1;

            molecule[mol_index].Natoms++;
    }

    int *n = new int[Nmolecules];
    for(i = 0; i < Nmolecules; i++)
    {
      molecule[i].allocate_memory(Nmolecules);
      molecule[i].sigma_displ = SIGMA_DISPL;
      molecule[i].sigma_rotat = SIGMA_ROTAT;
      molecule[i].sigma_breath = SIGMA_BREATH;
      molecule[i].sigma_collapse = SIGMA_COLLAPSE;
      n[i] = 0;
    }

    inFile.seekg (0, ios::beg);

    getline(inFile, foo);
    getline(inFile, foo);


    check = 0;
    shift = 0;
    for(i = 0; i < Natoms; i++)
    {

      /*inFile >> setw(5) >> mol_index\
                >> setw(5) >> mol_type\
                >> setw(5) >> local_label\
                >> setw(5) >> j\
                >> setw(8) >> setprecision(3) >> fixed >> x\
                >> setw(8) >> setprecision(3) >> fixed >> y\
                >> setw(8) >> setprecision(3) >> fixed >> z\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v\
                >> setw(8) >> setprecision(3) >> fixed >> v;
       */

      getline(inFile, foo);

      stream << foo;
      stream.str().copy(buffer,5,0);
      value.str(buffer);
      value >> mol_index;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,5,5);
      mol_type = buffer;
      trim(mol_type);
      memset(buffer, 0, 100);
      stream.str("");

      stream << foo;
      stream.str().copy(buffer,5,10);
      local_label = buffer;
      trim(local_label);
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,5,15);
      value.str(buffer);
      value >> temp;
      j = (int)temp;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,10,20);
      value.str(buffer);
      value >> x;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,10,30);
      value.str(buffer);
      value >> y;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);

      stream << foo;
      stream.str().copy(buffer,10,40);
      value.str(buffer);
      value >> z;
      value.clear();
      stream.str("");
      memset(buffer, 0, 100);



      if(check==0)
            {
                if(mol_index == 1) shift = 1;

                check = 1;

                if((mol_index!=0)&&(mol_index!=1))
                {
                    cout << "The molecules in the .gro file are not properly indexed! Exiting...\n";
                    exit(1);
                }
            }

            if(shift==1) mol_index = mol_index-1;

      atom_index = n[mol_index];

      molecule[mol_index].atom[atom_index].position.x() = x;
      molecule[mol_index].atom[atom_index].position.y() = y;
      molecule[mol_index].atom[atom_index].position.z() = z;
      molecule[mol_index].atom[atom_index].type = local_label;
      //cout << "atom type " << molecule[mol_index].atom[atom_index].type << endl;


      molecule[mol_index].index = mol_index;
      molecule[mol_index].type = mol_type;
      //cout << "index " << mol_index << " type " << molecule[mol_index].type << endl;

      n[mol_index] += 1;

    }

    inFile >> x >> y >> z;
    //cout << "BOX " << x << " " << y << " " << z << endl;

    inFile.close();

    /* copy in each molecule's box structure the system's box data */
    for(i = 0; i < Nmolecules; i++)
    {
        molecule[i].box.copy(box);

    }


  }






  void mocha_system::build_new_atoms()
  {

      int i, j;

      for (i = 0; i < Nmolecules; i++)
      {
          for (j = 0; j < molecule[i].Natoms; j++) molecule[i].new_atom[j].copy(molecule[i].atom[j]);
      }

  }









  void mocha_system::read_internal_interaction_file()
  {

    stringstream filename;
    filename << input_folder << "/internal_topology.dat";
    ifstream interactFile(filename.str().c_str(), ios::in);


    int i, j, a, b, Nmoltype, local_Natoms, local_Nbonds;
    double k, rzero, rmax;
    string type, typeA;

    interactFile >> Nmoltype;

    for(i = 0; i < Nmoltype; i++)
    {

      interactFile >> local_Natoms >> local_Nbonds >> interaction.internal[i].type;

      interaction.internal[i].Natoms = local_Natoms;
      interaction.internal[i].Nbonds = local_Nbonds;

      interaction.internal[i].allocate_memory();

      for(j = 0; j < local_Nbonds; j++)
      {
	interactFile >> a >> b >> k >> rzero >> rmax;
	interaction.internal[i].pairs[j].a = a;
	interaction.internal[i].pairs[j].b = b;
	interaction.internal[i].pairs[j].k = BETA*k;
	interaction.internal[i].pairs[j].rzero = rzero;
	interaction.internal[i].pairs[j].rmax = rmax;
      }


    }

    interactFile.close();


    for(i = 0; i < Nmolecules; i++)
    {

      for(j = 0; j < Nmoltype; j++)
      {

	type = molecule[i].type;
        trim(type);
	typeA = interaction.internal[j].type;
        trim(typeA);

	if(type.compare(typeA)==0)
	{
	  molecule[i].internal_int.copy(interaction.internal[j]);
          molecule[i].type_number = j;
	}

      }


    }

  }












  void mocha_system::read_bonded_intramol_interaction_file()
  {

    int i, k, a, b, l, m, n, Nbonds;
    double elast_k, rzero, rmax;

    stringstream filename;
    filename << input_folder << "/bonded_intramolecular_topology.dat";
    ifstream interactFile(filename.str().c_str(), ios::in);

    interactFile >> Nbonds;

    if(Nbonds == 0) cout << "There are no molecule pairs with\n    intra-molecular bonded interactions\n";
    if(Nbonds == 1) cout << "There is a pair of molecules with\n    intra-molecular bonded interactions\n";
    if(Nbonds > 1) cout << "There are " << Nbonds << " pairs of molecules with\n    intra-molecular bonded interactions\n";

    interaction.ext_bonded.allocate_memory(Nbonds);

    for(k = 0; k < Nbonds; k++) /** cycle over the pairs of interacting molecules **/
    {

      interactFile >> n; /** number of bonds between a given pair of molecules **/

      interaction.ext_bonded.pair_list[k].Nbonds = n;

      interaction.ext_bonded.pair_list[k].pairs = new struct pairlist [n];

      for(i = 0; i < n; i++)
      {
	interactFile >> a >> b >> l >> m >> elast_k >> rzero >> rmax;

	interaction.ext_bonded.pair_list[k].a = a; /** a molecule involved in the bond **/
	interaction.ext_bonded.pair_list[k].b = b; /** a molecule involved in the bond **/
	interaction.ext_bonded.pair_list[k].pairs[i].a = l; /** the atom from molecule a involved in the bond **/
	interaction.ext_bonded.pair_list[k].pairs[i].b = m; /** the atom from molecule b involved in the bond **/
	interaction.ext_bonded.pair_list[k].pairs[i].k = BETA*elast_k; /** strength of the bond **/
	interaction.ext_bonded.pair_list[k].pairs[i].rzero = rzero; /** rest length **/
	interaction.ext_bonded.pair_list[k].pairs[i].rmax = rmax; /** max length **/
      }


    }


  }






  void mocha_system::read_lj_parameters()
  {

    int Npairs, molNatoms, i, j, k, k1, l, q;
    int Erep, Eatt;
    double sigma, epsilon, cutoff, mass, charge;
    string type, typeA, typeB;
    string type_list[MAXINT];

    for(i = 0; i < MAXINT; i++) type_list[i] = "EMPTY";

    stringstream filename;
    filename << input_folder << "/lj_parameters.dat";
    ifstream ljFile(filename.str().c_str(), ios::in);

    ljFile >> Npairs;

    interaction.ext_nonbonded.Ninteractions = Npairs;
    interaction.ext_nonbonded.allocate_memory();

    for(i = 0; i < Npairs; i++)
    {

      ljFile >> typeA >> typeB >> Erep >> Eatt >> sigma >> epsilon >> cutoff >> mass >> charge;

      trim(typeA);
      trim(typeB);
      if( (typeA.compare("EMPTY")==0) || (typeB.compare("EMPTY")==0) )
      {
          cout << "Warning! 'EMPTY' is a protected character, it cannot be an atom/molecule type!\nExiting.....\n";
          exit(1);
      }

      type = typeA;
      type.append("_");
      type.append(typeB);

      interaction.ext_nonbonded.lj[i].type = type;
      interaction.ext_nonbonded.lj[i].label = i;
      interaction.ext_nonbonded.lj[i].Erep = Erep;
      interaction.ext_nonbonded.lj[i].Eatt = Eatt;
      interaction.ext_nonbonded.lj[i].sigma = sigma;
      interaction.ext_nonbonded.lj[i].epsilon = BETA*epsilon;
      interaction.ext_nonbonded.lj[i].cutoff = cutoff;
      interaction.ext_nonbonded.lj[i].mass = mass;
      interaction.ext_nonbonded.lj[i].charge = charge;
      if(cutoff >= box.min_hside)
      {
          cout << "ERROR! Cutoff of nonbonded interaction " << type << " has cutoff larger than box half side.\nExiting..." << endl;
          exit (1);
      }




      k = 0;
      k1 = 0;
      for(q = MAXINT-1; q >= 0; q--)
      {
          if(type_list[q].compare(typeA)==0) k1 = 1;
          if(type_list[q].compare("EMPTY")==0) k = q;
      }
      if(k1==0) type_list[k] = typeA;

      k = 0;
      k1 = 0;
      for(q = MAXINT-1; q >= 0; q--)
      {
          if(type_list[q].compare(typeB)==0) k1 = 1;
          if(type_list[q].compare("EMPTY")==0) k = q;
      }
      if(k1==0) type_list[k] = typeB;


    }


    ljFile.close();



    for(i = 0; i < Nmolecules; i++)
    {

      molNatoms = molecule[i].Natoms;

      k = 0;
      for(q = 0; q < MAXINT; q++)
      {

          if(molecule[i].type.compare(type_list[q])==0)
          {
              k = q;
              break;
          }
      }
      molecule[i].interaction_number = q;

      for(j = 0; j < molNatoms; j++)
      {

                k = 0;
                for(q = 0; q < MAXINT; q++)
                {
                    if (molecule[i].atom[j].type.compare(type_list[q]) == 0)
                    {
                        k = q;
                        break;
                    }
                }
                molecule[i].atom[j].interaction_number = q;



                k = 0;

                for(l = 0; l < Npairs; l++)
                {
                    typeA = molecule[i].atom[j].type;
                    typeB = typeA;
                    typeA.append("_");
                    typeA.append(typeB);
                    typeB = interaction.ext_nonbonded.lj[l].type;

		    //cout << "atom |" << typeA << "| atom |" << typeB << "|\n";

                    if(typeA.compare(typeB) == 0)
                    {
                        k = l;
                        break;
                    }
                }

                molecule[i].atom[j].mass = interaction.ext_nonbonded.lj[k].mass;
                molecule[i].atom[j].charge = interaction.ext_nonbonded.lj[k].charge;

      }



    }

    interaction.ext_nonbonded.initialize_shift();

    lj_table = new lj_table_struct [MAXINT];
    for(i = 0; i < MAXINT; i++) lj_table[i].row = new int [MAXINT];

    for(i = 0; i < MAXINT; i++)
    {
        for(j = 0; j < MAXINT; j++)
        {
            lj_table[i].row[j] = MAXINT;
        }
    }


    for(i = 0; i < MAXINT; i++)
    {
        for(j = 0; j < MAXINT; j++)
        {
            typeA = type_list[i];
            typeB = type_list[j];
            typeA.append("_");
            typeA.append(typeB);

            q = MAXINT;
            for(k = 0; k < Npairs; k++)
            {
                if(typeA.compare(interaction.ext_nonbonded.lj[k].type)==0)
                {
                    q = k;
                    break;
                }
            }

            typeA = type_list[j];
            typeB = type_list[i];
            typeA.append("_");
            typeA.append(typeB);

            for(k = 0; k < Npairs; k++)
            {
                if(typeA.compare(interaction.ext_nonbonded.lj[k].type)==0)
                {
                    q = k;
                    break;
                }
            }

            lj_table[i].row[j] = q;
	    //cout << "Interaction " << typeA << " number " << q << " i " << i << " j " << j << endl;

        }
    }


  }





  void mocha_system::total_energy()
  {

    int i, l;
    int n, q;
    double lambda, fecpos;
    double temp_erg;
    double e_diff;
    double int_energy, ext_NB_AA_energy, ext_BB_AA_energy, ext_NB_CG_energy;
    double int_virial, ext_NB_AA_virial, ext_BB_AA_virial, ext_NB_CG_virial;
    double unw_ext_NB_AA_energy, unw_ext_NB_CG_energy;
    double unw_ext_NB_AA_virial, unw_ext_NB_CG_virial;

    temp_erg = 0.0;

    for(i = 0; i < Nmolecules; i++)
    {

      e_diff = 0.0;
      int_energy = 0.0;
      ext_NB_AA_energy = 0.0;
      ext_BB_AA_energy = 0.0;
      ext_NB_CG_energy = 0.0;

      int_virial = 0.0;
      ext_NB_AA_virial = 0.0;
      ext_BB_AA_virial = 0.0;
      ext_NB_CG_virial = 0.0;

      unw_ext_NB_AA_energy = 0.0;
      unw_ext_NB_CG_energy = 0.0;
      unw_ext_NB_AA_virial = 0.0;
      unw_ext_NB_CG_virial = 0.0;

      internal_energy(i, 0, &int_energy, &int_virial, &e_diff);
      if((AAENERGY == true) && (PATHINT == false)) nonbonded_AA_molecular_energy(i, 0, &ext_NB_AA_energy, &ext_NB_AA_virial, &unw_ext_NB_AA_energy, &unw_ext_NB_AA_virial);
      bonded_AA_molecular_energy(i, 0, &ext_BB_AA_energy, &ext_BB_AA_virial);
      if((CGENERGY == true) && (PATHINT == false)) nonbonded_CG_molecular_energy(i, 0, &ext_NB_CG_energy, &ext_NB_CG_virial, &unw_ext_NB_CG_energy, &unw_ext_NB_CG_virial);
      if(PATHINT == true) nonbonded_PI_molecular_energy(i, 0, &ext_NB_AA_energy, &ext_NB_AA_virial, &unw_ext_NB_AA_energy, &unw_ext_NB_AA_virial);

      molecule[i].int_energy = int_energy;
      if(PATHINT == false)
      {
          molecule[i].ext_energy = lambda_AA*ext_NB_AA_energy + ext_BB_AA_energy + lambda_CG*ext_NB_CG_energy;
          molecule[i].virial = lambda_AA*ext_NB_AA_virial + int_virial + ext_BB_AA_virial + lambda_CG*ext_NB_CG_virial;
          molecule[i].kirkw_diff = ext_NB_AA_energy - ext_NB_CG_energy;
      }
      else
      {
          molecule[i].ext_energy = ext_NB_AA_energy + ext_BB_AA_energy;
          molecule[i].virial = unw_ext_NB_AA_virial + int_virial + ext_BB_AA_virial;
          molecule[i].kirkw_diff = unw_ext_NB_AA_energy + e_diff;
      }


      molecule[i].unweighted_energy_diff = unw_ext_NB_AA_energy - unw_ext_NB_CG_energy;
      molecule[i].unweighted_virial_diff = unw_ext_NB_AA_virial - unw_ext_NB_CG_virial;

      if(USEFEC == true)
      {
          if(ATOMLAM == true)
          {
                for(l = 0; l < molecule[i].Natoms; l++)
                {
			//cout << "Before, total energy\n";
                    // THIS IS A QUICK HACK FOR POSITION BASED APPLICATION OF THE FEC
		    fecpos = fabs(molecule[i].atom[l].position.x());
		    n = molecule[i].type_number;
		    if(fecpos < molecule[i].box.rmin - 0.5){  // 0.5 and 1.0 is hardcoded here...
			temp_erg = temp_erg - free_energy_comp[0].f[n] - free_energy_comp[0].p[n];
		    }
		    else if(fecpos > molecule[i].box.rmax + 0.5){
			temp_erg = temp_erg - free_energy_comp[fec_bins-1].f[n] - free_energy_comp[fec_bins -1].p[n];
		    }
		    else{
			q = floor((fec_bins - 1)*(fecpos - molecule[i].box.rmin + 0.5)/(molecule[i].box.rmax-molecule[i].box.rmin + 1.0));
			temp_erg = temp_erg - free_energy_comp[q].f[n] - free_energy_comp[q].p[n];
		    }
			//cout << "After, total energy\n";
                    /*lambda = molecule[i].atom[l].lambda;
                    q = floor((fec_bins - 1)*lambda);
                    n = molecule[i].type_number;
                    temp_erg = temp_erg - free_energy_comp[q].f[n] - free_energy_comp[q].p[n];*/
                }
          }
          else
          {
	      cout << "Only atomistic lambda implemented for position based FEC! (total energy)" << endl;
              exit(1);
              /*lambda = molecule[i].lambda;
              q = floor((fec_bins - 1)*lambda);
              n = molecule[i].type_number;
              temp_erg = temp_erg - free_energy_comp[q].f[n] - free_energy_comp[q].p[n];*/
          }
      }

      if(USEPOLYFEC==true)
      {
          lambda = molecule[i].lambda;
          n = molecule[i].type_number;
          temp_erg = temp_erg - poly_fec(n, lambda);
      }

      if(PATHINT == false)
      {
          temp_erg = temp_erg + int_energy + lambda_AA*ext_NB_AA_energy + ext_BB_AA_energy + lambda_CG*ext_NB_CG_energy;
      }
      else
      {
          temp_erg = temp_erg + int_energy + ext_NB_AA_energy + ext_BB_AA_energy;
      }

    }


    tot_energy = temp_erg;

  }





  // CORRECT!
  void mocha_system::calc_gyration()
  {
      int i, k, j;
      double mol_gyr_sqr, com_pos;
      mol_gyr_sqr = 0.0;
      com_pos = 0.0;

      for (i = 0; i < Nmolecules; i++)
      {
          mol_gyr_sqr = 0.0;

          for (k = 0; k < molecule[i].Natoms; k++)
          {
              mol_gyr_sqr += (pbcdist(molecule[i].atom[k].position, molecule[i].com, box)) * (pbcdist(molecule[i].atom[k].position, molecule[i].com, box));
          }
          mol_gyr_sqr = sqrt(mol_gyr_sqr/16.0);

          com_pos = molecule[i].com.getX()+box.hside.getX();
          j = (int)floor(com_pos * gyration->invDelta);
          gyration->counts[j]++;
          gyration->gyr[j] += mol_gyr_sqr;
      }

  }



  void mocha_system::mean_gyration()
  {
      int i, k;
      double temp_gyr_sqr, mol_gyr_sqr;

      temp_gyr_sqr = 0.0;
      mol_gyr_sqr = 0.0;

      for (i = 0; i < Nmolecules; i++)
      {
          mol_gyr_sqr = 0.0;

          for (k = 0; k < molecule[i].Natoms; k++)
          {
              mol_gyr_sqr += (pbcdist(molecule[i].atom[k].position, molecule[i].com, box)) * (pbcdist(molecule[i].atom[k].position, molecule[i].com, box));
          }

          temp_gyr_sqr += sqrt(mol_gyr_sqr/16.0);
      }

      rad_of_gyr = temp_gyr_sqr/(Nmolecules);
  }



  void mocha_system::calc_partQM()
  {
     // Works only for XSLAB at the moment.

     int i, n;

     n = 0;

     for (i = 0; i < Nmolecules; i++)
     {
         /*if (molecule[i].lambda_function( fabs( molecule[i].com.x() ) ) > 0.99999999 )*/
         if (molecule[i].box.rmin - 1.0 >= fabs( molecule[i].com.x() ) ) // CORRECTION WITH 0.5 AS FEC IS APPLIED ALSO IN 0.5nm WIDE PART OF QM REGION
         {
             //cout << molecule[i].com.x() << endl;
             n++;
         }
     }

     partQM = n;
  }



  void mocha_system::kirkwood_energy_difference(double *kirkw, double *virial)
  {

    int i, k;

    for (k = 0; k < Nmol_types; k++)
    {
        kirkw[k] = 0.0;
        virial[k] = 0.0;
    }

    for(i = 0; i < Nmolecules; i++)
    {
        k = molecule[i].type_number;
        kirkw[k] += molecule[i].kirkw_diff;
        virial[k] += molecule[i].virial;
        //cout << "molecule[i].kirkw_diff = " << molecule[i].kirkw_diff << endl;
    }

  }



  void mocha_system::number_of_molecules_per_type(void)
  {

    int i, k;

    molecules_per_type.number = new int[Nmol_types];
    molecules_per_type.atoms = new int[Nmol_types];
    molecules_per_type.type = new string[Nmol_types];

    for (k = 0; k < Nmol_types; k++)
    {
        molecules_per_type.number[k] = 0;
        molecules_per_type.atoms[k] = 0;
    }

    for(i = 0; i < Nmolecules; i++)
    {
        k = molecule[i].type_number;
        molecules_per_type.number[k] += 1;
        molecules_per_type.atoms[k] = molecule[i].Natoms;
        molecules_per_type.type[k] = molecule[i].type;
    }

    for (k = 0; k < Nmol_types; k++) cout << "There are " << molecules_per_type.number[k] << " molecules of type " << molecules_per_type.type[k] << endl;


  }





  void mocha_system::internal_energy(int mol_index, int new_coord, double *energy, double *virial, double *e_diff)
  {
    int i, j, l, lab, lmax, k, atoms;
    double x, vir, r_i, r_k, temp, rzero, rmax, kappa, temp_v, mu_m, lam_i, lam_k, temp_log, m_light, m_heavy, temp_ediff;
    mocha_distance d;
    vec diff, posit_i, posit_k;
    int z;

    atoms = molecule[mol_index].Natoms;

    m_light = 1.0;
    m_heavy = 100.0;

    temp_ediff = 0.0;
    temp = 0.0;
    temp_v = 0.0;
    temp_log = 0.0;

    d.r = 0.0;

    lam_i = 1.0;
    lam_k = 1.0;
    mu_m = 0.0;

    r_i = 0.0;
    r_k = 0.0;

    /** nonbonded part **/

    /*
    WHEN DOING PATH INTEGRAL SIMULATIONS THERE ARE NO INTERNAL NON-BONDED INTERACTIONS

    for(i = 0; i < molecule[mol_index].Natoms-1; i++)
    {
      for(j = i+1; j < molecule[mol_index].Natoms; j++)
      {

	if(new_coord == 0)
	{
            diff = pbcdist(molecule[mol_index].atom[i].position, molecule[mol_index].atom[j].position, box);
	}
	else
	{
            diff = pbcdist(molecule[mol_index].new_atom[i].position, molecule[mol_index].new_atom[j].position, box);
	}

	d.r = abs(diff);


	z = molecule[mol_index].internal_int.isbonded(i, j);

	if(z == 0)
	{

          lab = lj_table[molecule[mol_index].atom[i].interaction_number].row[molecule[mol_index].atom[j].interaction_number];

	  x  = 0.0;
          vir = 0.0;
          energy_lj(d, interaction.ext_nonbonded.lj[lab], &x, &vir);

	  temp += x;
          temp_v += vir;

	}


      }
    }
    */


    /** bonded part **/
    lmax = molecule[mol_index].internal_int.Nbonds;
    if(RIGID == true) lmax = 0;

    for(l = 0; l < lmax; l++)
    {

      i = molecule[mol_index].internal_int.pairs[l].a;
      j = molecule[mol_index].internal_int.pairs[l].b;

      kappa = molecule[mol_index].internal_int.pairs[l].k;
      rzero = molecule[mol_index].internal_int.pairs[l].rzero;
      rmax = molecule[mol_index].internal_int.pairs[l].rmax;

      if(new_coord == 0)
      {
          diff = pbcdist(molecule[mol_index].atom[i].position, molecule[mol_index].atom[j].position, box);
      }
      else
      {
          diff = pbcdist(molecule[mol_index].new_atom[i].position, molecule[mol_index].new_atom[j].position, box);
      }

      x = 0.0;
      vir = 0.0;

      if( (HADRESS_SGNL == true) && (ATOMLAM == true) )
      {
           if(new_coord == 0)
           {
               posit_i = molecule[mol_index].atom[i].position;
           }
           else
           {
               posit_i = molecule[mol_index].new_atom[i].position;
           }

           if(SPHERE_SGNL == true) r_i = abs(posit_i);
           if(XSLAB_SGNL == true) r_i = fabs(posit_i.x());

           lam_i = molecule[mol_index].lambda_function(r_i);
      }

      if( (KIRK_SGNL == true) || (CONST_LAMBDA == true) )
      {
          lam_i = lambda_AA;
      }

      mu_m = m_light*lam_i + (1.0-lam_i)*m_heavy;
      //mu_m = m_heavy;                                                           // GOES OUT AGAIN!!! IMP!!!
      /*
      d.r = fabs( abs(diff) - rzero );
      energy_fene(d, rmax, kappa, &x, &vir);
      */
      d.r = abs(diff);
      //energy_quartic(d, rzero, kappa, &x, &vir);
      energy_quadratic(d, rzero, kappa, &x, &vir);\
        //cout << "x = " << x << endl;
      temp_ediff += x*(m_light - m_heavy);
        //cout << "temp_ediff = " << temp_ediff << endl;
      x *= mu_m;
      vir *= mu_m;
      temp += x;
      temp_v += vir;

    }

    /** Adding the log term - only for Path Integral Simulations **/

    /* For constant-lambda simulation, the log term won't influence the simulation, since it corresponds to a constant global energy shift in this case.
     * Furthermore, when performing adaptive resolution simulations with FEC, we can correct for the log-term exactly and it won't have any influence.
     * Hence, it's easier to remove before directly. However, when performing adaptive resolution simulations without FEC, we need the term to obtain the
     * correct behavior. In this case, uncomment the following routine. */

    /*for(k = 0; k < atoms; k++)
    {

        if( (HADRESS_SGNL == true) && (ATOMLAM == true) )
        {
             if(new_coord == 0)
             {
                 posit_k = molecule[mol_index].atom[k].position;
             }
             else
             {
                 posit_k = molecule[mol_index].new_atom[k].position;
             }

             if(SPHERE_SGNL == true) r_k = abs(posit_k);
             if(XSLAB_SGNL == true) r_k = fabs(posit_k.x());

             lam_k = molecule[mol_index].lambda_function(r_k);
             //if( (lam_i == 0.0) && (lam_j == 0.0) ) continue;
        }

        if( (KIRK_SGNL == true) || (CONST_LAMBDA == true) )
        {
              lam_k = lambda_AA;
        }

        mu_m = m_light*lam_k + (1.0-lam_k)*m_heavy;
        temp_log += log (mu_m);
        temp_ediff += -0.5*(m_light - m_heavy)/mu_m;
    }

    temp_log *= -0.5;
    *energy = temp + temp_log;*/

    //cout << "Final temp = " << temp << endl;
    //cout << "Final temp_ediff = " << temp_ediff << endl;

    *energy = temp;
    *virial = temp_v;
    *e_diff = temp_ediff;

  }














  void mocha_system::bonded_AA_molecular_energy(int mol_index, int new_coord, double *energy, double *virial)
  {

    int i, j, ii, jj, m, n, a, b;
    double temp, x, vir, kappa, rzero, rmax, temp_v;
    vec diff;
    mocha_distance d;

    temp = 0.0;
    temp_v = 0.0;

    n = interaction.ext_bonded.Nmolecular_bonds;

    for(i = 0; i < n; i++)
    {

      m = interaction.ext_bonded.pair_list[i].Nbonds;

      a = interaction.ext_bonded.pair_list[i].a;
      b = interaction.ext_bonded.pair_list[i].b;

      if( (a == mol_index) || (b == mol_index) )
      {

	for(j = 0; j < m; j++)
	{

	  ii = interaction.ext_bonded.pair_list[i].pairs[j].a;
	  jj = interaction.ext_bonded.pair_list[i].pairs[j].b;
	  kappa = interaction.ext_bonded.pair_list[i].pairs[j].k;
	  rzero = interaction.ext_bonded.pair_list[i].pairs[j].rzero;
	  rmax = interaction.ext_bonded.pair_list[i].pairs[j].rmax;


	  if(new_coord == 0)
	  {
	      diff = pbcdist(molecule[a].atom[ii].position, molecule[b].atom[jj].position, box);
	  }
	  else
	  {
	      diff = pbcdist(molecule[a].new_atom[ii].position, molecule[b].atom[jj].position, box);
	      if(b == mol_index) diff = pbcdist(molecule[a].atom[ii].position, molecule[b].new_atom[jj].position, box);
	  }
	  d.r = fabs( abs(diff) - rzero );

          x = 0.0;
          vir = 0.0;
          energy_fene(d, rmax, kappa, &x, &vir);
	  temp += x;
          temp_v += vir;

	}

      }

    }

    /* the returned result V_i is 1/2 of the total energy on the molecule,
       so that \sum_i V_i = 1/2 \sum_ij V_ij  */
    *energy = 0.5 * temp;
    *virial = 0.5 * temp_v;

  }


  void mocha_system::nonbonded_AA_molecular_energy(int ii, int new_coord, double *energy, double *virial, double *unw_energy, double *unw_virial)
  {

    // shared
    int atoms_ii, n;
    double lambda_ii;

    // private
    int i, jj, l, m, lab, atoms_jj;
    double x, v, vir, w, lambda_jj, weight, q1, q2, r_i, r_j, lam_i, lam_j;
    bool boole_ii, boole_jj;
    vec diff, posit_i, posit_j;
    mocha_distance d;

    // reduced
    double temp, temp_v;

    double unw_v, unw_w;

    n = molecule[ii].Nneighbors;

    atoms_ii = molecule[ii].Natoms;

    if(new_coord == 0)
        lambda_ii = molecule[ii].lambda;
    else
        lambda_ii = molecule[ii].new_lambda;

    if(new_coord == 0)
          boole_ii = molecule[ii].ishybrid;
      else
          boole_ii = molecule[ii].new_ishybrid;

    temp = 0.0;
    temp_v = 0.0;
    unw_v = 0.0;
    unw_w = 0.0;

    r_i = 0.0;
    r_j = 0.0;
    lam_i = 0.0;
    lam_j = 0.0;

#ifdef PARALLEL
#pragma omp parallel for \
shared(ii, atoms_ii, n, lambda_ii, boole_ii, new_coord) private(i, jj, l, m, lab, atoms_jj, x, v, vir, w, lambda_jj, weight, boole_jj, diff, d, q1, q2) \
schedule(static) \
reduction(+:temp,temp_v)
#endif

    for(i = 0; i < n; i++)
    {

      jj = molecule[ii].neighbors[i];


      boole_jj = molecule[jj].ishybrid;

      weight = 1.0;
      atoms_jj = molecule[jj].Natoms;

      lambda_jj = molecule[jj].lambda;


      if(KIRK_SGNL == false)
      {
          if( (HADRESS_SGNL == true) && (!boole_ii) && (!boole_jj) && (ATOMLAM == false) )
          {
              continue;
          }

          if( (HADRESS_SGNL == true) && (ATOMLAM == false) )
              weight = 0.5*(lambda_ii + lambda_jj);
          else
              weight = 1.0;
      }
      else
      {
          weight = 1.0;
      }

      v = 0.0;
      w = 0.0;

      for(l = 0; l < atoms_ii; l++)
      {
	for(m = 0; m < atoms_jj; m++)
	{


          lab = lj_table[molecule[ii].atom[l].interaction_number].row[molecule[jj].atom[m].interaction_number];




	  if(new_coord == 0)
	  {
	    diff = pbcdist(molecule[ii].atom[l].position, molecule[jj].atom[m].position, box);
	  }
	  else
	  {
	    diff = pbcdist(molecule[ii].new_atom[l].position, molecule[jj].atom[m].position, box);
	  }



          if( (HADRESS_SGNL == true) && (ATOMLAM == true) )
          {
              lam_i = 1.0;
              lam_j = 1.0;

              if(new_coord == 0)
              {
                  posit_i = molecule[ii].atom[l].position;
              }
              else
              {
                  posit_i = molecule[ii].new_atom[l].position;
              }

              posit_j = molecule[jj].atom[m].position;

              if(SPHERE_SGNL == true) r_i = abs(posit_i);
              if(XSLAB_SGNL == true) r_i = fabs(posit_i.x());
              if(SPHERE_SGNL == true) r_j = abs(posit_j);
              if(XSLAB_SGNL == true) r_j = fabs(posit_j.x());

              lam_i = molecule[ii].lambda_function(r_i);
              lam_j = molecule[ii].lambda_function(r_j);
              if( (lam_i == 0.0) && (lam_j == 0.0) ) continue;
          }

	  d.r = abs(diff);

	  x  = 0.0;
          vir = 0.0;

          energy_lj(d, interaction.ext_nonbonded.lj[lab], &x, &vir);

          if( (HADRESS_SGNL == true) && (ATOMLAM == true) )
          {
              x  = x*0.5*(lam_i + lam_j);
              vir = vir*0.5*(lam_i + lam_j);
          }

          v += x;
          w += vir;

          if(ELECTRO_SGNL == true)
          {
              q1 = molecule[ii].atom[l].charge;
              q2 = molecule[jj].atom[m].charge;
              energy_electric(d, q1, q2, NB_CUTOFF, INV_NB_CUTOFF, BETA, &x, &vir);
              if( (HADRESS_SGNL == true) && (ATOMLAM == true) )
              {
                  x  = x*0.5*(lam_i + lam_j);
                  vir = vir*0.5*(lam_i + lam_j);
              }
              v += x;
              w += vir;
          }

	}
      }

      /* ************ TEST!!! ********** */
      unw_v += v;
      unw_w += w;

      v = weight*v;
      w = weight*w;

      temp += v;
      temp_v += w;

    }

    /* the returned result V_i is 1/2 of the total energy on the molecule,
       so that \sum_i V_i = 1/2 \sum_ij V_ij  */
    temp = 0.5 * temp;
    temp_v = 0.5 * temp_v;


    *energy = temp;
    *virial = temp_v;

    *unw_energy = 0.5*unw_v;
    *unw_virial = 0.5*unw_w;

  }










  void mocha_system::nonbonded_CG_molecular_energy(int ii, int new_coord, double *energy, double *virial, double *unw_energy, double *unw_virial)
  {

    int i, jj, n, lab;
    double x, vir, temp, v, temp_v, w;
    double lambda_ii, lambda_jj, weight;
    double unw_v, unw_w;
    bool boole_ii, boole_jj;
    vec diff;
    mocha_distance d;

    n = molecule[ii].Nneighbors;

    if(new_coord == 0)
        lambda_ii = molecule[ii].lambda;
    else
        lambda_ii = molecule[ii].new_lambda;

    if(new_coord == 0)
          boole_ii = molecule[ii].ishybrid;
      else
          boole_ii = molecule[ii].new_ishybrid;

    temp = 0.0;
    temp_v = 0.0;
    unw_v = 0.0;
    unw_w = 0.0;

#ifdef PARALLEL
#pragma omp parallel for \
shared(ii, n, lambda_ii, boole_ii, new_coord) private(i, jj, lab, x, v, vir, w, lambda_jj, weight, boole_jj, diff, d) \
schedule(static) \
reduction(+:temp,temp_v)
#endif

    for(i = 0; i < n; i++)
    {

      jj = molecule[ii].neighbors[i];

      boole_jj = molecule[jj].ishybrid;

      lambda_jj = molecule[jj].lambda;

      if(KIRK_SGNL == false)
      {
          if( (boole_ii) && (boole_jj) )
          {
              if( (HADRESS_SGNL == true) && (lambda_ii == 1.0) && (lambda_jj == 1.0) ) continue;
          }

          if(HADRESS_SGNL == true)
              weight = 1.0 - 0.5*(lambda_ii + lambda_jj);
          else
              weight = 1.0;
      }
      else
      {
          weight = 1.0;
      }


      lab = lj_table[molecule[ii].interaction_number].row[molecule[jj].interaction_number];

      if(new_coord == 0)
      {
	diff = pbcdist(molecule[ii].com, molecule[jj].com, box);
      }
      else
      {
	diff = pbcdist(molecule[ii].new_com, molecule[jj].com, box);
      }

      d.r = abs(diff);
      x  = 0.0;
      vir = 0.0;
      energy_lj_CG(d, interaction.ext_nonbonded.lj[lab], &x, &vir);


      /* ************ TEST!!! ********** */
      unw_v += x;
      unw_w += vir;

      v = weight*x;
      w = weight*vir;
      temp += v;
      temp_v += w;

    }

    /* the returned result V_i is 1/2 of the total energy on the molecule,
       so that \sum_i V_i = 1/2 \sum_ij V_ij  */
    temp = 0.5 * temp;
    temp_v = 0.5 * temp_v;

    *energy = temp;
    *virial = temp_v;

    *unw_energy = 0.5*unw_v;
    *unw_virial = 0.5*unw_w;

  }



  void mocha_system::nonbonded_PI_molecular_energy(int ii, int new_coord, double *energy, double *virial, double *unw_energy, double *unw_virial)
  {

    // shared
    int atoms_ii, n;
    //double lambda_ii;

    // private
    int i, jj, l, m, lab;//, atoms_jj;
    double x, v, vir, w, r_i, r_j, lam_i, lam_j, lam_ii, lam_jj;//, lambda_jj, weight, q1, q2;
    double x_aa, vir_aa, x_cg, vir_cg, tmp;
    //bool boole_ii, boole_jj;
    vec diff, posit_i, posit_j, posit_ii, posit_jj;
    mocha_distance d;

    // reduced
    double temp, temp_v;

    double unw_v, unw_w;

    n = molecule[ii].Nneighbors;


    atoms_ii = molecule[ii].Natoms;

    /*
    if(new_coord == 0)
        lambda_ii = molecule[ii].lambda;
    else
        lambda_ii = molecule[ii].new_lambda;

    if(new_coord == 0)
          boole_ii = molecule[ii].ishybrid;
      else
          boole_ii = molecule[ii].new_ishybrid;
    */

    temp = 0.0;
    temp_v = 0.0;
    unw_v = 0.0;
    unw_w = 0.0;

    r_i = 0.0;
    r_j = 0.0;
    lam_i = 1.0;
    lam_j = 1.0;
    lam_ii = 1.0;
    lam_jj = 1.0;

#ifdef PARALLEL
#pragma omp parallel for \
shared(ii, atoms_ii, n, lambda_ii, boole_ii, new_coord) private(i, jj, l, m, lab, atoms_jj, x, v, vir, w, lambda_jj, weight, boole_jj, diff, d, q1, q2) \
schedule(static) \
reduction(+:temp,temp_v)
#endif

    for(i = 0; i < n; i++)
    {

      jj = molecule[ii].neighbors[i];
      //atoms_jj = molecule[jj].Natoms;

      /*
      boole_jj = molecule[jj].ishybrid;
      weight = 1.0;
      atoms_jj = molecule[jj].Natoms;
      lambda_jj = molecule[jj].lambda;
      if(KIRK_SGNL == false)
      {
          if( (HADRESS_SGNL == true) && (!boole_ii) && (!boole_jj) && (ATOMLAM == false) )
          {
              continue;
          }

          if( (HADRESS_SGNL == true) && (ATOMLAM == false) )
              weight = 0.5*(lambda_ii + lambda_jj);
          else
              weight = 1.0;
      }
      else
      {
          weight = 1.0;
      }
      */

      v = 0.0;
      w = 0.0;


      if( HADRESS_SGNL == true )
      {
          // CALCULATE MOLECULAR LAMBDA TO GET POSITION
          if(new_coord == 0)
          {
              posit_ii = molecule[ii].com;
          }
          else
          {
              posit_ii = molecule[ii].new_com;
          }

          posit_jj = molecule[jj].com;

          if(SPHERE_SGNL == true) r_i = abs(posit_ii);
          if(XSLAB_SGNL == true) r_i = fabs(posit_ii.x());
          if(SPHERE_SGNL == true) r_j = abs(posit_jj);
          if(XSLAB_SGNL == true) r_j = fabs(posit_jj.x());

          lam_ii = molecule[ii].lambda_function(r_i);
          lam_jj = molecule[ii].lambda_function(r_j);
      }

      if((KIRK_SGNL == true) || (CONST_LAMBDA == true))
      {
          lam_ii = lambda_AA;
          lam_jj = lambda_AA;
      }

      // BOTH IN CLASSICAL REGION? IF YES, GO FOR EFFICIENT FORCE CALCULATION
      if( (lam_ii == 0.0) && (lam_jj == 0.0) ) //&& (HADRESS_SGNL == true))
      {
          //cout << "BOTH LAMBDA ZERO!" << "\n";
          lab = lj_table[molecule[ii].atom[1].interaction_number].row[molecule[jj].atom[1].interaction_number];

          if(new_coord == 0)
          {
            diff = pbcdist(molecule[ii].com, molecule[jj].com, box);
          }
          else
          {
            diff = pbcdist(molecule[ii].new_com, molecule[jj].com, box);
          }

          d.r = abs(diff);

          x  = 0.0;
          vir = 0.0;
          x_cg  = 0.0;
          vir_cg = 0.0;

          //if(AAENERGY == true) energy_lj(d, interaction.ext_nonbonded.lj[lab], &x_aa, &vir_aa);
          //if(AAENERGY == true) energy_SG(d, interaction.ext_nonbonded.lj[lab], BETA, &x_aa, &vir_aa);
          //if(CGENERGY == true)
          if((d.r <= 0.304144687))
              {
                auto t_start1 = std::chrono::high_resolution_clock::now();

                energy_WCA(d, interaction.ext_nonbonded.lj[lab], &x_cg, &vir_cg);
                //energy_IBI_LJ(d, interaction.ext_nonbonded.lj[lab], BETA, &x_cg, &vir_cg);

                auto t_end1 = std::chrono::high_resolution_clock::now();
                EnergyTime += std::chrono::duration<double, std::nano>(t_end1-t_start1).count();
              }

          //if(AAENERGY == true) energy_WCA(d, interaction.ext_nonbonded.lj[lab], &x_aa, &vir_aa);              // for CG simulations
          //if(CGENERGY == true) energy_SG(d, interaction.ext_nonbonded.lj[lab], BETA, &x_cg, &vir_cg);         // for CG simulations

          //if(CGENERGY == true) energy_SG(d, interaction.ext_nonbonded.lj[lab], BETA, &x_cg, &vir_cg);

          x  = x_cg;
          vir = vir_cg;

          /*
          if( (HADRESS_SGNL == true) && (ATOMLAM == true) )
          {
              x  = x*0.5*(lam_i + lam_j);
              vir = vir*0.5*(lam_i + lam_j);
          }
          */

          v += x;
          w += vir;

          unw_v -= x_cg;
          unw_w += vir;

      }
      // BOTH IN CLASSICAL REGION? IF NO, GO AS USUAL
      else
      {
          //cout << "UNEFFICIENT LOOP!" << "\n";
          for(l = 0; l < atoms_ii; l++)
          {
            //for(m = 0; m < atoms_jj; m++)
            //{

              m = l;


              lab = lj_table[molecule[ii].atom[l].interaction_number].row[molecule[jj].atom[m].interaction_number];




              if(new_coord == 0)
              {
                diff = pbcdist(molecule[ii].atom[l].position, molecule[jj].atom[m].position, box);
              }
              else
              {
                diff = pbcdist(molecule[ii].new_atom[l].position, molecule[jj].atom[m].position, box);
              }



              if( (HADRESS_SGNL == true) && (ATOMLAM == true) )
              {
                  lam_i = 1.0;
                  lam_j = 1.0;

                  if(new_coord == 0)
                  {
                      posit_i = molecule[ii].atom[l].position;
                  }
                  else
                  {
                      posit_i = molecule[ii].new_atom[l].position;
                  }

                  posit_j = molecule[jj].atom[m].position;

                  if(SPHERE_SGNL == true) r_i = abs(posit_i);
                  if(XSLAB_SGNL == true) r_i = fabs(posit_i.x());
                  if(SPHERE_SGNL == true) r_j = abs(posit_j);
                  if(XSLAB_SGNL == true) r_j = fabs(posit_j.x());

                  lam_i = molecule[ii].lambda_function(r_i);
                  lam_j = molecule[ii].lambda_function(r_j);
                  //if( (lam_i == 0.0) && (lam_j == 0.0) ) continue;
              }


              if((KIRK_SGNL == true) || (CONST_LAMBDA == true))
              {
                  lam_i = lambda_AA;
                  lam_j = lambda_AA;
              }

              d.r = abs(diff);

              x  = 0.0;
              vir = 0.0;
              x_aa  = 0.0;
              vir_aa = 0.0;
              x_cg  = 0.0;
              vir_cg = 0.0;
              tmp = 0.0;

              //if(AAENERGY == true) energy_lj(d, interaction.ext_nonbonded.lj[lab], &x_aa, &vir_aa);

              if( (lam_i == 1.0) && (lam_j == 1.0) ){
                  //cout << "BOTH LAMBDA ONE!" << "\n";
                  //if(AAENERGY == true)
                  if((d.r <= 0.9))
                  {
                      auto t_start2 = std::chrono::high_resolution_clock::now();

                      energy_SG(d, interaction.ext_nonbonded.lj[lab], BETA, &x_aa, &vir_aa);


                      //if(CGENERGY == true)
                      //energy_WCA(d, interaction.ext_nonbonded.lj[lab], &x_cg, &vir_cg);
                      //energy_IBI_LJ(d, interaction.ext_nonbonded.lj[lab], BETA, &x_cg, &vir_cg);
                      //cout << "if case, x_cg " << x_cg << "\n";

                      auto t_end2 = std::chrono::high_resolution_clock::now();
                      EnergyTime += std::chrono::duration<double, std::nano>(t_end2-t_start2).count();
                  }
              }
              else {
                  //cout << "I AM HYBRID!" << "\n";
                  //if(AAENERGY == true)
                  if((d.r <= 0.9))
                  {
                      auto t_start3 = std::chrono::high_resolution_clock::now();

                      energy_SG(d, interaction.ext_nonbonded.lj[lab], BETA, &x_aa, &vir_aa);
                      //if(CGENERGY == true)
                      energy_WCA(d, interaction.ext_nonbonded.lj[lab], &x_cg, &vir_cg);
                      //energy_IBI_LJ(d, interaction.ext_nonbonded.lj[lab], BETA, &x_cg, &vir_cg);
                      //cout << "else case, x_cg " << x_cg << "\n";

                      auto t_end3 = std::chrono::high_resolution_clock::now();
                      EnergyTime += std::chrono::duration<double, std::nano>(t_end3-t_start3).count();
                  }
              }

              //if(AAENERGY == true) energy_WCA(d, interaction.ext_nonbonded.lj[lab], &x_aa, &vir_aa);              // for CG simulations
              //if(CGENERGY == true) energy_SG(d, interaction.ext_nonbonded.lj[lab], BETA, &x_cg, &vir_cg);         // for CG simulations

              //if(CGENERGY == true) energy_SG(d, interaction.ext_nonbonded.lj[lab], BETA, &x_cg, &vir_cg);

              tmp = 0.5*(lam_i + lam_j);
              //cout << "tmp: " << tmp << "\n";
              x  = (x_aa*tmp + x_cg*(1.0 - tmp))/(double)atoms_ii;
              vir = (vir_aa*tmp + vir_cg*(1.0 - tmp))/(double)atoms_ii;
              //cout << "x: " << x << "\n";
              /*
              if( (HADRESS_SGNL == true) && (ATOMLAM == true) )
              {
                  x  = x*0.5*(lam_i + lam_j);
                  vir = vir*0.5*(lam_i + lam_j);
              }
              */

              v += x;
              w += vir;

              unw_v += (x_aa - x_cg)/(double)atoms_ii;
              unw_w += vir;


              /*
              if(ELECTRO_SGNL == true)
              {
                  q1 = molecule[ii].atom[l].charge;
                  q2 = molecule[jj].atom[m].charge;
                  energy_electric(d, q1, q2, NB_CUTOFF, INV_NB_CUTOFF, BETA, &x, &vir);
                  if( (HADRESS_SGNL == true) && (ATOMLAM == true) )
                  {
                      x  = x*0.5*(lam_i + lam_j);
                      vir = vir*0.5*(lam_i + lam_j);
                  }
                  v += x;
                  w += vir;
              }
              */

            //}
          }

      }

      /* ************ TEST!!! ********** */

      temp += v;
      temp_v += w;

    }

    /* the returned result V_i is 1/2 of the total energy on the molecule,
       so that \sum_i V_i = 1/2 \sum_ij V_ij  */
    temp = 0.5 * temp;
    temp_v = 0.5 * temp_v;


    *energy = temp;
    *virial = temp_v;

    *unw_energy = 0.5*unw_v;
    *unw_virial = 0.5*unw_w;

  }





  void mocha_system::update_neighbors()
  {

    int i, j, ii, jj, l, m, n;
    vec com_i, com_j, diff;
    double d;
    bool quit;
    //string type_i, type_j;




    n = 0;
    for(i = 0; i < Nmolecules; i++)
    {
        molecule[i].Nneighbors = 0;
        com_i = molecule[i].update_com();
    }


#ifdef PARALLEL_NB

#pragma omp parallel for \
private(i, j, ii, jj, l, m, com_i, com_j, diff, d, quit) \
schedule(static)

    for(i = 0; i < Nmolecules; i++)
    {

        ii = molecule[i].Natoms;

        for(j = 0; j < Nmolecules; j++)
        {
            if(i == j) continue;
            //ii = molecule[i].Natoms;
            jj = molecule[j].Natoms;
            quit = false;

            for(l = 0; l < ii; l++)
            {
                for(m = 0; m < jj; m++)
                {
                    com_i = molecule[i].atom[l].position;
                    com_j = molecule[j].atom[m].position;
                    diff = pbcdist(com_i, com_j, box);
                    d = abs(diff);
                    if(d <= NB_CUTOFF)
                    {
                        molecule[i].Nneighbors += 1;
                        n = molecule[i].Nneighbors - 1;
                        molecule[i].neighbors[n] = j;
                        quit = true;
                        break;
                    }
                }
                if(quit == true) break;
            }


        }

    }

#else



    for(i = 0; i < Nmolecules-1; i++)
    {

        ii = molecule[i].Natoms;

        for(j = i+1; j < Nmolecules; j++)
        {
            jj = molecule[j].Natoms;
            quit = false;

            for(l = 0; l < ii; l++)
            {
                for(m = 0; m < jj; m++)
                {
                    com_i = molecule[i].atom[l].position;
                    com_j = molecule[j].atom[m].position;

                    //type_i = molecule[i].atom[l].type;  // These modifications make sure that only particles of the same type are put together as pairs.
                    //type_j = molecule[j].atom[m].type;  // Very useful for Path Integral simulation.

                    diff = pbcdist(com_i, com_j, box);
                    d = abs(diff);
                    if(d <= NB_CUTOFF) //&& (type_i == type_j) )  // Only put particles together if in cutoff distance and if of same type
                    {
                        molecule[i].Nneighbors += 1;
                        molecule[j].Nneighbors += 1;
                        n = molecule[i].Nneighbors - 1;
                        molecule[i].neighbors[n] = j;
                        n = molecule[j].Nneighbors - 1;
                        molecule[j].neighbors[n] = i;
                        quit = true;
                        //cout << "i " << i << " j " << j << " dist " << d << endl;
                        break;
                    }
                }
                if(quit == true) break;
            }


        }

    }

#endif





  }



  void mocha_system::update_neighbors_PI()
  {

    int i, j, ii, jj, l, m, n;
    vec com_i, com_j, diff;
    double d;
    //bool quit;


    n = 0;
    for(i = 0; i < Nmolecules; i++)
    {
        molecule[i].Nneighbors = 0;
        com_i = molecule[i].update_com();
    }


    for(i = 0; i < Nmolecules-1; i++)
    {

        ii = molecule[i].Natoms;

        for(j = i+1; j < Nmolecules; j++)
        {
            jj = molecule[j].Natoms;
            //quit = false;

            for(l = 0; l < ii; l++)
            {
                com_i = molecule[i].atom[l].position;
                com_j = molecule[j].atom[l].position;

                diff = pbcdist(com_i, com_j, box);
                d = abs(diff);
                if(d <= NB_CUTOFF)
                {
                    molecule[i].Nneighbors += 1;
                    molecule[j].Nneighbors += 1;
                    n = molecule[i].Nneighbors - 1;
                    molecule[i].neighbors[n] = j;
                    n = molecule[j].Nneighbors - 1;
                    molecule[j].neighbors[n] = i;
                    //quit = true;
                    //cout << "i " << i << " j " << j << " dist " << d << endl;
                    break;
                }
            }


        }

    }

  }



  void mocha_system::com_shift_center()
  {

    int i, j, atoms;
    vec temp_av;

    temp_av = vec(0.0, 0.0, 0.0);

    for(i = 0; i < Nmolecules; i++)
    {
      atoms = molecule[i].Natoms;

      for(j = 0; j < atoms; j++)
      {
	temp_av += molecule[i].atom[j].position;
      }
    }

    temp_av = temp_av/Natoms;


    for(i = 0; i < Nmolecules; i++)
    {
      atoms = molecule[i].Natoms;

      for(j = 0; j < atoms; j++)
      {
	molecule[i].atom[j].position = pbcdist(molecule[i].atom[j].position, temp_av, box);

	if((molecule[i].atom[j].position.x() > box.hside.x()) || (molecule[i].atom[j].position.x() < -box.hside.x()) ) cout << "ARGH! X\n";
	if((molecule[i].atom[j].position.y() > box.hside.y()) || (molecule[i].atom[j].position.y() < -box.hside.y()) ) cout << "ARGH! Y\n";
	if((molecule[i].atom[j].position.z() > box.hside.z()) || (molecule[i].atom[j].position.z() < -box.hside.z()) ) cout << "ARGH! Z\n";
      }
    }

  }





  void mocha_system::com_shift_corner()
  {

    int i, j, atoms;

    for(i = 0; i < Nmolecules; i++)
    {
      atoms = molecule[i].Natoms;

      for(j = 0; j < atoms; j++)
      {
	molecule[i].atom[j].position += box.side;
      }
    }

  }







/*************/
  /* this one is NOT used!!! */
  /*
    void mocha_system::delta_energy(int i, double *delta_erg, double *int_energy, double *ext_energy, double *kirkw_diff, double *virial)
    {

        double int_erg, ext_NB_AA_energy, ext_BB_AA_energy, ext_NB_CG_energy, delta;
        double int_vir, ext_NB_AA_virial, ext_BB_AA_virial, ext_NB_CG_virial;
        double old_lambda, new_lambda;
        double unw_ext_NB_AA_energy, unw_ext_NB_AA_virial, unw_ext_NB_CG_energy, unw_ext_NB_CG_virial;
        int old_q, new_q, n;

        int_erg = 0.0;
        ext_NB_AA_energy = 0.0;
        ext_BB_AA_energy = 0.0;
        ext_NB_CG_energy = 0.0;

        int_vir = 0.0;
        ext_NB_AA_virial = 0.0;
        ext_BB_AA_virial = 0.0;
        ext_NB_CG_virial = 0.0;

        unw_ext_NB_AA_energy = 0.0;
        unw_ext_NB_AA_virial = 0.0;
        unw_ext_NB_CG_energy = 0.0;
        unw_ext_NB_CG_virial = 0.0;

        internal_energy(i, 1, &int_erg, &int_vir);
        if(AAENERGY == true) nonbonded_AA_molecular_energy(i, 1, &ext_NB_AA_energy, &ext_NB_AA_virial, &unw_ext_NB_AA_energy, &unw_ext_NB_AA_virial);
        bonded_AA_molecular_energy(i, 1, &ext_BB_AA_energy, &ext_BB_AA_virial);
        if(CGENERGY == true) nonbonded_CG_molecular_energy(i, 1, &ext_NB_CG_energy, &ext_NB_CG_virial, &unw_ext_NB_CG_energy, &unw_ext_NB_CG_virial);

        delta = int_erg - molecule[i].int_energy + 2.0 * (lambda_AA*ext_NB_AA_energy + ext_BB_AA_energy + lambda_CG*ext_NB_CG_energy - molecule[i].ext_energy);

        *int_energy = int_erg;
        *ext_energy = lambda_AA*ext_NB_AA_energy + ext_BB_AA_energy + lambda_CG*ext_NB_CG_energy;
        *kirkw_diff = ext_NB_AA_energy - ext_NB_CG_energy;
        *virial = lambda_AA*ext_NB_AA_virial + int_vir + ext_BB_AA_virial + lambda_CG*ext_NB_CG_virial;


        if(USEFEC == true)
        {
            old_lambda = molecule[i].lambda;
            new_lambda = molecule[i].new_lambda;

            if(fabs(new_lambda - old_lambda)>0.0)
            {
                old_q = floor((fec_bins - 1)*old_lambda);
                new_q = floor((fec_bins - 1)*new_lambda);
                if(new_q!=old_q)
                {
                    n = molecule[i].type_number;
                    delta = delta - ( free_energy_comp[new_q].f[n] - free_energy_comp[old_q].f[n]) - ( free_energy_comp[new_q].p[n] - free_energy_comp[old_q].p[n]);
                }
            }
        }

        if(USEPOLYFEC==true)
        {

            old_lambda = molecule[i].lambda;
            new_lambda = molecule[i].new_lambda;

            if(fabs(new_lambda - old_lambda)>0.0)
            {
                n = molecule[i].type_number;
                delta = delta - (poly_fec(n, new_lambda) - poly_fec(n, old_lambda));
            }

        }



        *delta_erg = delta;


    }
    */

    data_list_struct mocha_system::delta_energy_function(int i)
    {
        double e_diff;
        double int_erg, ext_NB_AA_energy, ext_BB_AA_energy, ext_NB_CG_energy, delta;
        double int_vir, ext_NB_AA_virial, ext_BB_AA_virial, ext_NB_CG_virial;
        double old_lambda, new_lambda, old_fecpos, new_fecpos;
        double unw_ext_NB_AA_energy, unw_ext_NB_AA_virial, unw_ext_NB_CG_energy, unw_ext_NB_CG_virial;
        int old_q, new_q, n, l;
        data_list_struct data_list;

        e_diff = 0.0;
        int_erg = 0.0;
        ext_NB_AA_energy = 0.0;
        ext_BB_AA_energy = 0.0;
        ext_NB_CG_energy = 0.0;

        int_vir = 0.0;
        ext_NB_AA_virial = 0.0;
        ext_BB_AA_virial = 0.0;
        ext_NB_CG_virial = 0.0;

        unw_ext_NB_AA_energy = 0.0;
        unw_ext_NB_AA_virial = 0.0;
        unw_ext_NB_CG_energy = 0.0;
        unw_ext_NB_CG_virial = 0.0;

        internal_energy(i, 1, &int_erg, &int_vir, &e_diff);
        if((AAENERGY == true) && (PATHINT == false)) nonbonded_AA_molecular_energy(i, 1, &ext_NB_AA_energy, &ext_NB_AA_virial, &unw_ext_NB_AA_energy, &unw_ext_NB_AA_virial);
        bonded_AA_molecular_energy(i, 1, &ext_BB_AA_energy, &ext_BB_AA_virial);
        if((CGENERGY == true) && (PATHINT == false)) nonbonded_CG_molecular_energy(i, 1, &ext_NB_CG_energy, &ext_NB_CG_virial, &unw_ext_NB_CG_energy, &unw_ext_NB_CG_virial);
        if(PATHINT == true) nonbonded_PI_molecular_energy(i, 1, &ext_NB_AA_energy, &ext_NB_AA_virial, &unw_ext_NB_AA_energy, &unw_ext_NB_AA_virial);



        if(PATHINT == false)
        {
            delta = int_erg - molecule[i].int_energy + 2.0 * (lambda_AA*ext_NB_AA_energy + ext_BB_AA_energy + lambda_CG*ext_NB_CG_energy - molecule[i].ext_energy);
            data_list.x_1 = int_erg;
            data_list.x_2 = lambda_AA*ext_NB_AA_energy + ext_BB_AA_energy + lambda_CG*ext_NB_CG_energy;
            data_list.x_3 = ext_NB_AA_energy - ext_NB_CG_energy;
            data_list.x_4 = lambda_AA*ext_NB_AA_virial + int_vir + ext_BB_AA_virial + lambda_CG*ext_NB_CG_virial;
        }
        else
        {
            delta = int_erg - molecule[i].int_energy + 2.0 * (ext_NB_AA_energy + ext_BB_AA_energy - molecule[i].ext_energy);
            data_list.x_1 = int_erg;
            data_list.x_2 = ext_NB_AA_energy + ext_BB_AA_energy;
            data_list.x_3 = unw_ext_NB_AA_energy + e_diff;
            data_list.x_4 = unw_ext_NB_AA_virial + int_vir + ext_BB_AA_virial;
            //cout << "NB energies " << 2.0 * (ext_NB_AA_energy + ext_BB_AA_energy - molecule[i].ext_energy) << endl;
            //cout << "        int_erg before " << molecule[i].int_energy << endl;
            //cout << "        int_erg after " << int_erg << endl;
            //cout << "delta1 " << delta << endl;
            //cout << "unw_ext_NB_AA_energy = " << unw_ext_NB_AA_energy << endl;
            //cout << "e_diff = " << e_diff << endl;
        }

        if(USEFEC == true)
        {
            if(ATOMLAM == true)
            {
                for(l = 0; l < molecule[i].Natoms; l++)
                {
			//cout << "Before, delta energy\n";
		    // THIS IS A QUICK HACK FOR POSITION BASED APPLICATION OF THE FEC
                    old_fecpos = fabs(molecule[i].atom[l].position.x());
                    new_fecpos = fabs(molecule[i].new_atom[l].position.x());

		    // Checks for old pos
                    if(old_fecpos < molecule[i].box.rmin - 0.5){  // 0.5 and 1.0 is hardcoded here...
                        //temp_erg = temp_erg - free_energy_comp[0].f[n] - free_energy_comp[0].p[n];
			old_q = 0;
                    }
                    else if(old_fecpos > molecule[i].box.rmax + 0.5){
                        //temp_erg = temp_erg - free_energy_comp[fec_bins-1].f[n] - free_energy_comp[fec_bins -1].p[n];
                        old_q = fec_bins -1;
		    }
                    else{
                        old_q = floor((fec_bins - 1)*(old_fecpos - molecule[i].box.rmin + 0.5)/(molecule[i].box.rmax-molecule[i].box.rmin + 1.0));
                        //temp_erg = temp_erg - free_energy_comp[q].f[n] - free_energy_comp[q].p[n];
                    }
		    // Checks for new pos
                    if(new_fecpos < molecule[i].box.rmin - 0.5){  // 0.5 and 1.0 is hardcoded here...
                        //temp_erg = temp_erg - free_energy_comp[0].f[n] - free_energy_comp[0].p[n];
                        new_q = 0;
                    }
                    else if(new_fecpos > molecule[i].box.rmax + 0.5){
                        //temp_erg = temp_erg - free_energy_comp[fec_bins-1].f[n] - free_energy_comp[fec_bins -1].p[n];
                        new_q = fec_bins -1;
                    }
                    else{
                        new_q = floor((fec_bins - 1)*(new_fecpos - molecule[i].box.rmin + 0.5)/(molecule[i].box.rmax-molecule[i].box.rmin + 1.0));
                        //temp_erg = temp_erg - free_energy_comp[q].f[n] - free_energy_comp[q].p[n];
                    }

                    if(new_q!=old_q)
                    {
                        n = molecule[i].type_number;
                        delta = delta - ( free_energy_comp[new_q].f[n] - free_energy_comp[old_q].f[n]) - ( free_energy_comp[new_q].p[n] - free_energy_comp[old_q].p[n]);

                        //cout << "    DELTA = " << ( free_energy_comp[new_q].f[n] - free_energy_comp[old_q].f[n]) - ( free_energy_comp[new_q].p[n] - free_energy_comp[old_q].p[n]) << endl;
                        //cout << "    free_energy_comp[new_q].f[n] " << free_energy_comp[new_q].f[n] << endl;
                        //cout << "    free_energy_comp[old_q].f[n] " << free_energy_comp[old_q].f[n] << endl;
                        //cout << "    free_energy_comp[new_q].p[n] " << free_energy_comp[new_q].p[n] << endl;
                        //cout << "    free_energy_comp[old_q].p[n] " << free_energy_comp[old_q].p[n] << endl;
                    }

			//cout << "After, delta energy\n";

                    /*old_lambda = molecule[i].atom[l].lambda;
                    new_lambda = molecule[i].new_atom[l].lambda;

                    if(fabs(new_lambda - old_lambda)>0.0)
                    {
                        old_q = floor((fec_bins - 1)*old_lambda);
                        new_q = floor((fec_bins - 1)*new_lambda);
                        if(new_q!=old_q)
                        {
                            n = molecule[i].type_number;
                            delta = delta - ( free_energy_comp[new_q].f[n] - free_energy_comp[old_q].f[n]) - ( free_energy_comp[new_q].p[n] - free_energy_comp[old_q].p[n]);

                            //cout << "    DELTA = " << ( free_energy_comp[new_q].f[n] - free_energy_comp[old_q].f[n]) - ( free_energy_comp[new_q].p[n] - free_energy_comp[old_q].p[n]) << endl;
                            //cout << "    free_energy_comp[new_q].f[n] " << free_energy_comp[new_q].f[n] << endl;
                            //cout << "    free_energy_comp[old_q].f[n] " << free_energy_comp[old_q].f[n] << endl;
                            //cout << "    free_energy_comp[new_q].p[n] " << free_energy_comp[new_q].p[n] << endl;
                            //cout << "    free_energy_comp[old_q].p[n] " << free_energy_comp[old_q].p[n] << endl;
                        }
                    }*/
                }
            }
            else
            {
		cout << "Only atomistic lambda implemented for position based FEC! (delta energy)" << endl;
                exit(1);
                /*old_lambda = molecule[i].lambda;
                new_lambda = molecule[i].new_lambda;

                if(fabs(new_lambda - old_lambda)>0.0)
                {
                    old_q = floor((fec_bins - 1)*old_lambda);
                    new_q = floor((fec_bins - 1)*new_lambda);
                    if(new_q!=old_q)
                    {
                        n = molecule[i].type_number;
                        delta = delta - ( free_energy_comp[new_q].f[n] - free_energy_comp[old_q].f[n]) - ( free_energy_comp[new_q].p[n] - free_energy_comp[old_q].p[n]);
                    }
                }*/
            }
        }

        if(USEPOLYFEC==true)
        {

            old_lambda = molecule[i].lambda;
            new_lambda = molecule[i].new_lambda;

            if(fabs(new_lambda - old_lambda)>0.0)
            {
                n = molecule[i].type_number;
                delta = delta - (poly_fec(n, new_lambda) - poly_fec(n, old_lambda));
            }

        }



        data_list.x_0 = delta;
        //cout << "delta2 " << delta << endl;

        return data_list;


    }



    /*void mocha_system::sweep()
    {

        int i, j, jj, NatomsLocal, l, m, q;
        double delta, x;
        data_list_struct data_list;
        NatomsLocal = molecule[0].Natoms;
        vec coord, new_coord, displ;

        for (jj = 0; jj < NatomsLocal*Nmolecules; jj++)
        {

            //ATTEMPTS += 1;
            j = distNmol(mt);
            //j = rand() % Nmolecules;
            //cout << "Molecule j chosen: j = " << j << endl;
            molecule[j].displace();
            //cout << "MoveCase " << molecule[j].MoveCase << endl;

            data_list.clear();
            data_list = delta_energy_function(j);
            //x = rand() / double(RAND_MAX);
            x = ddist2(mt);
            delta = data_list.x_0;

            if (molecule[j].MoveCase == 0)
                {ATTEMPTS[0] += 1;}
            else if (molecule[j].MoveCase == 1)
                {ATTEMPTS[1] += 1;}
            else if (molecule[j].MoveCase == 2)
                {ATTEMPTS[2] += 1;}
            else if (molecule[j].MoveCase == 3)
                {ATTEMPTS[3] += 1;}
            else
                {cout << "MoveCase fail... MoveCase = " << molecule[j].MoveCase << endl;
                exit(1);}

            if (x < exp(-delta))
            {
                //cout << "ACCEPTED!" << endl;
                //SUCCESSES += 1;

                //molecule[j].accept();
                //molecule[j].int_energy = data_list.x_1;
                //molecule[j].ext_energy = data_list.x_2;
                //molecule[j].kirkw_diff = data_list.x_3;
                //molecule[j].virial = data_list.x_4;

                if (molecule[j].MoveCase == 0)
                {
                    SUCCESSES[0] += 1;
                    i = molecule[j].MovedAtom;                                //REMOVE
                    molecule[j].atom[i].dist_trav += molecule[j].Move;

                    if (abs(molecule[j].atom[i].dist_trav) > NB_SKIN_HALF)
                    {
                        UpdateFlag = true;
                    }                                                         //REMOVE

                    /*for (q = 0; q < NatomsLocal; q++)
                    {
                        //coord = pbcdist(molecule[j].atom[q].position, molecule[j].com, box);
                        //new_coord = molecule[j].rotate(molecule[j].Move, coord, molecule[j].angle);
                        //displ = pbcdist(new_coord, coord, box);
                        displ = pbcdist(molecule[j].new_atom[q].position, molecule[j].atom[q].position, box);
                        //cout << "atomistic displacement, rotation " << displ << endl;
                        molecule[j].atom[q].dist_trav += displ;

                        if (abs(molecule[j].atom[q].dist_trav) > NB_SKIN_HALF)
                        {
                            UpdateFlag = true;
                            break;
                        }
                    }*/
                /*}
                else if (molecule[j].MoveCase == 1)
                {
                    SUCCESSES[1] += 1;
                    for (q = 0; q < NatomsLocal; q++)
                    {
                        molecule[j].atom[q].dist_trav += molecule[j].Move;
                        //cout << "atomistic displacement, CoM move " << molecule[j].Move << endl;
                        if (abs(molecule[j].atom[q].dist_trav) > NB_SKIN_HALF)
                        {
                            UpdateFlag = true;
                            break;
                        }
                    }
                }
                else if (molecule[j].MoveCase == 2)
                {
                    SUCCESSES[2] += 1;
                    for (q = 0; q < NatomsLocal; q++)
                    {
                        //coord = pbcdist(molecule[j].atom[q].position, molecule[j].com, box);
                        //new_coord = molecule[j].rotate(molecule[j].Move, coord, molecule[j].angle);
                        //displ = pbcdist(new_coord, coord, box);
                        displ = pbcdist(molecule[j].new_atom[q].position, molecule[j].atom[q].position, box);
                        //cout << "atomistic displacement, rotation " << displ << endl;
                        molecule[j].atom[q].dist_trav += displ;

                        if (abs(molecule[j].atom[q].dist_trav) > NB_SKIN_HALF)
                        {
                            UpdateFlag = true;
                            break;
                        }
                    }
                }
                else if (molecule[j].MoveCase == 3)
                {
                    SUCCESSES[3] += 1;
                    for (q = 0; q < NatomsLocal; q++)
                    {
                        displ = pbcdist(molecule[j].new_atom[q].position, molecule[j].atom[q].position, box);
                        //cout << "atomistic displacement, collapse/expand " << displ << endl;
                        molecule[j].atom[q].dist_trav += displ;

                        if (abs(molecule[j].atom[q].dist_trav) > NB_SKIN_HALF)
                        {
                            UpdateFlag = true;
                            break;
                        }
                    }
                }
                else
                {
                    cout << "ERROR in sweep routine! Wrong MoveCase" << endl;
                }

                molecule[j].accept();
                molecule[j].int_energy = data_list.x_1;
                molecule[j].ext_energy = data_list.x_2;
                molecule[j].kirkw_diff = data_list.x_3;
                molecule[j].virial = data_list.x_4;

                if (UpdateFlag == true)
                {
                    break;
                }

            }
        }

        if (UpdateFlag == true)
        {
            for(m = 0; m < Nmolecules; m++)
            {
                for (l = 0; l < NatomsLocal; l++)
                {
                        molecule[m].atom[l].dist_trav = vec(0.0,0.0,0.0);
                }
            }
        }

    }*/



    void mocha_system::sweep()
    {

        int i, j, jj, NatomsLocal, l, m, q;
        double delta, x;
        data_list_struct data_list;
        NatomsLocal = molecule[0].Natoms;
        vec coord, new_coord, displ;

        for (jj = 0; jj < NatomsLocal*Nmolecules; jj++)
        {

            //ATTEMPTS += 1;
            j = distNmol(mt);
            //j = rand() % Nmolecules;
            //cout << "Molecule j chosen: j = " << j << endl;
            molecule[j].displace();
            //cout << "MoveCase " << molecule[j].MoveCase << endl;

            data_list.clear();
            data_list = delta_energy_function(j);
            //x = rand() / double(RAND_MAX);
            x = ddist2(mt);
            delta = data_list.x_0;

            if (molecule[j].MoveCase == 0)
                {ATTEMPTS[0] += 1;}
            else if (molecule[j].MoveCase == 1)
                {ATTEMPTS[1] += 1;}
            else if (molecule[j].MoveCase == 2)
                {ATTEMPTS[2] += 1;}
            else if (molecule[j].MoveCase == 3)
                {ATTEMPTS[3] += 1;}
            else
                {cout << "MoveCase fail... MoveCase = " << molecule[j].MoveCase << endl;
                exit(1);}

            if (x < exp(-delta))
            {
                //cout << "ACCEPTED!" << endl;
                //SUCCESSES += 1;

                //molecule[j].accept();
                //molecule[j].int_energy = data_list.x_1;
                //molecule[j].ext_energy = data_list.x_2;
                //molecule[j].kirkw_diff = data_list.x_3;
                //molecule[j].virial = data_list.x_4;

                if (molecule[j].MoveCase == 0)
                {
                    SUCCESSES[0] += 1;
                    i = molecule[j].MovedAtom;                                //REMOVE
                    molecule[j].atom[i].dist_trav += molecule[j].Move;

                    if (abs(molecule[j].atom[i].dist_trav) > NB_SKIN_HALF)
                    {
                        UpdateFlag = true;
                    }                                                         //REMOVE

                    /*for (q = 0; q < NatomsLocal; q++)
                    {
                        //coord = pbcdist(molecule[j].atom[q].position, molecule[j].com, box);
                        //new_coord = molecule[j].rotate(molecule[j].Move, coord, molecule[j].angle);
                        //displ = pbcdist(new_coord, coord, box);
                        displ = pbcdist(molecule[j].new_atom[q].position, molecule[j].atom[q].position, box);
                        //cout << "atomistic displacement, rotation " << displ << endl;
                        molecule[j].atom[q].dist_trav += displ;

                        if (abs(molecule[j].atom[q].dist_trav) > NB_SKIN_HALF)
                        {
                            UpdateFlag = true;
                            break;
                        }
                    }*/
                }
                else if (molecule[j].MoveCase == 1)
                {
                    SUCCESSES[1] += 1;
                    for (q = 0; q < NatomsLocal; q++)
                    {
                        molecule[j].atom[q].dist_trav += molecule[j].Move;
                        //cout << "atomistic displacement, CoM move " << molecule[j].Move << endl;
                        if (abs(molecule[j].atom[q].dist_trav) > NB_SKIN_HALF)
                        {
                            UpdateFlag = true;
                            break;
                        }
                    }
                }
                else if (molecule[j].MoveCase == 2)
                {
                    SUCCESSES[2] += 1;
                    for (q = 0; q < NatomsLocal; q++)
                    {
                        //coord = pbcdist(molecule[j].atom[q].position, molecule[j].com, box);
                        //new_coord = molecule[j].rotate(molecule[j].Move, coord, molecule[j].angle);
                        //displ = pbcdist(new_coord, coord, box);
                        displ = pbcdist(molecule[j].new_atom[q].position, molecule[j].atom[q].position, box);
                        //cout << "atomistic displacement, rotation " << displ << endl;
                        molecule[j].atom[q].dist_trav += displ;

                        if (abs(molecule[j].atom[q].dist_trav) > NB_SKIN_HALF)
                        {
                            UpdateFlag = true;
                            break;
                        }
                    }
                }
                else if (molecule[j].MoveCase == 3)
                {
                    SUCCESSES[3] += 1;
                    for (q = 0; q < NatomsLocal; q++)
                    {
                        displ = pbcdist(molecule[j].new_atom[q].position, molecule[j].atom[q].position, box);
                        //cout << "atomistic displacement, collapse/expand " << displ << endl;
                        molecule[j].atom[q].dist_trav += displ;

                        if (abs(molecule[j].atom[q].dist_trav) > NB_SKIN_HALF)
                        {
                            UpdateFlag = true;
                            break;
                        }
                    }
                }
                else
                {
                    cout << "ERROR in sweep routine! Wrong MoveCase" << endl;
                }

                molecule[j].accept();
                molecule[j].int_energy = data_list.x_1;
                molecule[j].ext_energy = data_list.x_2;
                molecule[j].kirkw_diff = data_list.x_3;
                molecule[j].virial = data_list.x_4;

            }

            if (UpdateFlag == true)
            {
                update_neighbors_PI();
                for(m = 0; m < Nmolecules; m++)
                {
                    for (l = 0; l < NatomsLocal; l++)
                    {
                            molecule[m].atom[l].dist_trav = vec(0.0,0.0,0.0);
                    }
                }
                VerletListUpdates += 1;
                UpdateFlag = false;
            }

        }

    }



    /* THIS IS NEVER USED */
    void mocha_system::parallel_sweep()
    {

        int j, jj;
        double delta, x;

        int mol_a, mol_b, i, list[2];
        double dist;
        data_list_struct data_list[2];

        for (jj = 0; jj < Nmolecules; jj++)
        {

            list[0] = 0;
            list[1] = 0;

            if(NESTED==0)
            {
                mol_a = rand() % Nmolecules;
                mol_b = -1;
            }
            else
            {
                do
                {
                    mol_a = rand() % Nmolecules;
                    mol_b = rand() % Nmolecules;
                    dist = abs(pbcdist(molecule[mol_a].com, molecule[mol_b].com, box));
                }while((mol_a==mol_b)||(dist < 2.0*NB_CUTOFF));
            }

            list[0] = mol_a;
            list[1] = mol_b;

#ifdef PARALLEL_SWEEP
#pragma omp parallel for shared(list) private(i, j) schedule(static)
#endif
            for(i = 0; i <= NESTED; i++)
            {

                j = list[i];
                molecule[j].displace();

                data_list[i].clear();
                data_list[i] = delta_energy_function(j);
                x = rand() / double(RAND_MAX);

                delta = data_list[i].x_0;

                if (x < exp(-delta))
                {
                    molecule[j].accept();
                    molecule[j].int_energy = data_list[i].x_1;
                    molecule[j].ext_energy = data_list[i].x_2;
                    molecule[j].kirkw_diff = data_list[i].x_3;
                    molecule[j].virial = data_list[i].x_4;
                }
            }

        }



    }



    void mocha_system::initialize_gyration_structure(void)
    {
        int i;

        gyration = new gyration_struct;

        gyration->bins = GYR_BINS;
        gyration->invDelta = gyration->bins / box.side.getX();
        gyration->counts = new int [gyration->bins];
        gyration->gyr = new double [gyration->bins];

        for(i = 0; i < gyration->bins; i++)
        {
            gyration->counts[i] = 0.0;
            gyration->gyr[i] = 0.0;
        }
    }



    void mocha_system::initialize_density_structure(void)
    {

        int i, j, t;

        density = new density_struct [Nthreads];
        FEC_density = new density_struct [1];

        for(t = 0; t < Nthreads; t++)
        {

            density[t].Xbins = XDENS;
            density[t].Ybins = YDENS;
            density[t].Zbins = ZDENS;
            density[t].Rbins = RDENS;

            density[t].invXdelta = density[t].Xbins / box.side.getX();
            density[t].invYdelta = density[t].Ybins / box.side.getY();
            density[t].invZdelta = density[t].Zbins / box.side.getZ();
            density[t].invRdelta = density[t].Rbins / box.min_side;

            density[t].dens_vector = new density_vector_struct [Nmol_types];

            for(i = 0; i < Nmol_types; i++)
            {
                density[t].dens_vector[i].Xdens = new int [density[t].Xbins];
                density[t].dens_vector[i].Ydens = new int [density[t].Ybins];
                density[t].dens_vector[i].Zdens = new int [density[t].Zbins];
                density[t].dens_vector[i].Rdens = new int [density[t].Rbins];
            }

            for(j = 0; j < Nmol_types; j++)
            {
                for (i = 0; i < density[t].Xbins; i++) density[t].dens_vector[j].Xdens[i] = 0;
                for (i = 0; i < density[t].Ybins; i++) density[t].dens_vector[j].Ydens[i] = 0;
                for (i = 0; i < density[t].Zbins; i++) density[t].dens_vector[j].Zdens[i] = 0;
                for (i = 0; i < density[t].Rbins; i++) density[t].dens_vector[j].Rdens[i] = 0;
            }

        }

        FEC_density[0].Xbins = XDENS;
        FEC_density[0].invXdelta = FEC_density[0].Xbins / box.side.getX();
        FEC_density[0].dens_vector = new density_vector_struct [Nmol_types];
        for(i = 0; i < Nmol_types; i++) FEC_density[0].dens_vector[i].Xdens = new int [FEC_density[0].Xbins];
        for(j = 0; j < Nmol_types; j++) for (i = 0; i < FEC_density[0].Xbins; i++) FEC_density[0].dens_vector[j].Xdens[i] = 0;
    }



  void mocha_system::calc_density(int direction)
  {


      int i, j, t, n;
      double d;



#ifdef PARALLEL_DENS
#pragma omp parallel for shared(direction) private(i, j, t, d, n) schedule(static)
#endif
      for(i = 0; i < Nmolecules; i++)
      {

          //t = omp_get_thread_num();
          t = 0;

          n = molecule[i].type_number;

          switch ( direction )
          {
              case 0:
                  d = molecule[i].com.getX()+box.hside.getX();
                  j = (int)floor(d * density[t].invXdelta);
                  density[t].dens_vector[n].Xdens[j]++;
                  FEC_density[0].dens_vector[n].Xdens[j]++;
                  break;
              case 1:
                  d = molecule[i].com.getY()+box.hside.getY();
                  j = (int)floor(d * density[t].invYdelta);
                  density[t].dens_vector[n].Ydens[j]++;
                  break;
              case 2:
                  d = molecule[i].com.getZ()+box.hside.getZ();
                  j = (int)floor(d * density[t].invZdelta);
                  density[t].dens_vector[n].Zdens[j]++;
                  break;
              case 3:
                  d = abs(molecule[i].com);
                  if(d <= box.min_hside)
                  {
                      j = (int)floor(d * density[t].invRdelta);
                      density[t].dens_vector[n].Rdens[j]++;
                  }
                  break;
          }

      }




  }




  void mocha_system::print_density(int ticket)
  {

      int i, j, t;


      for(t = 1; t < Nthreads; t++)
      {
          for(j = 0; j < Nmol_types; j++)
          {
              for(i = 0; i < density[0].Xbins; i++) density[0].dens_vector[j].Xdens[i] += density[t].dens_vector[j].Xdens[i];
              for(i = 0; i < density[0].Ybins; i++) density[0].dens_vector[j].Ydens[i] += density[t].dens_vector[j].Ydens[i];
              for(i = 0; i < density[0].Zbins; i++) density[0].dens_vector[j].Zdens[i] += density[t].dens_vector[j].Zdens[i];
              for(i = 0; i < density[0].Rbins; i++) density[0].dens_vector[j].Rdens[i] += density[t].dens_vector[j].Rdens[i];
          }
      }



      stringstream name_dens;
      name_dens << "density_" << ticket << ".dat";
      ofstream densFile(name_dens.str().c_str(), ios::out);

      for(i = 0; i < density[0].Xbins; i++)
      {
          densFile << i / density[0].invXdelta << " ";
          for(j = 0; j < Nmol_types; j++) densFile << density[0].dens_vector[j].Xdens[i] / ((double)OUTER_SWEEPS * (double)Nmolecules / (double)density[0].Xbins) << " ";
          densFile  << endl;
      }

      densFile  << endl;
      for(i = 0; i < density[0].Ybins; i++)
      {
          densFile << i / density[0].invYdelta << " ";
          for(j = 0; j < Nmol_types; j++) densFile << density[0].dens_vector[j].Ydens[i] / ((double)OUTER_SWEEPS * (double)Nmolecules / (double)density[0].Ybins) << " ";
          densFile  << endl;
      }

      densFile  << endl;
      for(i = 0; i < density[0].Zbins; i++)
      {
          densFile << i / density[0].invZdelta << " ";
          for(j = 0; j < Nmol_types; j++) densFile << density[0].dens_vector[j].Zdens[i] / ((double)OUTER_SWEEPS * (double)Nmolecules / (double)density[0].Zbins) << " ";
          densFile  << endl;
      }

      densFile  << endl;
      for(i = 0; i < density[0].Rbins; i++)
      {
          densFile << i / density[0].invRdelta << " ";
          for(j = 0; j < Nmol_types; j++) densFile << density[0].dens_vector[j].Rdens[i] / ((double)OUTER_SWEEPS * (4.0 / 3.0) * PI * (double)Nmolecules * ( 3.0*i*i + 3.0*i + 1.0) / box.volume ) << " ";
          densFile  << endl;
      }

      densFile.close();

  }



  void mocha_system::print_gyration(int ticket)
  {

      int i, j, t;



      //for(i = 0; i < density[0].Xbins; i++) density[0].dens_vector[j].Xdens[i] += density[t].dens_vector[j].Xdens[i];





      stringstream name_gyr;
      name_gyr << "gyr_x_" << ticket << ".dat";
      ofstream gyrFile(name_gyr.str().c_str(), ios::out);

      for(i = 0; i < gyration->bins; i++)
      {
          gyrFile << i / gyration->invDelta << " ";
          gyrFile << gyration->gyr[i] / gyration->counts[i] << " ";
          gyrFile << endl;
      }

      gyrFile.close();

  }






  /****************************************************************************/
  /*****************************  I/O ROUTINES  *******************************/
  /****************************************************************************/

/*
  void mocha_system::get_threads()
  {

    int nthreads, tid;

    #pragma omp parallel private(tid)
    {

      tid = omp_get_thread_num();
      if (tid == 0)
      {
	nthreads = omp_get_num_threads();
	cout << "Number of threads = " << nthreads << endl;
        Nthreads = nthreads;
      }

    }


  }
 * */

  void mocha_system::dump_info()
  {
    /** print basic info about the system on screen **/
    cout << "System has " << Natoms << " atoms\n";
    cout << "System has " << Nmolecules << " molecules\n";
    cout << "System has " << Nmol_types << " molecule types\n";
  }




  void mocha_system::get_data_folder(int in_argc, char *in_argv[])
  {

      if(in_argc == 1)
      {
          cout << "The system data are expected to be in the input folder\n";
          input_folder = "input";
      }
      else
      {
          cout << "The system data are expected to be in the " << in_argv[1] << " folder\n";
          input_folder = in_argv[1];
      }

  }


  void mocha_system::read_params()
  {

      stringstream filename;
      filename << input_folder << "/parameters.dat";

      ifstream inFile(filename.str().c_str(), ios::in);
      string line, string, string2, name;
      stringstream ss;
      double x;
      int i;
      bool boole, error;

      error = false;

      NB_CUTOFF = -1;
      NB_SKIN_HALF = -1;
      KTI_BEGIN = 1.0;
      KTI_END = 0.0;
      LAMBDA_VALUE = -1;
      INNER_SWEEPS = -1;
      OUTER_SWEEPS = -1;
      EQUILIBRATE = -1;
      BETA = -1;
      KIRKW_STEPS = -1;
      RMIN = -1;
      RMAX = -1;
      SIGMA_DISPL = -1;
      SIGMA_ROTAT = -1;
      SIGMA_BREATH = -1;
      SIGMA_COLLAPSE = 1.0;
      CGENERGY = false;
      AAENERGY = false;
      USEFEC = false;
      USEPOLYFEC = false;
      XDENS = -1;
      YDENS = -1;
      ZDENS = -1;
      RDENS = -1;
      GYR_BINS = -1;
      PRINT_TRAJ = 1;

      cout << "------------------------------------------------\n";
      while (inFile.good())
      {

          getline(inFile, line);

          ss << line;

          if(ss.str().length()>0)
          {
              if(ss.str().compare(0,1,"#")==0)
              {

                  string = "#NB_CUTOFF";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\tvalue " << x << endl;
                      NB_CUTOFF = x;
                      cutoff_nb = x;
                      INV_NB_CUTOFF = 1.0/x;
                      goto next;
                  }

                  string = "#NB_SKIN_HALF";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\tvalue " << x << endl;
                      NB_SKIN_HALF = x;
                      goto next;
                  }

                  string = "#INNER_SWEEPS";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      INNER_SWEEPS = i;
                      goto next;
                  }

                  string = "#OUTER_SWEEPS";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      OUTER_SWEEPS = i;
                      goto next;
                  }

                  string = "#EQUILIBRATE";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      EQUILIBRATE = i;
                      goto next;
                  }

                  string = "#BETA";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\t\tvalue " << x << endl;
                      BETA = x;
                      goto next;
                  }

                  string = "#KIRKW_STEPS";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      KIRKW_STEPS = i;
                      goto next;
                  }

                  string = "#RMIN";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\t\tvalue " << x << endl;
                      RMIN = x;
                      goto next;
                  }

                  string = "#RMAX";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\t\tvalue " << x << endl;
                      RMAX = x;
                      goto next;
                  }

                  string = "#SIGMA_DISPL";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\tvalue " << x << endl;
                      SIGMA_DISPL = x;
                      goto next;
                  }

                  string = "#SIGMA_ROTAT";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\tvalue " << x << endl;
                      SIGMA_ROTAT = x;
                      goto next;
                  }

                  string = "#SIGMA_BREATH";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\tvalue " << x << endl;
                      SIGMA_BREATH = x;
                      goto next;
                  }

                  string = "#SIGMA_COLLAPSE";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\tvalue " << x << endl;
                      SIGMA_COLLAPSE = x;
                      goto next;
                  }

                  string = "#HADRESS_MC";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      HADRESS_SGNL = boole;
                      goto next;
                  }

                  string = "#CONST_LAMBDA";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      CONST_LAMBDA = boole;
                      goto next;
                  }

                  string = "#LAMBDA_VALUE";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\tvalue " << x << endl;
                      LAMBDA_VALUE = x;
                      goto next;
                  }

                  string = "#HADRESS_SPHERE";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      SPHERE_SGNL = boole;
                      goto next;
                  }

                  string = "#HADRESS_SLAB";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> string2;
                      cout << "Parameter name " << name << "\tvalue " << string2 << endl;

                      if(string2.compare("x")==0)
                      {
                          XSLAB_SGNL = true;
                          HASLAB_SGNL = true;
                      }
                      if(string2.compare("y")==0)
                      {
                          YSLAB_SGNL = true;
                          HASLAB_SGNL = true;
                      }
                      if(string2.compare("z")==0)
                      {
                          ZSLAB_SGNL = true;
                          HASLAB_SGNL = true;
                      }

                      cout << "               Xslab: " << XSLAB_SGNL << " Yslab: " << YSLAB_SGNL << " Zslab: " << ZSLAB_SGNL << endl;
                      goto next;
                  }

                  string = "#KIRKWOOD_TI";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      KIRK_SGNL = boole;
                      goto next;
                  }

                  string = "#KTI_BEGIN";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\tvalue " << x << endl;
                      KTI_BEGIN = x;
                      goto next;
                  }

                  string = "#KTI_END";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> x;
                      cout << "Parameter name " << name << "\tvalue " << x << endl;
                      KTI_END = x;
                      goto next;
                  }

                  string = "#CGENERGY";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      CGENERGY = boole;
                      goto next;
                  }

                  string = "#AAENERGY";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      AAENERGY = boole;
                      goto next;
                  }

                  string = "#USE_FEC";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\t\tvalue " << boole << endl;
                      USEFEC = boole;
                      goto next;
                  }

                  string = "#USE_POLY_FEC";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      USEPOLYFEC = boole;
                      POLYFEC_SGNL = boole;
                      goto next;
                  }

                  string = "#DENSITY_X";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      XDENS = i;
                      goto next;
                  }

                  string = "#GYR_BINS";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      GYR_BINS = i;
                      goto next;
                  }

                  string = "#DENSITY_Y";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      YDENS = i;
                      goto next;
                  }

                  string = "#DENSITY_Z";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      ZDENS = i;
                      goto next;
                  }

                  string = "#DENSITY_R";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      RDENS = i;
                      goto next;
                  }

                  string = "#DO_PRESSURE";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      PRESSURE_SGNL = boole;
                      goto next;
                  }

                  string = "#DO_GYRATION";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      GYR_SGNL = boole;
                      goto next;
                  }

                  string = "#PRINT_TRJ";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      PRINT_TRAJ = i;
                      goto next;
                  }

                  string = "#DO_ELECTRO";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      ELECTRO_SGNL = boole;
                      goto next;
                  }

                  string = "#RIGID";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      RIGID = boole;
                      goto next;
                  }

                  string = "#ATOMIC_LAMBDA";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\tvalue " << boole << endl;
                      ATOMLAM = boole;
                      goto next;
                  }

                  string = "#PATH_INTEGRAL";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> boole;
                      cout << "Parameter name " << name << "\t\tvalue " << boole << endl;
                      PATHINT = boole;
                      goto next;
                  }

                  string = "#NESTING";
                  if(line.compare(0, string.length(), string)==0)
                  {
                      ss >> name >> i;
                      cout << "Parameter name " << name << "\tvalue " << i << endl;
                      NESTED = i;

/*
#pragma omp parallel
                      {
                          if (omp_get_thread_num() == 0)
                          {
                              omp_set_nested(NESTED);
                              if(omp_get_nested()==0)
                                  cout << "Nested parallelization **disabled**\n";
                              else
                                  cout << "Nested parallelization **enabled**\n";
                          }
                      }
 * */
                      goto next;
                  }

              }
          }

          next:
          ss.str("");
          ss.clear();


      }
      cout << "------------------------------------------------\n\n\n";

      inFile.close();

      cout << "Error checklisting...\n";

      if( (KIRK_SGNL == true) )
      {
          if( (PATHINT == false) && ( (CGENERGY == false) || (AAENERGY == false) ) )
          {
              cout << "Kirkwood integration was set up but AA or CG energy is off! Exiting...\n";
              error = true;
          }

          if(USEFEC == true)
          {
              cout << "Usage of free energy compensation is not compatible with Kirkwood integration! Exiting...\n";
              error = true;
          }

      }

      if( (CGENERGY == false) && (AAENERGY == false) )
      {
          cout << "Both AA and CG energies are off! Exiting...\n";
          error = true;
      }

      if( (HADRESS_SGNL == true) && (CONST_LAMBDA == true) )
      {
          cout << "H-AdResS simulation and constant-lambda simulation cannot coexist! Exiting...\n";
          error = true;
      }

      if( (KIRK_SGNL == true) && (CONST_LAMBDA == true) )
      {
          cout << "Kirkwood TI and constant-lambda simulation cannot coexist! Exiting...\n";
          error = true;
      }

      if(KIRK_SGNL == true)
      {
          if( KTI_BEGIN <= KTI_END )
          {
              cout << "KTI_BEGIN smaller than or equal KTI_END! Exiting...\n";
              error = true;
          }
      }

      if(KIRK_SGNL == true)
      {
          if( (KTI_BEGIN < 0) || (KTI_BEGIN > 1) || (KTI_END < 0) || (KTI_END > 1) )
          {
              cout << "KTI_BEGIN or/and KTI_END out of [0,1] interval! Exiting...\n";
              error = true;
          }
      }

      if( ((LAMBDA_VALUE < 0) || (LAMBDA_VALUE > 1)) && (CONST_LAMBDA == true) )
      {
          cout << "Need to specify lambda value for constant lambda simulation! Exiting...\n";
          error = true;
      }

      if( (HADRESS_SGNL == true) && (KIRK_SGNL == true) )
      {
          cout << "H-AdResS simulation and Kirkwood TI cannot coexist! Exiting...\n";
          error = true;
      }

      /*if( (HADRESS_SGNL == false) && (SPHERE_SGNL == true) )
      {
          cout << "H-AdResS simulation flag is off but sphere geometry is on! Exiting...\n";
          error = true;
      }

      if( (HADRESS_SGNL == false) && (HASLAB_SGNL == true) )
      {
          cout << "H-AdResS simulation flag is off but slab geometry is on! Exiting...\n";
          error = true;
      }*/

      if( (HASLAB_SGNL == true) && (SPHERE_SGNL == true) )
      {
          cout << "Shere and slab H-AdResS cannot coexist! Exiting...\n";
          error = true;
      }

      if( (HADRESS_SGNL == false) && (USEFEC == true) )
      {
          cout << "Usage of free energy compensation is not compatible with non-HAdResS simulations! Exiting...\n";
          error = true;
      }

      if( (HADRESS_SGNL == false) && (USEPOLYFEC == true) )
      {
          cout << "Usage of free energy compensation is not compatible with non-HAdResS simulations! Exiting...\n";
          error = true;
      }

      if( (USEFEC == true) && (USEPOLYFEC == true) )
      {
          cout << "Choose one type of FEC only! Exiting...\n";
          error = true;
      }

      if(HADRESS_SGNL == true)
      {
          if( (HASLAB_SGNL == false) && (SPHERE_SGNL == false) )
          {
              cout << "H-AdResS simulation flag is on but both sphere and slab geometry flags are off! Exiting...\n";
              error = true;
          }
      }

      if(NB_CUTOFF < 0)
      {
          cout << "Neighbor list cutoff not specified! Exiting...\n";
          error = true;
      }

      if(NB_SKIN_HALF < 0)
      {
          cout << "Skin not specified! Exiting...\n";
          error = true;
      }

      if(INNER_SWEEPS < 0)
      {
          cout << "Number of inner loops not specified! Exiting...\n";
          error = true;
      }

      if(OUTER_SWEEPS < 0)
      {
          cout << "Number of outer loops not specified! Exiting...\n";
          error = true;
      }

      if(BETA < 0)
      {
          cout << "Inverse temperature not specified! Exiting...\n";
          error = true;
      }

      if(EQUILIBRATE < 0)
      {
          cout << "Number of equilibration sweeps not specified! Exiting...\n";
          error = true;
      }

      if(KIRKW_STEPS < 0)
      {
          cout << "Number of Kirkwood TI steps not specified! Exiting...\n";
          error = true;
      }

      if(RMIN < 0)
      {
          cout << "RMIN not specified! Exiting...\n";
          error = true;
      }

      if(RMAX < 0)
      {
          cout << "RMAX not specified! Exiting...\n";
          error = true;
      }

      if(SIGMA_DISPL < 0)
      {
          cout << "SIGMA_DISPL not specified! Exiting...\n";
          error = true;
      }

      if(SIGMA_ROTAT < 0)
      {
          cout << "SIGMA_ROTAT not specified! Exiting...\n";
          error = true;
      }

      if(SIGMA_BREATH < 0)
      {
          cout << "SIGMA_BREATH not specified! Exiting...\n";
          error = true;
      }

      if( (SIGMA_COLLAPSE < 0) || (SIGMA_COLLAPSE > 1) )
      {
          cout << "SIGMA_COLLAPSE not in [0,1] interval! Exiting...\n";
          error = true;
      }

      if(RMAX <= RMIN)
      {
          cout << "RMAX is smaller then or equal to RMIN. Exiting...\n";
          error = true;
      }

      if( (NESTED != 0) && (NESTED != 1) )
      {
          cout << "The NESTING paramtere can be 0 or 1. Exiting...\n";
          error = true;
      }

      if(error == true)
      {
          exit(1);
      }
      else
      {
          cout << "None found. Proceeding.\n\n";
      }

      if (KIRK_SGNL == true)
      {
          cout << "Performing Kirkwood thermodynamic integration.\n\n";
      }
      else
      {
          if (HADRESS_SGNL == false)
          {
              cout << "Performing all-atom simulation.\n\n";
          }
          else
          {
              cout << "Performing H-AdResS simulation ";
              if(SPHERE_SGNL == true)
              {
                  cout << "with spherical geometry.\n\n";
              }
              else
              {
                  if(XSLAB_SGNL == true) cout << "with slab geometry in X direction.\n\n";
                  if(YSLAB_SGNL == true) cout << "with slab geometry in Y direction.\n\n";
                  if(ZSLAB_SGNL == true) cout << "with slab geometry in Z direction.\n\n";
              }
          }
      }



      /** assign the values of global lambdas:
       * they are always unity unless Kirkwood TI is performed **/
      lambda_AA = 1.0;
      lambda_CG = 1.0;
      if(KIRK_SGNL == true)
      {
          lambda_AA = KTI_BEGIN;
          lambda_CG = 1.0 - KTI_BEGIN;
      }
      if(CONST_LAMBDA == true)
      {
          lambda_AA = LAMBDA_VALUE;
          lambda_CG = 1.0 - LAMBDA_VALUE;
      }



  }




  void mocha_system::append_GRO_conf_unfolded(const char * filename, int t)
  {

    int i, j, k, Nat;
    vec coord;

    ofstream outFile(filename, ios::app);

    /*
    residue number (5 positions, integer)
    residue name (5 characters)
    atom name (5 characters)
    atom number (5 positions, integer)
    position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
     */

    outFile << "MoCHA simulation, t= " << t << endl;
    outFile << Natoms << endl;

    k = 0;
    for(i = 0; i < Nmolecules; i++)
    {

      Nat = molecule[i].Natoms;

      for(j = 0; j < Nat; j++)
      {

	coord = molecule[i].atom[j].position_unfolded;

	outFile << setw(5) << i\
                << setw(5) << molecule[i].type\
                << setw(5) << molecule[i].atom[j].type\
                << setw(5) << k\
                << setw(8) << setprecision(3) << fixed << coord.x()\
                << setw(8) << setprecision(3) << fixed << coord.y()\
                << setw(8) << setprecision(3) << fixed << coord.z()\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0 << endl;

        k++;

      }


    }

    outFile << box.side.x() << " " << box.side.y() << " " << box.side.z() << endl;

    outFile.close();

  }






  void mocha_system::append_XYZ_conf_unfolded(const char * filename, int t)
  {

    int i, j, k, Nat;
    vec coord;

    ofstream outFile(filename, ios::app);

    /*
     * Go for a super plain format for the local rdf evaluation tool.
     */

    //outFile << "MoCHA simulation, t= " << t << endl;
    //outFile << Natoms << endl;

    k = 0;
    for(i = 0; i < Nmolecules; i++)
    {

      Nat = molecule[i].Natoms;

      for(j = 0; j < Nat; j++)
      {

	coord = molecule[i].atom[j].position_unfolded;

	/*outFile << setw(5) << i\
                << setw(5) << k\
                << setw(8) << setprecision(3) << fixed << coord.x()\
                << setw(8) << setprecision(3) << fixed << coord.y()\
                << setw(8) << setprecision(3) << fixed << coord.z() << endl;*/

        outFile << i << " "\
                << j << " "\
                << setprecision(5) << fixed << coord.x() << " "\
                << setprecision(5) << fixed << coord.y() << " "\
                << setprecision(5) << fixed << coord.z() << endl;

        k++;

      }


    }

    //outFile << box.side.x() << " " << box.side.y() << " " << box.side.z() << endl;

    outFile.close();

  }






  void mocha_system::dump_GRO_conf_unfolded(const char * filename, int t)
  {

    int i, j, k, Nat;
    vec coord;

    ofstream outFile(filename, ios::out);

    /*
    residue number (5 positions, integer)
    residue name (5 characters)
    atom name (5 characters)
    atom number (5 positions, integer)
    position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
     */

    outFile << "MoCHA simulation, t= " << t << endl;
    outFile << Natoms << endl;

    k = 0;
    for(i = 0; i < Nmolecules; i++)
    {

      Nat = molecule[i].Natoms;

      for(j = 0; j < Nat; j++)
      {

	coord = molecule[i].atom[j].position_unfolded;

	outFile << setw(5) << i\
                << setw(5) << molecule[i].type\
                << setw(5) << molecule[i].atom[j].type\
                << setw(5) << k\
                << setw(8) << setprecision(3) << fixed << coord.x()\
                << setw(8) << setprecision(3) << fixed << coord.y()\
                << setw(8) << setprecision(3) << fixed << coord.z()\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0 << endl;

        k++;

      }


    }

    outFile << box.side.x() << " " << box.side.y() << " " << box.side.z() << endl;

    outFile.close();

    }






  void mocha_system::dump_GRO_conf(const char * filename, int t)
  {

    int i, j, k, Nat;
    vec coord;

    ofstream outFile(filename, ios::out);

    /*
    residue number (5 positions, integer)
    residue name (5 characters)
    atom name (5 characters)
    atom number (5 positions, integer)
    position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
     */

    outFile << "MoCHA simulation, t= " << t << endl;
    outFile << Natoms << endl;

    k = 0;
    for(i = 0; i < Nmolecules; i++)
    {

      Nat = molecule[i].Natoms;

      for(j = 0; j < Nat; j++)
      {

	coord = molecule[i].atom[j].position;

	outFile << setw(5) << i\
                << setw(5) << molecule[i].type\
                << setw(5) << molecule[i].atom[j].type\
                << setw(5) << k\
                << setw(8) << setprecision(3) << fixed << coord.x()\
                << setw(8) << setprecision(3) << fixed << coord.y()\
                << setw(8) << setprecision(3) << fixed << coord.z()\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0 << endl;

        k++;

      }


    }

    outFile << box.side.x() << " " << box.side.y() << " " << box.side.z() << endl;

    outFile.close();

    }




  void mocha_system::dump_GRO_CG_conf(const char * filename, int t)
  {

    int i;
    vec coord;

    ofstream outFile(filename, ios::out);

    /*
    residue number (5 positions, integer)
    residue name (5 characters)
    atom name (5 characters)
    atom number (5 positions, integer)
    position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
     */

    outFile << "MoCHA simulation, t= " << t << endl;
    outFile << Nmolecules << endl;

    for(i = 0; i < Nmolecules; i++)
    {

	coord = molecule[i].com;

	outFile << setw(5) << i\
                << setw(5) << molecule[i].type\
                << setw(5) << molecule[i].type\
                << setw(5) << i\
                << setw(8) << setprecision(3) << fixed << coord.x()\
                << setw(8) << setprecision(3) << fixed << coord.y()\
                << setw(8) << setprecision(3) << fixed << coord.z()\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0 << endl;

    }

    outFile << box.side.x() << " " << box.side.y() << " " << box.side.z() << endl;

    outFile.close();

    }




    void mocha_system::dump_MOC_conf(const char * filename, int t)
    {

    int i, j, k, Nat;
    vec coord;

    ofstream outFile(filename, ios::out);

    /*
    residue number (5 positions, integer)
    residue name (5 characters)
    atom name (5 characters)
    atom number (5 positions, integer)
    position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
     */

    outFile << "MoCHA simulation, t= " << t << endl;
    outFile << Natoms << endl;

    k = 1;
    for(i = 0; i < Nmolecules; i++)
    {

      Nat = molecule[i].Natoms;

      for(j = 0; j < Nat; j++)
      {

	coord = molecule[i].atom[j].position;

	outFile << setw(5) << i\
                << setw(5) << molecule[i].type\
                << setw(5) << molecule[i].atom[j].type\
                << setw(5) << k\
                << setw(10) << setprecision(5) << fixed << coord.x()\
                << setw(10) << setprecision(5) << fixed << coord.y()\
                << setw(10) << setprecision(5) << fixed << coord.z()\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0\
                << setw(8) << setprecision(3) << fixed << 0.0 << endl;

        //k++;

      }


    }

    outFile << box.side.x() << " " << box.side.y() << " " << box.side.z() << endl;

    outFile.close();

    }






void mocha_system::test_func()
{
    int SIZE = 4000;
    int nthreads, tid;
    int i, j, k;
    double c = 0.0;
    double e = 0.0;

    class aa_class
    {
    public:
        double x;
    };
    aa_class *aa = new aa_class [SIZE];

    /*
#pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            cout << "Number of threads = " << nthreads << endl;
        }
    }
     */

    for(i = 0; i < SIZE; i++) aa[i].x = i;

    void start_clock();

    static struct tms st_cpu;
    clock_t start = times(&st_cpu);
    double f = 0.0;

/*
#pragma omp parallel
    {

#pragma omp for private(i,j,k,c,f) schedule(dynamic) reduction(+:e)
        for(i = 0; i < SIZE; i++)
        {
            c = 0.0;
            f = 0.0;
            j = i;

            for(k = 0; k < 1000000; k++) c += sqrt(aa[j].x);
            for(k = 0; k < 1000000; k++) f += sqrt(0.9);
            e += c;
        }
    }
 * */

    double diff = (times(&st_cpu) - start) / (double) CLOCKS_PER_SEC;
    cout << "time " << diff << endl;

    if(nthreads == 1)
        cout << "Serial e is " << setprecision(10) << e << endl;
    else
        cout << "Parallel e is " << setprecision(10) << e << endl;

}










}

#endif























