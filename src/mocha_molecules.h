#ifndef _MOLS
#define	_MOLS
#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>

#include "vec.h"
#include "mocha_types.h"
//#include "random.h"


namespace mocha_molecules_nmspc
{
  
  using namespace std;
  using namespace mocha_interaction_nmspc;
  using namespace mocha_types_nmspc;
  using namespace mocha_energy_nmspc;
  using namespace mocha_myrandom_nmspc;
  
  
  
  class mocha_molecule
  {
    
    public:
      int Natoms;
      int index;
      int Nneighbors;
      int *neighbors;
      int type_number;
      int interaction_number;
      int MovedAtom;
      int MoveCase;
      bool ishybrid;
      bool new_ishybrid;
      bool isdeepCG;
      bool new_isdeepCG;
      bool is_in_HR;
      bool new_is_in_HR;
      double int_energy;
      double ext_energy;
      double kirkw_diff;
      double unweighted_virial_diff;
      double unweighted_energy_diff;
      double virial;
      double mass;
      double lambda;
      double new_lambda;
      double sigma_displ;
      double sigma_rotat;
      double sigma_breath;
      double sigma_collapse;
      //double angle;
      vec com;
      vec new_com;
      vec Move;
      string type;
      mocha_atom *atom;
      mocha_atom *new_atom;
      mocha_internal_interaction internal_int;
      mocha_box box;
      
      void allocate_memory(int n);
      void displace();
      void accept();
      void calc_mass();
      void get_lambda(int new_coord);
      double gyration_radius();
      double average_atom_distance();
      vec update_com();
      vec update_new_com();
      double lambda_function(double r);
      vec rotate(vec axis, vec position, double psi);
      double randrand_gaussian(double w)
      {
          return w*(rand() / double(RAND_MAX)-0.5);
      }
      
      
  };
  
  
  
  
  
  /** **************************************************** **/  
  
  
  
  void mocha_molecule::allocate_memory(int n)
  {

    int i;
    index = 0;
    com = vec(0.0, 0.0, 0.0);
    new_com = vec(0.0, 0.0, 0.0);
    box.initialize();
    atom = new mocha_atom [Natoms];
    new_atom = new mocha_atom [Natoms];
    Nneighbors = n;
    neighbors = new int [n];
    type = "";
    type_number = 0.0;
    int_energy = 0.0;
    ext_energy = 0.0;
    kirkw_diff = 0.0;
    unweighted_virial_diff = 0.0;
    unweighted_energy_diff = 0.0;
    virial = 0.0;
    mass = 0.0;
    lambda = 1.0;
    new_lambda = 1.0;
    ishybrid = true;
    new_ishybrid = true;
    isdeepCG = true;
    new_isdeepCG = true;
    is_in_HR = true;
    new_is_in_HR = true;
    Move = vec(0.0,0.0,0.0);
    MovedAtom = 0;
    MoveCase = 0;
    //angle = 0.0;
    for(i = 0; i < Natoms; i++) atom[i].initialize();
    
  }
  
  
  void mocha_molecule::get_lambda(int new_coord)
  {

        if (HADRESS_SGNL == true)
        {
            /* use one of the H-AdResS geometries */
            /* define the distance r according to the H-AdResS geometry */
            double r = -1.0; // the distance from the box center
            double lamMol = 1.0; // the calculated value of lambda
            int l = 0;
            
            if (SPHERE_SGNL == true)
            {
                if (new_coord == 0)
                    r = abs(com);
                else
                    r = abs(new_com);
            }

            if (XSLAB_SGNL == true)
            {
                if (new_coord == 0)
                    r = fabs(com.x());
                else
                    r = fabs(new_com.x());
            }

            if (YSLAB_SGNL == true)
            {
                if (new_coord == 0)
                    r = fabs(com.y());
                else
                    r = fabs(new_com.y());
            }

            if (ZSLAB_SGNL == true)
            {
                if (new_coord == 0)
                    r = fabs(com.z());
                else
                    r = fabs(new_com.z());
            }

            /* calculate lambda; assign the appropriate value to the boolean variable ishybrid */
            if (r >= box.rmax + cutoff_nb)
            {
                lamMol = 0.0;
                if (new_coord == 0)
                {
                    ishybrid = false;
                    is_in_HR = false;
                }
                else
                {
                    new_ishybrid = false;
                    new_is_in_HR = false;
                }
                
                // Set atomic lambdas
                if (ATOMLAM == true)
                {
                    if (new_coord == 0)
                    {
                        for(l = 0; l < Natoms; l++)
                        {
                                atom[l].lambda = 0.0;
                        }
                    }
                    else
                    {
                        for(l = 0; l < Natoms; l++)
                        {
                                new_atom[l].lambda = 0.0;
                        }
                    }
                }                                
            }
            
            // *****  Apparently not used  ***** //
            /*if (r >= box.rmax + cutoff_nb)
            {
                if (new_coord == 0)
                {
                    isdeepCG = true;
                }
                else
                {
                    new_isdeepCG = true;
                }
            }
            else
            {
                if (new_coord == 0)
                {
                    isdeepCG = false;
                }
                else
                {
                    new_isdeepCG = false;
                }
            }*/    

            if (r <= box.rmin - cutoff_nb)
            {
                lamMol = 1.0;
                if (new_coord == 0)
                {
                    ishybrid = true;
                    is_in_HR = false;
                }
                else
                {
                    new_ishybrid = true;
                    new_is_in_HR = false;
                }
                
                // Set atomic lambdas
                if (ATOMLAM == true)
                {
                    if (new_coord == 0)
                    {
                        for(l = 0; l < Natoms; l++)
                        {
                                atom[l].lambda = 1.0;
                        }
                    }
                    else
                    {
                        for(l = 0; l < Natoms; l++)
                        {
                                new_atom[l].lambda = 1.0;
                        }
                    }
                }                
            }
            if ((r > box.rmin - cutoff_nb) && (r < box.rmax + cutoff_nb))
            {
                if (new_coord == 0)
                {
                    ishybrid = true;
                    is_in_HR = true;
                }
                else
                {
                    new_ishybrid = true;
                    new_is_in_HR = true;
                }
                
                //lam = 1.0 - (r - box.rmin) * box.invhregsize;
                lamMol = lambda_function(r);
                
                // Set atomic lambdas
                if (ATOMLAM == true)
                {
                    if (new_coord == 0)
                    {
                        for(l = 0; l < Natoms; l++)
                        {
                                if (SPHERE_SGNL == true) {atom[l].lambda = lambda_function(abs(atom[l].position));}
                                if (XSLAB_SGNL == true) {atom[l].lambda = lambda_function(fabs(atom[l].position.x()));}
                                if (YSLAB_SGNL == true) {atom[l].lambda = lambda_function(fabs(atom[l].position.y()));}
                                if (ZSLAB_SGNL == true) {atom[l].lambda = lambda_function(fabs(atom[l].position.z()));}                        
                                //atom[l].lambda = lambda_function(atom[l].position);
                        }
                    }
                    else
                    {
                        for(l = 0; l < Natoms; l++)
                        {
                                if (SPHERE_SGNL == true) {new_atom[l].lambda = lambda_function(abs(new_atom[l].position));}
                                if (XSLAB_SGNL == true) {new_atom[l].lambda = lambda_function(fabs(new_atom[l].position.x()));}
                                if (YSLAB_SGNL == true) {new_atom[l].lambda = lambda_function(fabs(new_atom[l].position.y()));}
                                if (ZSLAB_SGNL == true) {new_atom[l].lambda = lambda_function(fabs(new_atom[l].position.z()));}    
                                //new_atom[l].lambda = lambda_function(new_atom[l].position);
                        }
                    }
                }                
            }

            if (new_coord == 0)
                lambda = lamMol;
            else
                new_lambda = lamMol;

        }
        else
        {
            /* performing a standard all-atom simulation */
            if (new_coord == 0)
                lambda = 1.0;
            else
                new_lambda = 1.0;
            
            ishybrid = true;
            new_ishybrid = true;
            isdeepCG = false;
            new_isdeepCG = false;
        }


  }




  void mocha_molecule::calc_mass()
  {
    
    int i;
    double m;
    
    m = 0.0;
    for(i = 0; i < Natoms; i++) m += atom[i].mass;
    
    mass = m;
    
  }
  
  
  
  
  
  
  double mocha_molecule::gyration_radius()
  {
    
    int i;
    double s, r;
    vec temp;
    
    r = 0.0;
    
    for(i = 0; i < Natoms; i++)
    {
      temp = pbcdist(atom[i].position, com, box);
      s = abs(temp);
      r += s * s;
    }
    
    return r / Natoms;
    
  }
  
  
  double mocha_molecule::average_atom_distance()
  {

    int i, j, n;
    double s, r;
    vec temp;

    r = 0.0;
    n = 0;

    for(i = 0; i < Natoms-1; i++)
    {
        for(j = i+1; j < Natoms; j++)
        {
            temp = pbcdist(atom[i].position, atom[j].position, box);
            s = abs(temp);
            r += s;
            n++;
        }
    }

    return r / n;

  }
  
  
  
  vec mocha_molecule::update_com()
  {
    
    int i, flag_x, flag_y, flag_z;
    vec temp_av, *temp;
    
    temp = new vec [Natoms];
    
    for(i = 0; i < Natoms; i++) temp[i] = atom[i].position;
    
    flag_x = 0;
    flag_y = 0;
    flag_z = 0;
    
    for(i = 1; i < Natoms; i++)
    {
      if( fabs(atom[0].position.x() - atom[i].position.x()) >= box.hside.x() ) flag_x = 1;
      if( fabs(atom[0].position.y() - atom[i].position.y()) >= box.hside.y() ) flag_y = 1;
      if( fabs(atom[0].position.z() - atom[i].position.z()) >= box.hside.z() ) flag_z = 1;
    }
    
    
    for(i = 1; i < Natoms; i++)
    {
      if((flag_x==1)&&(atom[0].position.x() - atom[i].position.x() > box.hside.x())) temp[i].x() += box.side.x();
      if((flag_x==1)&&(atom[0].position.x() - atom[i].position.x() < -box.hside.x())) temp[i].x() -= box.side.x();
      if((flag_y==1)&&(atom[0].position.y() - atom[i].position.y() > box.hside.y())) temp[i].y() += box.side.y();
      if((flag_y==1)&&(atom[0].position.y() - atom[i].position.y() < -box.hside.y())) temp[i].y() -= box.side.y();
      if((flag_z==1)&&(atom[0].position.z() - atom[i].position.z() > box.hside.z())) temp[i].z() += box.side.z();
      if((flag_z==1)&&(atom[0].position.z() - atom[i].position.z() < -box.hside.z())) temp[i].z() -= box.side.z();
    }
    
    temp_av = vec(0.0, 0.0, 0.0);
    for(i = 0; i < Natoms; i++)
    {
      temp_av += temp[i]*atom[i].mass;
      atom[i].position_unfolded = temp[i];
    }

    delete [] temp;
    
    com = temp_av/mass;
    com = pbcdist(com, vec(0.0, 0.0, 0.0), box);

    return com;
    
  }
  
  
  
  
  
  vec mocha_molecule::update_new_com()
  {

    int i, flag_x, flag_y, flag_z;
    vec temp_av, *temp;

    temp = new vec [Natoms];

    for(i = 0; i < Natoms; i++) temp[i] = new_atom[i].position;

    flag_x = 0;
    flag_y = 0;
    flag_z = 0;

    for(i = 1; i < Natoms; i++)
    {
      if( fabs(new_atom[0].position.x() - new_atom[i].position.x()) >= box.hside.x() ) flag_x = 1;
      if( fabs(new_atom[0].position.y() - new_atom[i].position.y()) >= box.hside.y() ) flag_y = 1;
      if( fabs(new_atom[0].position.z() - new_atom[i].position.z()) >= box.hside.z() ) flag_z = 1;
    }


    for(i = 1; i < Natoms; i++)
    {
      if((flag_x==1)&&(new_atom[0].position.x() - new_atom[i].position.x() > box.hside.x())) temp[i].x() += box.side.x();
      if((flag_x==1)&&(new_atom[0].position.x() - new_atom[i].position.x() < -box.hside.x())) temp[i].x() -= box.side.x();
      if((flag_y==1)&&(new_atom[0].position.y() - new_atom[i].position.y() > box.hside.y())) temp[i].y() += box.side.y();
      if((flag_y==1)&&(new_atom[0].position.y() - new_atom[i].position.y() < -box.hside.y())) temp[i].y() -= box.side.y();
      if((flag_z==1)&&(new_atom[0].position.z() - new_atom[i].position.z() > box.hside.z())) temp[i].z() += box.side.z();
      if((flag_z==1)&&(new_atom[0].position.z() - new_atom[i].position.z() < -box.hside.z())) temp[i].z() -= box.side.z();
    }

    temp_av = vec(0.0, 0.0, 0.0);
    for(i = 0; i < Natoms; i++)
    {
      temp_av += temp[i]*new_atom[i].mass;
      new_atom[i].position_unfolded = temp[i];
    }

    delete [] temp;

    new_com = temp_av/mass;
    new_com = pbcdist(new_com, vec(0.0, 0.0, 0.0), box);

    return new_com;

  }

  

  
  void mocha_molecule::displace()
  {

    int i, j, move, direction;
    double omega_norm, r, psi, scaler, shift, norm, randomv, decider;
    vec coord[Natoms], delta, omega, s, displ, int_def[Natoms], scaled_coord[Natoms];

    //move = rand() % (Natoms + 3);
    //move = rand() % (Natoms + 2);
    //move = rand() % 13; //15 13 11                                 REMOVE!
    move = dist13(mt);
    //move = SYSTEM
    if(move > 2) move = 2;                                 //REMOVE!
    //move = 2;          // (JUST FOR FULL CLASSICAL RUN!!)
    //move = 1;
    //move = rand() % 4;                                       //REMOVE!                                             
    //if(move > 0) {move = 3;} else {move = 0;}    
    //if(RIGID == true) move = 1 + rand() % 2;
    
    
    switch(move)
    {        
                /* random individual move of random atom */
                /**************************************************************/
            case 2:
                
                /* random individual moves of all atoms */
                // !! TEST !! //
                /*for (i = 0; i < Natoms; i++) {
                    norm = 2.0;
                    do{
                    s.x() = 2.0*(rand() / double(RAND_MAX)) - 1.0;
                    s.y() = 2.0*(rand() / double(RAND_MAX)) - 1.0;
                    s.z() = 2.0*(rand() / double(RAND_MAX)) - 1.0;
                    norm = abs(s);
                    }while(norm > 1.0);

                    randomv =  sigma_breath*gasdev();
                    s = randomv*s/norm;

                    new_atom[i].position = pbcdist(atom[i].position, -s, box);
                }
                MoveCase = 0;*/
                // !! TEST !! //
                
                /* random individual move of random atom */
                                                              //REMOVE!
                //j = rand() % Natoms;
                j = dist16(mt);
                //j = 0;
                /*s.x() = sigma_breath*gasdev();
                s.y() = sigma_breath*gasdev();
                s.z() = sigma_breath*gasdev();*/              //REMOVE!
                norm = 2.0;
                do{
                s.x() = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                s.y() = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                s.z() = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                norm = abs(s);
                }while(norm > 1.0);
                
                randomv =  sigma_breath*gasdev();
                //cout << "sigma_breath = " << sigma_breath << endl;
                //cout << "randomv = " << randomv << endl;
                s = randomv*s/norm;
                //cout << "s = " << s << endl;
                //cout << "norm = " << norm << endl << "Changing s to vec(0.0,0.0,0.0) ..." << endl;
                //s = vec(0.0,0.0,0.0);
                //cout << "s = " << s << endl << "NEXT STEP" << endl;
                for (i = 0; i < Natoms; i++) 
                {    
                    if ( i == j)
                    {
                    new_atom[i].position = pbcdist(atom[i].position, -s, box);
                    }
                    else
                    {
                    new_atom[i].position = atom[i].position;    
                    }
                }
                //cout << "pbcdist(atom[j].position, new_atom[j].position, box) = " << pbcdist(atom[j].position, new_atom[j].position, box) << endl;
                MovedAtom = j;
                Move = s;
                MoveCase = 0;
                                                              //REMOVE!
                
                
                
                
                
                /*
                for (i = 0; i < Natoms; i++) {
                    if (i == j) {
                    int_def[i] = s * atom[i].mass;
                    }
                    else {
                    int_def[i] = vec(0.0, 0.0, 0.0);    
                    }
                }
                delta = int_def[j] / mass;
                
                for (i = 0; i < Natoms; i++) {

                    displ = int_def[i] - delta;
                    new_atom[i].position = pbcdist(atom[i].position, -displ, box);
                }
                 */
                break;        
       
        
                /* breath */
                /**************************************************************//*
            case 0:
                delta = vec(0.0, 0.0, 0.0);
                for (i = 0; i < Natoms; i++) {
                    s.x() = sigma_breath*gasdev();
                    s.y() = sigma_breath*gasdev();
                    s.z() = sigma_breath*gasdev();
                    int_def[i] = s * atom[i].mass;
                    delta += int_def[i];
                }

                delta = delta / mass;

                for (i = 0; i < Natoms; i++) {

                    displ = int_def[i] - delta;
                    new_atom[i].position = pbcdist(atom[i].position, -displ, box);
                }
                break;*/


                /* displace CoM */
                /**************************************************************/
            case 1:
                norm = 2.0;
                do{
                s.x() = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                s.y() = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                s.z() = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                norm = abs(s);
                }while(norm > 1.0);
                
                randomv =  sigma_displ*gasdev();
                //cout << "randomv = " << randomv << endl;
                s = randomv*s/norm;
                //cout << "s = " << s << endl;
                //cout << "norm = " << norm << endl;
                
                for (i = 0; i < Natoms; i++) new_atom[i].position = pbcdist(atom[i].position, -s, box);
                //cout << "pbcdist(atom[i].position, new_atom[i].position, box) = " << pbcdist(atom[i].position, new_atom[i].position, box) << endl;}
                Move = s;
                MoveCase = 1;
                
                /* Only 3 orthogonal directions implementation */
                //direction = rand() % 3;
                
                //switch (direction) {
                        /* move along X */
                        /**************************************************************/
                    /*case 0:
                        s = vec(0.0, 0.0, 0.0);
                        s.x() = sigma_displ*gasdev();
                        for (i = 0; i < Natoms; i++) new_atom[i].position = pbcdist(atom[i].position, -s, box);
                        Move = s;
                        MoveCase = 1;
                        break;

                        /* move along Y */
                        /**************************************************************/
                    /*case 1:
                        s = vec(0.0, 0.0, 0.0);
                        s.y() = sigma_displ*gasdev();
                        for (i = 0; i < Natoms; i++) new_atom[i].position = pbcdist(atom[i].position, -s, box);
                        Move = s;
                        MoveCase = 1;
                        break;

                        /* move along Z */
                        /**************************************************************/
                    /*case 2:
                        s = vec(0.0, 0.0, 0.0);
                        s.z() = sigma_displ*gasdev();
                        for (i = 0; i < Natoms; i++) new_atom[i].position = pbcdist(atom[i].position, -s, box);
                        Move = s;
                        MoveCase = 1;
                        break;

                }*/
                
                break;


                /* rotate */
                /**************************************************************/
            case 0:
                // coordinates in the COM reference frame
                for (i = 0; i < Natoms; i++) coord[i] = pbcdist(atom[i].position, com, box);

                // rotation vector
                r = 2.0;
                do {
                    omega.x() = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                    omega.y() = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                    omega.z() = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                    r = abs(omega);
                } while (r > 1.0);

                // normalize the rotation vector and assign the right angle
                omega_norm = abs(omega);
                omega = omega / omega_norm;

                // get rotation angle
                psi = sigma_rotat*gasdev();

                // rotate molecule
                for (i = 0; i < Natoms; i++) {
                    displ = rotate(omega, coord[i], psi);
                    new_atom[i].position = pbcdist(displ, -com, box); 
                    //cout << "pbcdist(atom[i].position, new_atom[i].position, box) = " << pbcdist(atom[i].position, new_atom[i].position, box) << endl;
                }
                //Move = omega;
                //angle = psi;        
                MoveCase = 2;
                break;
                
                
                /* ring collapse/expansion or breath move */
                /**************************************************************/
            case 3:
                
                // #####   The Scaling Approach   ##### //
                
                // Value by which the bead-CoM vector is scaled
                
                //Equally distributed random number between 1.0-sigma_collapse and 1.0.
                //scaler = sigma_collapse*(rand() / double(RAND_MAX)) + 1.0-sigma_collapse;
                scaler = sigma_collapse*ddist2(mt) + 1.0-sigma_collapse;
                /*//Gaussian distributed random number between 0.0 and 1.0 with peak on 1.0 and width sigma_collapse.
                //scaler = 1.0-fabs(sigma_collapse*gasdev());                
                do{
                    scaler = 1.0-fabs(sigma_collapse*gasdev());
                    if( scaler > 1.0 )
                    {
                        cerr << "Error in CoM-bead scaler: scaler > 1.0. Scaler = " << scaler << endl;
                        exit(1);
                    }
                } while (scaler < 0.0);*/
                
                if( (scaler > 1.0) || (scaler < 1.0 -sigma_collapse) )
                {
                    cerr << "Error in CoM-bead scaler: scaler > 1.0 or < 1.0 - sigma_collapse, scaler = " << scaler << endl;
                    exit(1);
                }
                
                // Decides whether we scale with x or 1/x (To keep detailed balance we need a probability distribution with equal probability for x and 1/x)
                //decider = rand()%2;
                decider =ddist2(mt);// rand() / double(RAND_MAX);
                
                // Scale bead-CoM vectors with scaler or 1/scaler according to decider
                if(decider >= 0.5)
                {    
                    for (i = 0; i < Natoms; i++)
                    { 
                        scaled_coord[i] = pbcdist(atom[i].position, com, box)*scaler;
                    }
                    //cout  << "Scaler " << scaler << "\n";
                }
                else //if(decider < 0.5)
                {
                    for (i = 0; i < Natoms; i++)
                    { 
                        scaled_coord[i] = pbcdist(atom[i].position, com, box)/scaler;
                    }
                    //cout << "Scaler " << 1/scaler << "\n";
                }
                /*else
                {
                    cerr << "Error in decider: neither 1, nor 0, decider = " << decider << endl;
                    exit(1);
                }*/
                
                // Put back the CoM vector and move atoms
                for (i = 0; i < Natoms; i++) {
                    new_atom[i].position = pbcdist(scaled_coord[i], -com, box);
                    //cout << "pbcdist(atom[i].position, new_atom[i].position, box) = " << pbcdist(atom[i].position, new_atom[i].position, box) << endl;
                }
                
                
                // #####   The Scaling Approach With Partially Flat Probability Distribution   ##### //
                
                /* // Decides whether we scale with x or 1/x (To keep detailed balance we need a probability distribution with equal probability for x and 1/x)
                decider = rand() / double(RAND_MAX);
                
                // Scale bead-CoM vectors with scaler or 1/scaler according to decider
                if(decider >= 0.5)
                {                  
                    // Value by which the bead-CoM vector is scaled   
                    scaler = sigma_collapse*(rand() / double(RAND_MAX)) + 1.0-sigma_collapse;
                 
                    for (i = 0; i < Natoms; i++)
                    { 
                        scaled_coord[i] = pbcdist(atom[i].position, com, box)*scaler;
                    }
                }
                else
                {
                    // Value by which the bead-CoM vector is scaled
                    scaler = (sigma_collapse/(1.0-sigma_collapse))*(rand() / double(RAND_MAX)) + 1.0;
                 
                    for (i = 0; i < Natoms; i++)
                    { 
                        scaled_coord[i] = pbcdist(atom[i].position, com, box)*scaler;
                    }
                }
                
                // Put back the CoM vector and move atoms
                for (i = 0; i < Natoms; i++) {
                    new_atom[i].position = pbcdist(scaled_coord[i], -com, box);
                }*/
                
                
                // #####   The Scaling Approach With Fully Flat Probability Distribution   ##### //
                
                // Value by which the bead-CoM vector is scaled                
                // Fully flat distribution between 1.0-sigma_collapse and 1.0/(1.0-sigma_collapse) 
                /*scaler = (1.0/(1.0-sigma_collapse)-(1.0-sigma_collapse)) * (rand() / double(RAND_MAX)) + 1.0-sigma_collapse;
                
                for (i = 0; i < Natoms; i++)
                {
                    // Scale coordinates with scaler
                    scaled_coord[i] = pbcdist(atom[i].position, com, box)*scaler;
                    
                    // Put back the CoM vector and move atoms
                    new_atom[i].position = pbcdist(scaled_coord[i], -com, box);
                }*/                

                               
                // #####   The Shifting Approach   ##### //
                
                // Value by which the bead-CoM vector is shifted
                
                //Gaussian distributed random number
                //shift = sigma_collapse*gasdev();
                
                //Equally distributed random number between -sigma_collapse and +sigma_collapse.
                /*shift = sigma_collapse * (2.0*(rand() / double(RAND_MAX)) - 1.0);
                
                if( (shift > sigma_collapse) || (shift < -sigma_collapse) )
                {
                    cerr << "Error in CoM-bead shift: shift > sigma_collapse or < -sigma_collapse, shift = " << shift << endl;
                    exit(1);
                }
                
                // Go to CoM reference frame, calculate unit vector, add scaled unit vector, 
                for (i = 0; i < Natoms; i++) {
                    coord[i] = pbcdist(atom[i].position, com, box);
                    scaled_coord[i] = coord[i] + (coord[i]/abs(coord[i]))*shift;
                    new_atom[i].position = pbcdist(scaled_coord[i], -com, box);
                }*/
                
                
                MoveCase = 3;
                break;
                


        }
    
    displ = update_new_com();
    new_com = displ;
    get_lambda(1);

  }
  
  
  void mocha_molecule::accept()
  {
    
    int i;
    
    for(i = 0; i < Natoms; i++)
    {
        atom[i].position = new_atom[i].position;
        atom[i].lambda = new_atom[i].lambda;
    }
    com = new_com;
    lambda = new_lambda;
    ishybrid = new_ishybrid;
    is_in_HR = new_is_in_HR;
    isdeepCG = new_isdeepCG;
    
  }
  
  
  double mocha_molecule::lambda_function(double r)
  {
      //cout << "Lambda function called!\n";
      double lam;

      if (r >= box.rmax)
      {
          lam = 0.0;
      }
      if (r <= box.rmin)
      {
          lam = 1.0;
      }
      
      if((r > box.rmin) && (r < box.rmax))
      {
          //lam = 1.0 - (r - box.rmin) * box.invhregsize;
          //lam = ( 1.0 - ((r - box.rmin) * box.invhregsize) * ((r - box.rmin) * box.invhregsize) )*( 1.0 - ((r - box.rmin) * box.invhregsize) * ((r - box.rmin) * box.invhregsize) );
          lam = pow(cos(box.invhregsize * M_PI_2 * (r - box.rmin)),2.0);
          //lam = (1.0/99.0) * ( 100.0 - ( 1.0/( pow(pow(cos(box.invhregsize * M_PI_2 * (r - box.rmin)),2.0) + (1.0/10.0)*pow(sin(box.invhregsize * M_PI_2 * (r - box.rmin)),2.0)  , 2.0 ) ) ) );
          //lam = (1.0 / 99.0) * (100.0 - (1.0 / (pow(pow(cos(box.invhregsize * M_PI_2 * (r - box.rmin)), 2.0) + (0.6309573444801932) * pow(sin(box.invhregsize * M_PI_2 * (r - box.rmin)), 2.0), 10.0))));
      }

      return lam;

  }

  
  
  
  
  vec mocha_molecule::rotate(vec axis, vec position, double psi)
  {

        double q0, q1, q2, q3;
        vec rot_0, rot_1, rot_2, result;

        /* use quaternions */
        q0 = cos(psi / 2.0);
        q1 = axis.getX() * sin(psi / 2);
        q2 = axis.getY() * sin(psi / 2);
        q3 = axis.getZ() * sin(psi / 2);


        rot_0.x() = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
        rot_0.y() = 2 * (q1 * q2 - q0 * q3);
        rot_0.z() = 2 * (q0 * q2 + q1 * q3);

        rot_1.x() = 2 * (q0 * q3 + q1 * q2);
        rot_1.y() = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
        rot_1.z() = 2 * (q2 * q3 - q0 * q1);

        rot_2.x() = 2 * (q1 * q3 - q0 * q2);
        rot_2.y() = 2 * (q0 * q1 + q2 * q3);
        rot_2.z() = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

        result.x() = rot_0 * position;
        result.y() = rot_1 * position;
        result.z() = rot_2 * position;

        return result;

    }

  
  
  
  
}

#endif


