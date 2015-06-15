
/********************************/
/* MoCHA - Monte Carlo H-AdResS */
/********************************/

bool KIRK_SGNL = false;
bool HADRESS_SGNL = false;
bool SPHERE_SGNL = false;
bool HASLAB_SGNL = false;
bool XSLAB_SGNL = false;
bool YSLAB_SGNL = false;
bool ZSLAB_SGNL = false;
bool PRESSURE_SGNL = false;
bool GYR_SGNL = false;
bool ELECTRO_SGNL = false;
bool POLYFEC_SGNL = false;
bool RIGID = false;
bool ATOMLAM = false;
bool PATHINT = false;
bool CONST_LAMBDA = false;
int NESTED = 0;
double DIELECTRIC = 138.9354426312346;
double PI = 3.14159265358979323846;
double inverse_frame = 1.0;
double cutoff_nb = 0.0;
//int Num_fec_poly_params = 6;
//double epsilon_param_rho = 10.0;
#define OUTPUTFREQ (1)


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iomanip>
//#include <omp.h>
#include <math.h>
#include "src/headers.h"
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <chrono>
#include <ctime>
using namespace std;
using namespace mocha_interaction_nmspc;
using namespace mocha_types_nmspc;
using namespace mocha_system_nmspc;
using namespace mocha_numbers_nmspc;
using namespace mocha_molecules_nmspc;
using namespace mocha_matrixInv_nmspc;

using namespace mocha_myrandom_nmspc;

int main(int argc, char *argv[])
{   
    cout << "Welcome in MoCHA\n";
    //srand(time(NULL));
    srand(getpid());
    mocha_system System;                                    /** Build the system **/
    //System.get_threads();                                   /** Acquire number of available threads **/
    System.get_data_folder(argc, argv);                     /** Read folder in which input data are located **/
    System.read_params();                                   /** read simulation parameters from file **/
    //System.read_GRO_input_file();                           /** get # atoms, # molecules, # molecule types - NOT COORDINATES **/
    System.read_MOC_input_file();                           /** get # atoms, # molecules, # molecule types - NOT COORDINATES **/
    System.dump_info();                                     /** print basic info on screen **/
    System.initialize_attempt_success();                        /** allocate memory for attemps/success arrays **/
    System.initialize_interaction();                        /** allocate as many interaction structures as molecule types are **/
    System.initialize_molecules();                          /** allocate molecule structures and initialize # of atoms in each molecule to zero **/
    //System.fill_GRO_molecules();                            /** Read the input conf file and put coordinates in the atom structures **/
    System.fill_MOC_molecules();                            /** Read the input conf file and put coordinates in the atom structures **/
    System.com_shift_center();                              /** Center the whole system in the coordinate origin **/
    System.read_internal_interaction_file();                /** Read the internal topology file **/
    System.read_bonded_intramol_interaction_file();         /** Read the intermolecular topology file **/
    System.read_lj_parameters();                            /** read lj interaction parameters **/
    System.get_molecule_masses();                           /** Calculate the mass of each molecule from the atoms' mass **/
    System.update_molecule_com();                           /** Calculate the molecules' CoM and fold it according to the PBC **/
    //System.update_neighbors();                            /** Build neighbor list **/
    System.update_neighbors_PI();                           /** Build neighbor list - PATH INTEGRALS **/
    System.UpdateFlag = false;
    System.VerletListUpdates = 0;
    System.EnergyTime = 0.0;
    System.build_new_atoms();                               /** Copy in the new_atom fields the dta of the atom structures **/
    System.number_of_molecules_per_type();                  /** Calculate number of molecules of each type and allocate it in a vector **/
    System.initialize_density_structure();                  /** Allocate and initialize the arrays containing the density **/
    System.initialize_gyration_structure();                  /** Allocate and initialize the gyration structure **/
    if(System.USEFEC==true) System.read_and_allocate_fec(); /** Make use of free energy compensation **/
    if(System.USEPOLYFEC==true) System.initialize_fec_poly_params(); /** Make use of free energy compensation **/
    System.total_energy();                                  /** Calculate total energy of the initial configuration **/



    int i, j, ticket;
    stringstream name_log;
    stringstream name_gyr;
    stringstream name_GROtraj;
    stringstream name_XYZtraj;
    stringstream name_out;
    stringstream name_unfolded_out;
    stringstream name_CGout;
    stringstream name_MOCout;
    stringstream name_particleQM;

    if(CONST_LAMBDA == true)
    {
        cout << "This is a const-lambda simulation: lambda_AA = " << System.lambda_AA << ". lambda_CG = " << System.lambda_CG << ".\n\n";
    }
    
    ticket = getpid();
    cout << "Ticket is " << ticket << endl;

    name_log << "log_" << ticket << ".dat";
    name_gyr << "gyr_" << ticket << ".dat";
    name_GROtraj << "traj_" << ticket << ".gro";
    name_XYZtraj << "traj_xyz_" << ticket << ".xyz";
    name_out << "output_" << ticket << ".gro";
    name_unfolded_out << "output_unfolded_" << ticket << ".gro";
    name_CGout << "CGoutput" << ticket << ".gro";
    name_MOCout << "MOCoutput" << ticket << ".moc";
    name_particleQM << "particleQM_" << ticket << ".dat";
    ofstream logFile(name_log.str().c_str(), ios::out);
    ofstream gyrFile(name_gyr.str().c_str(), ios::out);
    ofstream partFile(name_particleQM.str().c_str(), ios::out);
    

    /************************ Kirkwood TI variables ***************************/
    int k, n;
    double *kfed, *temp_K, *temp_V, *virial;
    double helmholtz, pressure;
    stringstream name_KTI, name_POLYFEC;
    name_KTI << "Kirkwood_free_energy_diff_" << ticket << ".dat";
    ofstream KirkwFile(name_KTI.str().c_str(), ios::out);
    name_POLYFEC << "FEC_params_" << ticket << ".dat";
    ofstream PolyFecFile(name_POLYFEC.str().c_str(), ios::out);
    kfed = new double[System.Nmol_types];
    temp_K = new double[System.Nmol_types];
    temp_V = new double[System.Nmol_types];
    virial = new double[System.Nmol_types];

    if( (KIRK_SGNL == true) || (PRESSURE_SGNL == true) )
    {
        for (k = 0; k < System.Nmol_types; k++)
        {
            kfed[k] = 0.0;
            virial[k] = 0.0;
            temp_K[k] = 0.0;
            temp_V[k] = 0.0;
        }
    }

    if(KIRK_SGNL == true)
    {
        n = 0;
        KirkwFile << "# L ";
        for (k = 0; k < System.Nmol_types; k++) KirkwFile <<\
                " <DV>_" << System.interaction.internal[k].type <<\
                " DP_" << System.interaction.internal[k].type;
        KirkwFile << endl;
    }
    /**************************************************************************/

    
    
  

    /*************************** PRINT INPUT **********************************/
    System.append_GRO_conf_unfolded(name_GROtraj.str().c_str(), 0);
    System.append_XYZ_conf_unfolded(name_XYZtraj.str().c_str(), 0);
    /************************* END PRINT INPUT ********************************/    
    
    

    /************************** EQUILIBRATION *********************************/ 
    for(i = 0; i < System.EQUILIBRATE; i++)
    {

        for(j = 0; j < System.INNER_SWEEPS; j++)
        {
                System.sweep();
                //cout << "Inner sweep done." << endl;
                /*if (System.UpdateFlag == true) 
                {
                    System.UpdateFlag = false;
                    //System.update_neighbors();
                    System.update_neighbors_PI();  // For PATH INTEGRAL SIMULATIONS
                    //cout << "Updated neighborlist." << endl;
                }*/
        }

        //System.update_neighbors();
        System.total_energy();
        
        if(GYR_SGNL == true)
        {
            System.mean_gyration();
            gyrFile << "RADIUSOFGYRATION " << i + 1 << " " << setprecision(9) << System.rad_of_gyr << endl;
        }
        
        System.calc_partQM();
        partFile << "PARTICLES_QM " << i + 1 << " " << System.partQM << endl;        
        logFile << "ENERGY_EQ " << i + 1 << " " << setprecision(9) << System.tot_energy/System.BETA << endl;
        if(i%OUTPUTFREQ==0) 
        {
            cout << "Equilibration sweep " << i << " done. Acceptance ratios: " << endl;
            cout << "Acceptance ratio Breath Moves: " << (double)System.SUCCESSES[0]/(double)System.ATTEMPTS[0] << " for " << System.ATTEMPTS[0] << " attempted moves." << endl;
            cout << "Acceptance ratio Center of Mass Moves: " << (double)System.SUCCESSES[1]/(double)System.ATTEMPTS[1] << " for " << System.ATTEMPTS[1] << " attempted moves." << endl;
            cout << "Acceptance ratio Rotation Moves: " << (double)System.SUCCESSES[2]/(double)System.ATTEMPTS[2] << " for " << System.ATTEMPTS[2] << " attempted moves." << endl;
            cout << "Acceptance ratio Collapse/Expand Moves: " << (double)System.SUCCESSES[3]/(double)System.ATTEMPTS[3] << " for " << System.ATTEMPTS[3] << " attempted moves." << endl;
        }

        if(i%System.PRINT_TRAJ==0)
        {
            System.append_GRO_conf_unfolded(name_GROtraj.str().c_str(), i + 1);
            System.append_XYZ_conf_unfolded(name_XYZtraj.str().c_str(), i + 1);
        }

    }
    /********************** END EQUILIBRATION *********************************/


    
    /************** ALL ATOM OR H-AdResS OR CONST-LAMBDA RUN ******************/
    if(KIRK_SGNL == false)
    {
        
        // Start the timer
        cout << "Starting Timer." << endl;
        auto t_start = std::chrono::high_resolution_clock::now();
        
        for (i = 0; i < System.OUTER_SWEEPS; i++)
        {            
            
            inverse_frame = 1.0/(i+1);
            
            for (j = 0; j < System.INNER_SWEEPS; j++)
            {
                System.sweep();
                //cout << "Inner Sweep Done." << endl;
                /*if (System.UpdateFlag == true)
                {
                    System.UpdateFlag = false;
                    //System.update_neighbors();
                    System.update_neighbors_PI();  // For PATH INTEGRAL SIMULATIONS
                    //cout << "Updated neighborlist." << endl;
                }*/                
                
                if(System.USEPOLYFEC==true) System.poly_fec_coefficient_update();
            }

            //System.update_neighbors();
            System.total_energy();
            
            if(GYR_SGNL == true)
            {
                System.mean_gyration();
                gyrFile << "RADIUSOFGYRATION " << i + 1 + System.EQUILIBRATE << " " << setprecision(9) << System.rad_of_gyr << endl;
            }

            if(System.XDENS > 0) System.calc_density(0);
            if(System.YDENS > 0) System.calc_density(1);
            if(System.ZDENS > 0) System.calc_density(2);
            if(System.RDENS > 0) System.calc_density(3);
            if(System.GYR_BINS > 0) System.calc_gyration();

            System.calc_partQM();
            
            if(PRESSURE_SGNL == false)
            {
                partFile << "PARTICLES_QM " << i + 1 + System.EQUILIBRATE << " " << System.partQM << endl;
                logFile << "ENERGY " << i + 1 + System.EQUILIBRATE << " " << setprecision(9) << System.tot_energy/System.BETA << endl;
            }
            else
            {
                partFile << "PARTICLES_QM " << i + 1 + System.EQUILIBRATE << " " << System.partQM << endl;
                logFile << "ENERGY " << i + 1 + System.EQUILIBRATE << " " << setprecision(9) << System.tot_energy/System.BETA;
                System.kirkwood_energy_difference(kfed, virial);
                for (k = 0; k < System.Nmol_types; k++)
                {
                    pressure = (System.molecules_per_type.atoms[k] * System.molecules_per_type.number[k] - virial[k] / 3.0 ) / ( System.BETA * System.box.volume);
                    logFile << " P_" << k << " " << pressure;
                }
                logFile << endl;
            }
            
            //if(i%OUTPUTFREQ==0) cout << "Run sweep " << i << " done\n";
            if(i%OUTPUTFREQ==0)
            {   
                cout << "Run sweep " << i << " done. Acceptance ratios: " << endl;
                cout << "Acceptance ratio Breath Moves: " << (double)System.SUCCESSES[0]/(double)System.ATTEMPTS[0] << " for " << System.ATTEMPTS[0] << " attempted moves." << endl;
                cout << "Acceptance ratio Center of Mass Moves: " << (double)System.SUCCESSES[1]/(double)System.ATTEMPTS[1] << " for " << System.ATTEMPTS[1] << " attempted moves." << endl;
                cout << "Acceptance ratio Rotation Moves: " << (double)System.SUCCESSES[2]/(double)System.ATTEMPTS[2] << " for " << System.ATTEMPTS[2] << " attempted moves." << endl;
                cout << "Acceptance ratio Collapse/Expand Moves: " << (double)System.SUCCESSES[3]/(double)System.ATTEMPTS[3] << " for " << System.ATTEMPTS[3] << " attempted moves." << endl;
            }
        
            if(i%(3*OUTPUTFREQ)==0)
            {
                if(System.USEPOLYFEC==true)
                {
                    System.poly_fec_coefficient_calc();
                    for(k = 0; k < System.Nmol_types; k++)
                    {
                        PolyFecFile << i << " ";
                        for(j = 0 ; j < System.Num_fec_poly_params; j++) PolyFecFile << System.fec_poly_params[k].a[j] + System.fec_poly_params_rho[k].a[j] << " " ;
                        PolyFecFile << "MOLNUM_" << k << endl;
                    }
                }
            }

            if(i%System.PRINT_TRAJ==0)
            {
                System.append_GRO_conf_unfolded(name_GROtraj.str().c_str(), i + System.EQUILIBRATE + 1);
                System.append_XYZ_conf_unfolded(name_XYZtraj.str().c_str(), i + System.EQUILIBRATE + 1);
            }

        }
        
        // End the timer
        auto t_end = std::chrono::high_resolution_clock::now();
        cout << "Timer stopped. Wall clock time passed: " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
        logFile << "Wall clock time passed: " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms\n";
        cout << "Number of Verlet List Updates: " << System.VerletListUpdates << endl;
        logFile << "Number of Verlet List Updates: " << System.VerletListUpdates << "\n";
        cout << "Time spend for WCA or SG energy calculations: " << System.EnergyTime << " ns" << endl;
        logFile << "Time spend for WCA or SG energy calculations: " << System.EnergyTime << " ns\n";
        
    }
    /****************** END ALL ATOM OR H-AdResS RUN **************************/






    /**************************** KIRKWOOD TI *********************************/
    /******  WARNING! The free energy is printed out multiplied by BETA  ******/
    /**************************************************************************/
    if (KIRK_SGNL == true)
    {


        for (i = 0; i <= System.OUTER_SWEEPS + System.KIRKW_STEPS; i++)
        {



            if ((i % System.KIRKW_STEPS) == 0)
            {

                if (i > 0)
                {                    
                    KirkwFile << System.lambda_AA;
                    for (k = 0; k < System.Nmol_types; k++)
                    {
                        helmholtz = kfed[k] / ( n * System.molecules_per_type.number[k]);
                        pressure = (System.molecules_per_type.atoms[k] - virial[k] / ( n * 3.0 * System.molecules_per_type.number[k] ));
                        KirkwFile << " " << helmholtz << " " << pressure;
                    }
                    KirkwFile << endl;
                }
                
                for (k = 0; k < System.Nmol_types; k++)
                {
                    kfed[k] = 0.0;
                    virial[k] = 0.0;
                }
                n = 0;


                System.lambda_AA = System.KTI_BEGIN - ((double) i / (double) (System.OUTER_SWEEPS)) * (System.KTI_BEGIN - System.KTI_END);
                System.lambda_CG = 1.0 - System.KTI_BEGIN + ((double) i / (double) (System.OUTER_SWEEPS)) * (System.KTI_BEGIN - System.KTI_END);
                
                //System.lambda_AA = 1.0 - (double) i / (double) (System.OUTER_SWEEPS);
                //System.lambda_CG = (double) i / (double) (System.OUTER_SWEEPS);
                
                cout << "Next KTI step, set new lambda's: lambda_AA = " << System.lambda_AA << " and lambda_CG = " << System.lambda_CG << ".\n";

            }


            for (j = 0; j < System.INNER_SWEEPS; j++)
            {

                System.sweep();
                //cout << "Inner Sweep Done." << endl;
                /*if (System.UpdateFlag == true)
                {
                    System.UpdateFlag = false;
                    //System.update_neighbors();
                    System.update_neighbors_PI();  // For PATH INTEGRAL SIMULATIONS
                    //cout << "Updated neighborlist." << endl;
                }*/
                
                System.kirkwood_energy_difference(temp_K, temp_V);

                for (k = 0; k < System.Nmol_types; k++)
                {
                    kfed[k] += temp_K[k];
                    virial[k] += temp_V[k];
                }                
                n++;
                
            }

            //System.update_neighbors();
            System.calc_partQM();
            System.total_energy();
            partFile << "PARTICLES_QM " << i + 1 + System.EQUILIBRATE << " " << System.partQM << endl;
            logFile << "ENERGY " << i + 1 + System.EQUILIBRATE << " " << setprecision(9) << System.tot_energy/System.BETA << endl;
            
            if(GYR_SGNL == true)
            {
                System.mean_gyration();
                gyrFile << "RADIUSOFGYRATION " << i + 1 + System.EQUILIBRATE << " " << setprecision(9) << System.rad_of_gyr << endl;
            }
            
            if(i%OUTPUTFREQ==0) 
            {   
                cout << "TI sweep " << i << " done. Acceptance ratios: " << endl;
                cout << "Acceptance ratio Breath Moves: " << (double)System.SUCCESSES[0]/(double)System.ATTEMPTS[0] << " for " << System.ATTEMPTS[0] << " attempted moves." << endl;
                cout << "Acceptance ratio Center of Mass Moves: " << (double)System.SUCCESSES[1]/(double)System.ATTEMPTS[1] << " for " << System.ATTEMPTS[1] << " attempted moves." << endl;
                cout << "Acceptance ratio Rotation Moves: " << (double)System.SUCCESSES[2]/(double)System.ATTEMPTS[2] << " for " << System.ATTEMPTS[2] << " attempted moves." << endl;
                cout << "Acceptance ratio Collapse/Expand Moves: " << (double)System.SUCCESSES[3]/(double)System.ATTEMPTS[3] << " for " << System.ATTEMPTS[3] << " attempted moves." << endl;
            }
            
            if(i%System.PRINT_TRAJ==0) 
            {
                System.append_GRO_conf_unfolded(name_GROtraj.str().c_str(), i + System.EQUILIBRATE + 1);
                System.append_XYZ_conf_unfolded(name_XYZtraj.str().c_str(), i + System.EQUILIBRATE + 1);
            }

        }


    }
    /************************ END KIRKWOOD TI *********************************/


    if( (System.XDENS > 0) || (System.YDENS > 0) || (System.ZDENS > 0) || (System.RDENS > 0) ) System.print_density(ticket);
    if(System.GYR_BINS > 0) System.print_gyration(ticket);
    
    System.dump_GRO_conf(name_out.str().c_str(), 0);
    System.dump_MOC_conf(name_MOCout.str().c_str(), 0);
    System.dump_GRO_CG_conf(name_CGout.str().c_str(), 0);
    System.dump_GRO_conf_unfolded(name_unfolded_out.str().c_str(), 0);


    logFile.close();
    KirkwFile.close();
    PolyFecFile.close();
    partFile.close();




    return 0;

}






















