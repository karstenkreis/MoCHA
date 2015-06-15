#ifndef _MYRAND
#define	_MYRAND

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
//#include "mocha_system.h"

namespace mocha_myrandom_nmspc
{

    using namespace std;
    
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int> dist13(0,12);
    uniform_int_distribution<int> dist16(0,15);              // Hardcoded for 16 atoms...
    uniform_real_distribution<double> ddist(-1.0,1.0);
    uniform_real_distribution<double> ddist2(0.0,1.0);
    uniform_int_distribution<int> distNmol(0,4963);          // Hardcoded for 4964 molecules...
    //uniform_int_distribution<int> distNmol(0,6618);
    //uniform_int_distribution<int> distNmol(0,8272);
    //uniform_int_distribution<int> distNmol(0,9927);
    //uniform_int_distribution<int> distNmol(0,359);
    //uniform_int_distribution<int> distNmol(0,1356);
    //uniform_int_distribution<int> distNmol(0,827);
    //uniform_int_distribution<int> distNmol(0,9099);
    
    double gasdev()
    {
        static int iset = 0;
        static double gset;
        double fac, rsq, v1, v2;

        if (iset == 0) {
            do {
                v1 = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                v2 = ddist(mt);//2.0*(rand() / double(RAND_MAX)) - 1.0;
                if( (v1 < -1.0) || (v1 > 1.0) ||(v2 < -1.0) || (v2 > 1.0) )
                {
                    cerr << "Error in gasdev: v1 = " << v1 << " v2 = " << v2 << endl;
                    exit(1);
                }
                rsq = v1 * v1 + v2*v2;
            } while (rsq >= 1.0 || rsq == 0.0);
            fac = sqrt(-2.0 * log(rsq) / rsq);
            gset = v1*fac;
            iset = 1;
            return v2*fac;
        } else {
            iset = 0;
            return gset;
        }
    }


}

#endif
