
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

#ifndef MAXWELL_CONSTRUCTION_H
#define MAXWELL_CONSTRUCTION_H

#include "interactionPotentialExtended.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>



namespace olb {

template <typename T, typename ipType>
class MaxwellConstruction {
public:
    // Constructor
    MaxwellConstruction( ipType & IP, T lowerBound=1e-1, T upperBound=1e20 );
    ~MaxwellConstruction();
    
    // Calculate densities from temperature ratio
    bool operator() ( T rho[], T tr );

    // Set new accuracy
    void setAccuracy(double newAccuracy);

    void setMaxIter(int newMaxIter);

    void setIntegrationSize(int newIntSize);

private:

    int maxIter = 200;
    double accuracy = 1e-6;

    ipType & _ip;
    T _lowerBound;
    T _upperBound;
    
    const gsl_min_fminimizer_type * _gsl_fmin_type;
    gsl_min_fminimizer * _gsl_fminimizer;
    gsl_function _gsl_fnMin;

    const gsl_root_fsolver_type * _gsl_root_type;
    gsl_root_fsolver * _gsl_root_fsolver;
    gsl_function _gsl_fnRoot;

    
    int _gsl_int_size = 1000;
    //gsl_function _gsl_fnInt;

    // GSL functions
    static double pressureOffset(double volume, void * param);
    struct pressureOffsetParameters {
        gsl_function gsl_func;
        double Poffset;
    };
    static double areaDifference(double pressure, void * param);
    struct areaDifferenceParameters {
        gsl_function gsl_func;
        double VminX, VmaxX;
        double lowerBound, upperBound;
        double accuracy;
        int maxIter;
        int _gsl_int_size;
        double intersectLeft = 0;
        double intersectRight = 0;
    };
};


} // end namespace olb
#endif
