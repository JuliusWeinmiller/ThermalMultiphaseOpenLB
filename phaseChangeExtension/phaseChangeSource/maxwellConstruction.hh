
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/** \file
 * Allows the calculation of the densities for vapour and liquid
 * via the maxwell construction of any EOS
 */

#ifndef MAXWELL_CONSTRUCTION_HH
#define MAXWELL_CONSTRUCTION_HH

#include "maxwellConstruction.h"


// Use GSL_min -> gsl_min_fminimizer_brent https://www.gnu.org/software/gsl/doc/html/min.html
// Use GSL_roots -> gsl_root_fsolver_brent   https://www.gnu.org/software/gsl/doc/html/roots.html
// fminimizer is the same as root of derivative
// Use gsl_integration https://www.gnu.org/software/gsl/doc/html/integration.html


namespace olb {

template <typename T, typename ipType>
MaxwellConstruction<T,ipType>::MaxwellConstruction( ipType & IP, T lowerBound, T upperBound ):
_ip(IP), _lowerBound(lowerBound), _upperBound(upperBound) 
{ 

    //_gsl_fmin_type = gsl_min_fminimizer_brent;
    _gsl_fmin_type = gsl_min_fminimizer_goldensection;
    _gsl_fminimizer = gsl_min_fminimizer_alloc(_gsl_fmin_type);

    _gsl_root_type = gsl_root_fsolver_bisection;
    _gsl_root_fsolver = gsl_root_fsolver_alloc(_gsl_root_type);
    
}

template <typename T, typename ipType>
MaxwellConstruction<T,ipType>::~MaxwellConstruction()
{
    gsl_min_fminimizer_free(_gsl_fminimizer);
    gsl_root_fsolver_free(_gsl_root_fsolver);
}

template <typename T, typename ipType>
void MaxwellConstruction<T,ipType>::setAccuracy(double newAccuracy)
{
    accuracy = newAccuracy;
}

template <typename T, typename ipType>
void MaxwellConstruction<T,ipType>::setMaxIter(int newMaxIter)
{
    maxIter = newMaxIter;
}

template <typename T, typename ipType>
void MaxwellConstruction<T,ipType>::setIntegrationSize(int newIntSize)
{
    _gsl_int_size = newIntSize;
}

template <typename T, typename ipType>
double MaxwellConstruction<T,ipType>::pressureOffset(double vol, void * parameters)
{
    struct pressureOffsetParameters * params = (pressureOffsetParameters *) parameters;
    //gsl_function _gsl_func = params.gsl_func;
    return GSL_FN_EVAL(&params->gsl_func, vol) - params->Poffset;
}

template <typename T, typename ipType>
double MaxwellConstruction<T,ipType>::areaDifference(double pressureGuess, void * parameters)
{
    olb::OstreamManager clout( std::cout,"MaxwellConstruction::areaDifference" );
    struct areaDifferenceParameters * params = (areaDifferenceParameters *) parameters;

  // Setting up GSL solvers
    const gsl_root_fsolver_type * _local_gsl_root_type;
    gsl_root_fsolver * _local_gsl_root_fsolver;
    gsl_function _local_gsl_fnRoot;

    _local_gsl_root_type = gsl_root_fsolver_brent;
    _local_gsl_root_fsolver = gsl_root_fsolver_alloc(_local_gsl_root_type);

    gsl_integration_workspace * _local_gsl_int_workspace;
    _local_gsl_int_workspace = gsl_integration_workspace_alloc(params->_gsl_int_size);
  // Above code is just setting up GSL solvers

  // unpacking parameters
    struct MaxwellConstruction<T,ipType>::pressureOffsetParameters offsetParameters;
    offsetParameters.gsl_func = params->gsl_func;
    offsetParameters.Poffset = pressureGuess;


// -- Finding the left and right intersect
    // Those intersects are used at integration bounds
    // Usage of a generic root solver to find both intersects
    auto performRootSolverGSL = [&] () -> double {
            int iter = 0; // reset iter
            int status;
            double a,b;
            do { // perform max search
                iter ++;
                status = gsl_root_fsolver_iterate(_local_gsl_root_fsolver);
                if(status == GSL_EBADFUNC || status == GSL_FAILURE){
                    throw status;
                }
                a = gsl_root_fsolver_x_lower (_local_gsl_root_fsolver);
                b = gsl_root_fsolver_x_upper (_local_gsl_root_fsolver);
                status = gsl_root_test_interval (a, b, 0.0, params->accuracy);
                
            }
            while( status == GSL_CONTINUE && iter < params->maxIter);
            return gsl_root_fsolver_root(_local_gsl_root_fsolver);
        };

    // Setting up intersect root solver
    _local_gsl_fnRoot.function = &pressureOffset;
    _local_gsl_fnRoot.params = &offsetParameters;


    // Finding left intersect
    gsl_root_fsolver_set( _local_gsl_root_fsolver, &_local_gsl_fnRoot, params->lowerBound, params->VminX );
    params->intersectLeft = performRootSolverGSL();
   
    // Finding right intersect
    gsl_root_fsolver_set( _local_gsl_root_fsolver, &_local_gsl_fnRoot, params->VmaxX, params->upperBound );
    params->intersectRight = performRootSolverGSL();
    
// -- Integrating the EOS pressure function with offset from left to right intersect    
    double areaIntegration, errorIntegration;
    // calculate area enclosed by line Pguess, intersect 1 and interset 3
    // clout << "Integration" << std::endl;
    gsl_integration_qags(&_local_gsl_fnRoot, params->intersectLeft, params->intersectRight, 0.0, params->accuracy, params->_gsl_int_size, _local_gsl_int_workspace, &areaIntegration, &errorIntegration);
    
// Cleanup and return
    gsl_root_fsolver_free(_local_gsl_root_fsolver);
    gsl_integration_workspace_free(_local_gsl_int_workspace);
    
    return areaIntegration;
}

template <typename T, typename ipType>
bool MaxwellConstruction<T,ipType>::operator() ( T rho[2], T tr)
{
    olb::OstreamManager clout( std::cout,"MaxwellConstruction" );
    //clout << "Starting MaxwellConstruction" << std::endl;
    // Get point to the left and right of min max
    T Vc = 1./_ip.getCriticalDensity();
    //clout << "Critical Volume: "<< Vc << std::endl;
    struct staticPressureParameters IPparams = _ip.getParameters();
    IPparams.t = IPparams.tc * tr;


    // convergence lambda
    auto performMinSolverGSL = [&,this] () {
        int iter = 0; // reset iter
        int status;
        double a,b;
        do { // perform max search
            iter ++;
            status = gsl_min_fminimizer_iterate(_gsl_fminimizer);
            if(status == GSL_EBADFUNC || status == GSL_FAILURE){
                throw status;
            }
            a = gsl_min_fminimizer_x_lower (_gsl_fminimizer);
            b = gsl_min_fminimizer_x_upper (_gsl_fminimizer);
            status = gsl_min_test_interval (a, b, accuracy, 0.0); 
        }
        while( status == GSL_CONTINUE && iter < maxIter);
    };

    
    //clout << "Finding Minimum" << std::endl;
    _gsl_fnMin.function = &_ip.computePressure;
// Get minimum X
    _gsl_fnMin.params = &IPparams;
    
    //clout << "Left bound : " << _lowerBound << " : " << GSL_FN_EVAL(&_gsl_fnMin, _lowerBound)  << std::endl;
    //clout << "Init Guess : " << Vc * 0.9 << " : " << GSL_FN_EVAL(&_gsl_fnMin, Vc * 0.9)  << std::endl;
    //clout << "Right bound: " << Vc << " : " << GSL_FN_EVAL(&_gsl_fnMin, Vc)  << std::endl;

    gsl_min_fminimizer_set(_gsl_fminimizer, &_gsl_fnMin, Vc*0.9, _lowerBound, Vc );
    performMinSolverGSL();
    double VminX  = gsl_min_fminimizer_x_minimum (_gsl_fminimizer);
    double PressureMin = gsl_min_fminimizer_f_minimum (_gsl_fminimizer);

    // clout << "Finding Maximum" << std::endl;
// Get maximum x
    IPparams.minimize = false;
    _gsl_fnMin.params = &IPparams;
    // clout << "Left bound : " << VminX << " : " << GSL_FN_EVAL(&_gsl_fnMin, VminX)  << std::endl;
    // clout << "Init Guess : " << _upperBound/10 << " : " << GSL_FN_EVAL(&_gsl_fnMin, _upperBound/10)  << std::endl;
    // clout << "Right bound: " << _upperBound << " : " << GSL_FN_EVAL(&_gsl_fnMin, _upperBound)  << std::endl;
    gsl_min_fminimizer_set(_gsl_fminimizer, &_gsl_fnMin, _upperBound/10, VminX, _upperBound );
    // clout << "Perform Maximum" << std::endl;
    performMinSolverGSL();
    double VmaxX  = gsl_min_fminimizer_x_minimum (_gsl_fminimizer);
    double PressureMax = -gsl_min_fminimizer_f_minimum (_gsl_fminimizer);  
    
    // clout << "VminX : " << VminX << std::endl;
    // clout << "Vc : " << Vc << std::endl;
    // clout << "VmaxX : " << VmaxX << std::endl;

// Perform Maxwell Construction
    // - Find equilibrium Pressure such that 
    //       integral between intersects 1&2 is the same as 
    //       integral between intersects 2&3


    struct areaDifferenceParameters parametersArea;
    IPparams.minimize = true;
    _gsl_fnMin.params = &IPparams;
    parametersArea.gsl_func = _gsl_fnMin; 
    parametersArea.VminX = VminX;
    parametersArea.VmaxX = VmaxX;
    parametersArea.lowerBound = _lowerBound;
    parametersArea.upperBound = _upperBound;
    parametersArea.maxIter = maxIter;
    parametersArea.accuracy = accuracy;
    parametersArea._gsl_int_size = _gsl_int_size;

    double lowP = PressureMin < 0 ? 1e-20 : PressureMin;
    //double initP = (lowP+PressureMax)/2;


// This was the original method, using bisectional search to minimize area
// IF REUSING THIS METHOD, ENSURE TO PUT THE std::fabs(area) BACK
    // // clout << "AreaDifference Minimization" << std::endl;
    // gsl_function _gsl_fnArea;
    // _gsl_fnArea.function = &areaDifference;
    // _gsl_fnArea.params = &parametersArea;
    // // clout << "Init AreaDifference" << std::endl;
    // clout << "Left bound : " << lowP << " : " << GSL_FN_EVAL(&_gsl_fnArea, lowP)  << std::endl;
    // clout << "Init Guess : " << initP << " : " << GSL_FN_EVAL(&_gsl_fnArea, initP)  << std::endl;
    // clout << "Right bound: " << PressureMax << " : " << GSL_FN_EVAL(&_gsl_fnArea, PressureMax)  << std::endl;
    // clout << "Perform AreaDifference" << std::endl;
    // gsl_min_fminimizer_set(_gsl_fminimizer, &_gsl_fnArea, initP , lowP , PressureMax );
    // clout << "Performing AreaDifference" << std::endl;
    // performMinSolverGSL();


// Setting up GSL solvers
    

    // clout << "AreaDifference Root finder" << std::endl;
// -- Finding the area intersect -> areaDifference = 0
    // Usage of a generic root solver
    auto performRootSolverGSL = [&] () -> double {
            int iter = 0; // reset iter
            int status;
            double a,b;
            do { // perform max search
                iter ++;
                status = gsl_root_fsolver_iterate(_gsl_root_fsolver);
                if(status == GSL_EBADFUNC || status == GSL_FAILURE){
                    throw status;
                }
                a = gsl_root_fsolver_x_lower (_gsl_root_fsolver);
                b = gsl_root_fsolver_x_upper (_gsl_root_fsolver);
                status = gsl_root_test_interval (a, b, 0.0, parametersArea.accuracy);
                
            }
            while( status == GSL_CONTINUE && iter < parametersArea.maxIter);
            return gsl_root_fsolver_root(_gsl_root_fsolver);
        };

    // Setting up intersect root solver
    _gsl_fnRoot.function = &areaDifference;
    _gsl_fnRoot.params = &parametersArea;
    // Perform area solving
    gsl_root_fsolver_set( _gsl_root_fsolver, &_gsl_fnRoot, lowP, PressureMax );
    performRootSolverGSL();

    // Get left and right intersect from params
    rho[0] = 1./ parametersArea.intersectLeft;
    rho[1] = 1./ parametersArea.intersectRight;
    return true;

}



}
#endif
