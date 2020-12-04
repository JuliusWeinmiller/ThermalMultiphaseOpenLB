
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

#ifndef INTERACTION_POTENTIAL_EXTEND_HH
#define INTERACTION_POTENTIAL_EXTEND_HH


#include "interactionPotentialExtended.h"


namespace olb {

template <typename T, typename S>
InteractionPotentialExtended<T,S>::InteractionPotentialExtended(T G, T a, T b, T R) : AnalyticalF1D<T,S>(1), _G(G), _a(a), _b(b), _R(R)
{ 

}

template <typename T, typename S>
T InteractionPotentialExtended<T,S>::getCriticalDensity(){
  return _rhoc;
}

template <typename T, typename S>
T InteractionPotentialExtended<T,S>::getCriticalTemperature(){
  return _tc;
}

template <typename T, typename S>
T InteractionPotentialExtended<T,S>::getG(){
  return _G;
}

template <typename T, typename S>
void InteractionPotentialExtended<T, S>::write(std::string const& title)
{
  std::string dataFile = singleton::directories().getLogOutDir() + title + ".dat";

  if (singleton::mpi().isMainProcessor())
  {
    std::ofstream fout;
    fout.open(dataFile.c_str(), std::ios::trunc);
    struct staticPressureParameters param = this->getParameters();

    fout << "----------------- InteractionPotentialExtended information -----------------" << std::endl;
    fout << "-- Parameters:" << std::endl;
    fout << "Equation of state name:            name=                       " << this->getName() << std::endl;
    fout << "Attraction factor:                 G=                          " << getG() << std::endl;
    fout << "Dimensionless gas constant:        R=                          " << param.R << std::endl;
    fout << "EoS parameter a:                   a=                          " << param.a << std::endl;
    fout << "EoS parameter b:                   b=                          " << param.b << std::endl;
    fout << "EoS parameter c:                   c=                          " << param.c << std::endl;
    fout << "EoS acentric Factor:               acentricFactor=             " << param.acentricFactor << std::endl;


    fout.close();
  }
}




template <typename T, typename S>
PengRobinsonExt<T,S>::PengRobinsonExt(T G, T acentricFactor, T a, T b, T tr, T R) : InteractionPotentialExtended<T,S>(G, a, b, R), _acentricFactor(acentricFactor)
{
  this->_R = R;
  //a=0.45724*R*R*tc*tc/pc;
  //b=0.0778*R*tc/pc;
  this->_tc = 0.0778/0.45724*this->_a/this->_b/this->_R;
  this->_pc = 0.0778*this->_R*this->_tc/this->_b;
  this->_rhoc = this->_pc/0.307/this->_R/this->_tc;
  // From original PR 0.307
  // From python file: 0.3074
  this->_t = this->_tc*tr;
  //Zc=0.307 Tc=0.072922004 pc=0.059569985 rhoc=2.6609121
  _c = (0.37464+1.54226*_acentricFactor-0.26992*_acentricFactor*_acentricFactor);
  _alpha = 1. + _c*(1.-sqrt(tr));
  _alpha = _alpha*_alpha;
  this->getName() = "PengRobinsonExt";
}

template <typename T, typename S>
bool PengRobinsonExt<T,S>::operator()(T psi[], const S rho[])
{
  T p = (rho[0]*this->_R*this->_t/(1.-this->_b*rho[0]))-(this->_a*_alpha*rho[0]*rho[0]/(1.+2.*this->_b*rho[0]-this->_b*this->_b*rho[0]*rho[0]));
  psi[0] = sqrt(6.*(p-rho[0]/3.)/this->_G);
  return true;
}

// second operator allows to incorporate temperature changes
template <typename T, typename S>
bool PengRobinsonExt<T,S>::operator()(T psi[], const S rho[], const S tr[])
{
  this->_t = tr[0]*this->_tc;
  _alpha = 1. + _c*(1.-sqrt(this->_t/this->_tc));
  _alpha = _alpha*_alpha;
  T p = (rho[0]*this->_R*this->_t/(1.-this->_b*rho[0]))-(this->_a*_alpha*rho[0]*rho[0]/(1.+2.*this->_b*rho[0]-this->_b*this->_b*rho[0]*rho[0]));
  psi[0] = sqrt(6.*(p-rho[0]/3.)/this->_G);
  return true;
}

template <typename T, typename S>
bool PengRobinsonExt<T,S>::computeDerivativeConstRho( T deriv[], const S rho[], const S tr[] )
{
  T t = tr[0]*this->_tc;
  T dalpha = _c*_c/this->_tc - _c*(1. + _c) /sqrt(t*this->_tc);
  //       = c/sqrt(t[0]*_tc)* (-1 - c + c*sqrt(t[0]/_tc) );

  deriv[0] = rho[0]*this->_R/(1. - this->_b*rho[0]) - this->_a * rho[0]*rho[0]*dalpha/(1.+2.*this->_b*rho[0]-this->_b*this->_b*rho[0]*rho[0]);
  // deriv[0] = rho[0]*this->_R/(1. - this->_b*rho[0]);
  return true;
}

//template <typename T, typename S>
//bool PengRobinsonExt<T,S>::computeDerivativeConstT( T deriv[], const S rho[], const S tr[] )
//{
//  T t = tr[0]*_tc;
//  T alpha = 1. + _c*(1.-sqrt(tr[0]));
//  alpha = alpha*alpha;
//
//  deriv[0] = -rho[0]*rho[0]*_R*t/(1. - _b*rho[0]) 
//             - _a * alpha * rho[0]*rho[0]/(1.+2.*_b*rho[0]-_b*_b*rho[0]*rho[0]);
//  return true;
//}

template <typename T, typename S>
double PengRobinsonExt<T,S>::computePressure ( double volume, void * parameters )
{
    struct staticPressureParameters * params = (staticPressureParameters *) parameters;

    double rho = 1./volume;
    double alpha = 1. + params->c*(1.-sqrt(params->t / params->tc ) );
    alpha = alpha*alpha;
    double pressure = (rho*params->R*params->t/(1.-params->b*rho))-(params->a*alpha*rho*rho/(1.+2.*params->b*rho-params->b*params->b*rho*rho));
    int sign = params->minimize ? 1 : -1;
    return sign * pressure;
}

template <typename T, typename S>
struct staticPressureParameters PengRobinsonExt<T,S>::getParameters ( )
{
  struct staticPressureParameters p = {};
  p.a=this->_a;
  p.b=this->_b;
  p.c=this->_c;
  p.R=this->_R;
  p.t=this->_t;
  p.tc=this->_tc;
  p.acentricFactor = _acentricFactor;
  return p;
}



template <typename T, typename S>
CarnahanStarlingExt<T,S>::CarnahanStarlingExt(T G, T a, T b, T tr, T R) : InteractionPotentialExtended<T,S>(G, a, b, R)
{
  this->_R = R;
  //a=0.4963*tc*tc*R*R/pc;
  //b=0.18727*R*tc/pc;
  this->_tc = 0.18727/0.4963*this->_a/this->_b/this->_R;
  this->_pc = 0.18727*this->_R*this->_tc/this->_b;
  this->_rhoc = this->_pc/0.35930763/this->_R/this->_tc;
  this->_t = this->_tc*tr;
  //Zc=0.35930763 Tc=0.094333065 pc=0.0044164383 rhoc=0.13029921
  this->getName() = "CarnahanStarlingExt";
}

template <typename T, typename S>
bool CarnahanStarlingExt<T,S>::operator()(T psi[], const S rho[])
{
  T c = this->_b*rho[0]/4.;
  T p = rho[0]*this->_R*this->_t*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))-this->_a*rho[0]*rho[0];
  psi[0] = sqrt(6.*(p-rho[0]/3.)/this->_G);
  return true;
}

// second operator allows to incorporate temperature changes
template <typename T, typename S>
bool CarnahanStarlingExt<T,S>::operator()(T psi[], const S rho[], const S tr[])
{
  this->_t = tr[0]*this->_tc;
  T c = this->_b*rho[0]/4.;
  T p = rho[0]*this->_R*this->_t*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))-this->_a*rho[0]*rho[0];
  psi[0] = sqrt(6.*(p-rho[0]/3.)/this->_G);
  return true;
}


template <typename T, typename S>
bool CarnahanStarlingExt<T,S>::computeDerivativeConstRho( T deriv[], const S rho[], const S tr[] )
{
  // This particular derivative is temperature independent
  T c = this->_b*rho[0]/4.;
  deriv[0] = rho[0]*this->_R*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c)) ;
  return true;
}

template <typename T, typename S>
double CarnahanStarlingExt<T,S>::computePressure ( double volume, void * parameters )
{
  struct staticPressureParameters * params = (staticPressureParameters *) parameters;
  double rho = 1./volume;
  double c = params->b*rho/4.;
  double pressure = rho*params->R*params->t*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))-params->a*rho*rho;
  int sign = params->minimize ? 1 : -1;
  return sign * pressure;
}

template <typename T, typename S>
struct staticPressureParameters CarnahanStarlingExt<T,S>::getParameters ( )
{
  struct staticPressureParameters p = {};
  p.a=this->_a;
  p.b=this->_b;
  p.R=this->_R;
  p.t=this->_t;
  p.tc=this->_tc;
  return p;
}


} // end namespace olb

#endif
