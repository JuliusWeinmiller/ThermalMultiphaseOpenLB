
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

#ifndef INTERACTION_POTENTIAL_EXTEND_H
#define INTERACTION_POTENTIAL_EXTEND_H


#include "dynamics/interactionPotential.h"


namespace olb {

struct staticPressureParameters {
  double a=0;
  double b=0;
  double c=0;
  double R=0;
  double t=0;
  double tc=0;
  double acentricFactor=0;
  bool minimize = true;
};

template <typename T, typename S>
class InteractionPotentialExtended: public AnalyticalF1D<T,S> {
public:
  InteractionPotentialExtended(T G, T a, T b, T R=1);
  // first operator allows to calculate psi from rho
  virtual bool operator() (T psi[], const S rho[]) = 0;
  // second operator allows to incorporate temperature changes
  virtual bool operator() (T psi[], const S rho[], const S tr[]) = 0;
  // special function to get derivative at constant rho for thermal coupling 
  virtual bool computeDerivativeConstRho( T deriv[], const S rho[], const S tr[]) = 0;
  // special function to get derivative at constant temp for maxwell constrcution 
  //virtual bool computeDerivativeConstT( T deriv[], const S rho[], const S tr[]) = 0;
  virtual struct staticPressureParameters getParameters();

  // function to get _rhoc
  T getCriticalDensity();
  T getCriticalTemperature();

  // function to save to .dat file
  void write(std::string const& title = "InteractionPotential");

private:
  T getG();

protected:
  T _G;
  T _a;
  T _b;
  T _R;
  T _t;
  T _tc;
  T _pc;
  T _rhoc;

};


template <typename T, typename S>
class PengRobinsonExt : public InteractionPotentialExtended<T,S> {
private:
  T _acentricFactor;
  T _c;
  T _alpha;
public:
  PengRobinsonExt(T G, T acentricFactor=0.334, T a=2./49., T b=2./21., T tr=.7, T R=1.);
  bool operator() (T psi[], const S rho[]) override;
  // second operator allows to incorporate temperature changes
  bool operator() (T psi[], const S rho[], const S tr[]) override;
  // special function to get derivative at constant rho for thermal coupling 
  bool computeDerivativeConstRho( T deriv[], const S rho[], const S tr[]) override;
  // special function to get derivative at constant temp for maxwell constrcution 
  //bool computeDerivativeConstT( T deriv[], const S rho[], const S traits[]) override;
  // function to compute pressure, used to calculate Maxwell Construction
  static double computePressure( double volume, void * params );
  struct staticPressureParameters getParameters();

};


template <typename T, typename S>
class CarnahanStarlingExt : public InteractionPotentialExtended<T,S> {
public:
  CarnahanStarlingExt(T G, T a=1., T b=4., T tr=.7, T R=1.);
  bool operator() (T psi[], const S rho[]) override;
  // second operator allows to incorporate temperature changes
  bool operator() (T psi[], const S rho[], const S tr[]) override;
  // special function to get derivative at constant rho for thermal coupling 
  bool computeDerivativeConstRho( T deriv[], const S rho[], const S tr[]) override;
  // special function to get derivative at constant temp for maxwell constrcution 
  //bool computeDerivativeConstT( T deriv[], const S rho[], const S tr[]) override;
  // function to compute pressure, used to calculate Maxwell Construction
  static double computePressure( double volume, void * params );
  struct staticPressureParameters getParameters();
};

} // end olb namespace

#endif
