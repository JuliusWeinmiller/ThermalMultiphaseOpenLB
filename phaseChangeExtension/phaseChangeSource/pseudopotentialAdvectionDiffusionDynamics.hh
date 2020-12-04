
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/** \file
 * An extension to the collection of dynamics classes with 
 * which a Cell object can be instantiated -- generic implementation.
 */
#ifndef LB_PSEUDOPOTENTIAL_ADVECTION_DIFFUSION_DYNAMICS_HH
#define LB_PSEUDOPOTENTIAL_ADVECTION_DIFFUSION_DYNAMICS_HH

#include <type_traits>
#include <numeric>
#include "dynamics/advectionDiffusionDynamics.h"
#include "pseudopotentialAdvectionDiffusionDynamics.h"
#include "core/cell.h"


namespace olb{

////////////////////// Class PseudopotentialAdvectionDiffusionBGKdynamics /////////////////////////

/** \param omega relaxation parameter, related to the thermal diffusivity
 *  \param momenta Momenta object to know how to compute thermal momenta
 */

template<typename T, typename DESCRIPTOR>
PseudopotentialAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::PseudopotentialAdvectionDiffusionBGKdynamics (
  T omega, Momenta<T, DESCRIPTOR>& momenta )
  : BasicDynamics<T, DESCRIPTOR>( momenta ),
    _omega(omega), _omegaMod(1. - 0.5 * omega)
{ }

template<typename T, typename DESCRIPTOR>
T PseudopotentialAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::computeEquilibrium( int iPop, T rho,
    const T u[DESCRIPTOR::d], T uSqr ) const
{
  return lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, typename DESCRIPTOR>
void PseudopotentialAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::collide( Cell<T, DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  
  const T* u = cell.template getFieldPointer<olb::descriptors::VELOCITY>();      // get pointer to the first velocity array location
  const T phi = *cell.template getFieldPointer<olb::descriptors::PHI_THERMAL>(); // since it isn't an array. get the value at the pointer location
  T temperature = this->_momenta.computeRho( cell ) ;
  
  const T* Tv0 = cell.template getFieldPointer<olb::descriptors::PREV_T_V>(); // get pointer to previous timestep Tv
  const T phi0 = *cell.template getFieldPointer<olb::descriptors::PREV_PHI>(); // get pointer to previous timestep Phi

  T Tv1[ DESCRIPTOR::d];
  T derivTv[ DESCRIPTOR::d];
  for (int iD=0; iD< DESCRIPTOR::d; ++iD){
    Tv1[iD] = temperature * u[iD];              // calculate current Tv
    derivTv[iD] = Tv1[iD] - Tv0[iD];                // calculate derivative dTv <- Tv () (time step dt is 1)
  }

  T derivPhi = phi - phi0;


  // For the D2Q5 and the D3Q7, the bgkCollision function already uses the first order equilibrium calculations
  // Thus no changes are required
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision( cell, temperature, u, _omega );

  T cdTv;
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    T weight = descriptors::t<T,DESCRIPTOR>(iPop);
  
  //  Correct for phase change
    // Perform corrections according to 10.1103 / 96.063303
    cell[iPop] += weight * phi;
  
  //  Correct for dt phi
    // Perform corrections according to 10.1103 / 96.063303
    cell[iPop] += weight * 0.5 * derivPhi ;
  
  //  Correct for dt Tv
    // Perform corrections according to 10.1103 / 96.063303
    cdTv = 0.0;
    for (int iD=0; iD< DESCRIPTOR::d; ++iD){
      cdTv += descriptors::c<DESCRIPTOR>(iPop,iD) * derivTv[iD];
    }
    cell[iPop] += _omegaMod * weight * cdTv * descriptors::invCs2<T,DESCRIPTOR>();

  }

  cell.template defineField<olb::descriptors::PREV_T_V>(Tv1); // set Field to current timestep Tv
  cell.template defineField<olb::descriptors::PREV_PHI>(&phi); // set Field to current timestep phi

  statistics.incrementStats( temperature, uSqr );
}

template<typename T, typename DESCRIPTOR>
T PseudopotentialAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>::computeRho(Cell<T,DESCRIPTOR> const& cell) const
{
  return this->_momenta.computeRho( cell ) ;
}

template<typename T, typename DESCRIPTOR>
void PseudopotentialAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>::computeRhoU (
  Cell<T,DESCRIPTOR> const& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  this->_momenta.computeRhoU( cell, rho, u );
}


template<typename T, typename DESCRIPTOR>
T PseudopotentialAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void PseudopotentialAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::setOmega( T omega )
{
  _omega = omega;
  _omegaMod = (1. - 0.5 * omega);
}



//==================================================================//
//================= MRT Model for Advection diffusion ==============//
//==================================================================//


template<typename T, typename DESCRIPTOR>
PseudopotentialAdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::PseudopotentialAdvectionDiffusionMRTdynamics(
  T omega, Momenta<T, DESCRIPTOR>& momenta) :
  BasicDynamics<T, DESCRIPTOR>(momenta), _omega(omega)
{
  
  T rt[DESCRIPTOR::q]; // relaxation times vector.
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    rt[iPop] = descriptors::s<T,DESCRIPTOR>(iPop);
  }
  for (int iPop = 0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
    rt[descriptors::shearViscIndexes<DESCRIPTOR>(iPop)] = omega;
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      invM_S[iPop][jPop] = T();
      for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
        if (kPop == jPop) {
          invM_S[iPop][jPop] += descriptors::invM<T,DESCRIPTOR>(iPop,kPop) * rt[kPop];
        }
      }
    }
  }

  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      invM_correct[iPop][jPop] = T();
      for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
        if (kPop == jPop) {
          invM_correct[iPop][jPop] += descriptors::invM<T,DESCRIPTOR>(iPop,kPop) * (1 - rt[kPop]/2. );
        }
      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
T PseudopotentialAdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::computeEquilibrium(int iPop, T rho,
    const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder(iPop, rho, u);
}

template<typename T, typename DESCRIPTOR>
void PseudopotentialAdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::collide(Cell<T, DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics)
{
  // OstreamManager clout( std::cout,"MRT collide" );

  T temperature = this->_momenta.computeRho(cell);
  const T* u = cell.template getFieldPointer<descriptors::VELOCITY>();
  const T phi = *cell.template getFieldPointer<descriptors::PHI_THERMAL>();

  const T* Tv0 = cell.template getFieldPointer<olb::descriptors::PREV_T_V>(); // get pointer to previous timestep Tv
  const T phi0 = *cell.template getFieldPointer<olb::descriptors::PREV_PHI>(); // get pointer to previous timestep Phi

  T Tv1[ DESCRIPTOR::d];
  T derivTv[ DESCRIPTOR::d];
  for (int iD=0; iD< DESCRIPTOR::d; ++iD){
    Tv1[iD] = temperature * u[iD];                // calculate current Tv
    derivTv[iD] = Tv1[iD] - Tv0[iD];              // calculate derivative dTv <- Tv () (time step dt is 1)
  }
  T derivPhi = phi - phi0;


  T uSqr = lbHelpers<T, DESCRIPTOR>::mrtCollision(cell, temperature, u, invM_S);
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop){
    cell[iPop] += descriptors::invM<T, DESCRIPTOR>(iPop,0)*(phi + 0.5 * derivPhi);
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop){
    for (int iShear = 0; iShear < descriptors::shearIndexes<DESCRIPTOR>(); ++iShear) {
      cell[iPop] += invM_correct[iPop][descriptors::shearViscIndexes<DESCRIPTOR>(iShear)] * derivTv[iShear] ;
    }
  }


  cell.template defineField<olb::descriptors::PREV_T_V>(Tv1); // set Field to current timestep Tv
  cell.template defineField<olb::descriptors::PREV_PHI>(&phi); // set Field to current timestep phi

  statistics.incrementStats(temperature, uSqr);
}

template<typename T, typename DESCRIPTOR>
T PseudopotentialAdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::getOmega() const
{   return _omega; }

template<typename T, typename DESCRIPTOR>
void PseudopotentialAdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::setOmega(T omega)
{   _omega = omega; }



} // end namespace

#endif
