
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
#ifndef LB_PHASECHANGE_EXTENSION_DYNAMICS_HH
#define LB_PHASECHANGE_EXTENSION_DYNAMICS_HH

#include <type_traits>
#include <numeric>
#include "modifiedGuoBGKDynamics.h"
#include "core/cell.h"


namespace olb{

////////////////////// Class ForcedModifiedGouBGKdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta Momenta object to know how to compute velocity momenta
 *  \param sigma correction factor, related to the thermodynamic consistancy
 */
template<typename T, typename DESCRIPTOR>
ForcedModifiedGuoBGKdynamics<T,DESCRIPTOR>::ForcedModifiedGuoBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta, T sigma )
  :  BasicDynamics<T,DESCRIPTOR>(momenta), _sigma(sigma), _omega(omega), _prefactor(_sigma / (1./_omega - 0.5) )
{
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::PSI_PSEUDO_RHO>() );
}

template<typename T, typename DESCRIPTOR>
void ForcedModifiedGuoBGKdynamics<T,DESCRIPTOR>::computeU (Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  const T* force = cell.template getFieldPointer<descriptors::FORCE>();
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedModifiedGuoBGKdynamics<T,DESCRIPTOR>::computeRhoU (Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  const T* force = cell.template getFieldPointer<descriptors::FORCE>();
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
}


template<typename T, typename DESCRIPTOR>
void ForcedModifiedGuoBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  // Prepping of collision
  T rho, u[DESCRIPTOR::d];                  
  this->_momenta.computeRhoU(cell, rho, u);

  const T* force = cell.template getFieldPointer<descriptors::FORCE>();
  const T* extForce = cell.template getFieldPointer<descriptors::EXTERNAL_FORCE>();
  const T psi   = cell.template getFieldPointer<descriptors::PSI_PSEUDO_RHO>()[0];

  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }

  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, _omega);
  
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += _prefactor * (force[iVel] - extForce[iVel]) *rho / psi / psi;
  }
  
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, _omega, rho);
  
  
  // closure of collision
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ForcedModifiedGuoBGKdynamics<T,DESCRIPTOR>::getSigma() const
{
  return _sigma;
}

template<typename T, typename DESCRIPTOR>
void ForcedModifiedGuoBGKdynamics<T,DESCRIPTOR>::setSigma(T sigma)
{
  _sigma = sigma;
  _prefactor = _sigma / (1./_omega - 0.5);
}

template<typename T, typename DESCRIPTOR>
T ForcedModifiedGuoBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void ForcedModifiedGuoBGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
  _prefactor = _sigma / (1./_omega - 0.5);
}



} // end namespace

#endif
