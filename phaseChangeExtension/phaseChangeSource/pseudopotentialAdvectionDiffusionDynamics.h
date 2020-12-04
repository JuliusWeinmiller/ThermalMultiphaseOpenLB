
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/** \file
 * An extention to the collection of dynamics classes with 
 * which a Cell object can be instantiated -- header file.
 */
#ifndef LB_PSEUDOPOTENTIAL_ADVECTION_DIFFUSION_DYNAMICS_H
#define LB_PSEUDOPOTENTIAL_ADVECTION_DIFFUSION_DYNAMICS_H

#include "dynamics/latticeDescriptors.h"
#include "core/util.h"
#include "core/postProcessing.h"
#include "core/latticeStatistics.h"


namespace olb{

/// Implementation of the BGK collision step with force correction factor
// Affect of the forcing term in the pseudopotential lattice Boltzmann meodeling of thermal flows
// Q. Li and K. Luo
template<typename T, typename DESCRIPTOR>
class PseudopotentialAdvectionDiffusionBGKdynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  PseudopotentialAdvectionDiffusionBGKdynamics( T omega, Momenta<T, DESCRIPTOR>& momenta );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
 
  T computeRho(Cell<T,DESCRIPTOR> const& cell) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d]) const override;
 
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;

protected:
  T _omega;  ///< relaxation parameter
  T _omegaMod;  // (1. - 0.5 * omega)

};


// ========= the MRT advection diffusion dynamics ========//
/// This approach is based on the multi-distribution LBM model.
template<typename T, typename DESCRIPTOR>
class PseudopotentialAdvectionDiffusionMRTdynamics : public BasicDynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  PseudopotentialAdvectionDiffusionMRTdynamics( T omega, Momenta<T, DESCRIPTOR>& momenta );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
private:
  T _omega;  ///< relaxation parameter
protected:
  T invM_S[DESCRIPTOR::q][DESCRIPTOR::q];                 ///< inverse relaxation times matrix
  T invM_correct[DESCRIPTOR::q][DESCRIPTOR::q];           ///< inverse M * (1-S/2)
};


} // end namespace

#endif
