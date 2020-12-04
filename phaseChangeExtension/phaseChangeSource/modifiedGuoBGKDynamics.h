
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
#ifndef LB_PHASECHANGE_EXTENSION_DYNAMICS_H
#define LB_PHASECHANGE_EXTENSION_DYNAMICS_H

#include "dynamics/latticeDescriptors.h"
#include "dynamics/dynamics.h"
#include "core/util.h"
#include "core/postProcessing.h"
#include "core/latticeStatistics.h"



namespace olb{


/// Implementation of the BGK collision step with modified Guo force
// Forcing scheme in pseudopotential lattice Boltzmann model for multiphase flows
// Q. Li, K. Luo and X. Li
// 10.1103/PhysRevE.86.016709
template<typename T, typename DESCRIPTOR>
class ForcedModifiedGuoBGKdynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ForcedModifiedGuoBGKdynamics(T omega, Momenta<T,DESCRIPTOR>& momenta, T sigma);
  ///  Compute fluid velocity on the cell.
  void computeU (
    Cell<T,DESCRIPTOR> const& cell,
    T u[DESCRIPTOR::d] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,DESCRIPTOR> const& cell,
    T& rho, T u[DESCRIPTOR::d]) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;

  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override; 

  T getSigma() const;
  void setSigma(T sigma);

protected:
  T _sigma;  ///< correction factor
  T _omega;
  T _prefactor; // 
};








} // end namespace

#endif
