
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/** \file
 * A modification to the navier stokes into advection diffusion coupling 
 * -- header file.
 */
#ifndef COMBINED_SEMI_HYBRID_THERMAL_MRT_COUPLING_POST_PROCESSOR_2D_H
#define COMBINED_SEMI_HYBRID_THERMAL_MRT_COUPLING_POST_PROCESSOR_2D_H

#include <cmath>

#include "interactionPotentialExtended.h"
#include "fluidProperties.h"

#include "core/spatiallyExtendedObject2D.h"
#include "core/postProcessing.h"
#include "core/blockLattice2D.h"

namespace olb{


template<typename T, typename DESCRIPTOR>
class CombinedShanChenSemiHybridThermalMRTCouplingPostProcessor2D :
  public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  CombinedShanChenSemiHybridThermalMRTCouplingPostProcessor2D(
    int x0_, int x1_, int y0_, int y1_,
    T G_, T Tomega_,
    T averageRho_, T gravity_, std::vector<T> dir_, 
    FluidProperties<T>& fluidProperties_,
    InteractionPotentialExtended<T,T>& interactionPotential_,
    std::vector<SpatiallyExtendedObject2D* > partners_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1;
  T G, Tomega;
  T averageRho, gravity;
  std::vector<T> dir;
  
  T constant_k;
  T derivConductivity;
  T gravityForcePrefactor[L::d];
  
  FluidProperties<T>& fluidProperties;
  InteractionPotentialExtended<T,T>& interactionPotential;

  std::vector<SpatiallyExtendedObject2D*> partners;

  BlockLattice2D<T,descriptors::D2Q5<descriptors::tag::MRT,
                                     descriptors::VELOCITY,
                                     descriptors::PHI_THERMAL,
                                     descriptors::PREV_T_V,
                                     descriptors::PREV_PHI>> *tPartner;

};

template<typename T, typename DESCRIPTOR>
class CombinedShanChenSemiHybridThermalMRTCouplingGenerator2D :
  public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  CombinedShanChenSemiHybridThermalMRTCouplingGenerator2D(
    int x0_, int x1_, int y0_, int y1_,
    T G_, T Tomega_,
    T averageRho_, T gravity_, std::vector<T> dir_,
    FluidProperties<T>& fluidProperties,
    InteractionPotentialExtended<T,T>& interactionPotential_ );
  PostProcessor2D<T,DESCRIPTOR>* generate(
    std::vector<SpatiallyExtendedObject2D* > partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;

private:
  T G, Tomega;
  T averageRho, gravity;
  std::vector<T> dir;
  FluidProperties<T>& fluidProperties;
  InteractionPotentialExtended<T,T>& interactionPotential;
};


} // end namespace

#endif
