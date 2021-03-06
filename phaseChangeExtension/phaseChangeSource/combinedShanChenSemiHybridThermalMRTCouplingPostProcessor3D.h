
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
#ifndef COMBINED_SEMI_HYBRID_THERMAL_MRT_COUPLING_POST_PROCESSOR_3D_H
#define COMBINED_SEMI_HYBRID_THERMAL_MRT_COUPLING_POST_PROCESSOR_3D_H

#include <cmath>

#include "interactionPotentialExtended.h"
#include "fluidProperties.h"

#include "core/spatiallyExtendedObject3D.h"
#include "core/postProcessing.h"
#include "core/blockLattice3D.h"


namespace olb{


template<typename T, typename DESCRIPTOR>
class CombinedShanChenSemiHybridThermalMRTCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  CombinedShanChenSemiHybridThermalMRTCouplingPostProcessor3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T G_, T Tomega_,
    T averageRho_, T gravity_, std::vector<T> dir_, 
    FluidProperties<T>& fluidProperties_,
    InteractionPotentialExtended<T,T>& interactionPotential_,
    std::vector<SpatiallyExtendedObject3D* > partners_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1, z0, z1;
  T G, Tomega;
  T averageRho, gravity;
  std::vector<T> dir;

  T constant_k;
  T derivConductivity;
  T gravityForcePrefactor[L::d];
  
  FluidProperties<T>& fluidProperties;
  InteractionPotentialExtended<T,T>& interactionPotential;

  std::vector<SpatiallyExtendedObject3D*> partners;

  BlockLattice3D<T,descriptors::D3Q7<descriptors::tag::MRT,
                                     descriptors::VELOCITY,
                                     descriptors::PHI_THERMAL,
                                     descriptors::PREV_T_V,
                                     descriptors::PREV_PHI>> *tPartner;
};

template<typename T, typename DESCRIPTOR>
class CombinedShanChenSemiHybridThermalMRTCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  CombinedShanChenSemiHybridThermalMRTCouplingGenerator3D(
    int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_,
    T G_, T Tomega_,
    T averageRho_, T gravity_, std::vector<T> dir_,
    FluidProperties<T>& fluidProperties,
    InteractionPotentialExtended<T,T>& interactionPotential_ );
  PostProcessor3D<T,DESCRIPTOR>* generate(
    std::vector<SpatiallyExtendedObject3D* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;

private:
  T G, Tomega;
  T averageRho, gravity;
  std::vector<T> dir;
  FluidProperties<T>& fluidProperties;
  InteractionPotentialExtended<T,T>& interactionPotential;
};


} // end namespace

#endif
