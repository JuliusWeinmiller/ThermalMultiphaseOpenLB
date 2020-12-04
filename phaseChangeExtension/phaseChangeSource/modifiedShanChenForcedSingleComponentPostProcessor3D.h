
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/* Extention of the shan chen forced single component post processor
 * This post processor will include coupling to a thermal lattice
 *
*/

#ifndef MODIFIED_SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_3D_H
#define MODIFIED_SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_3D_H

#include "core/spatiallyExtendedObject3D.h"
#include "core/postProcessing.h"
#include "core/blockLattice3D.h"


namespace olb {

/**
* Multiphysics class for coupling between different lattices.
*/

// =========================================================================//
// ===========Shan Chen coupling without wall interaction===================//
// =========================================================================//

template<typename T, typename DESCRIPTOR>
class ModifiedShanChenForcedSingleComponentPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ModifiedShanChenForcedSingleComponentPostProcessor3D (
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T G_, std::vector<T> rho0_,
    AnalyticalF1D<T,T>& iP_, std::vector<SpatiallyExtendedObject3D*> partners_);
  ModifiedShanChenForcedSingleComponentPostProcessor3D (
    T G_, std::vector<T> rho0_,
    AnalyticalF1D<T,T>& iP_, std::vector<SpatiallyExtendedObject3D*> partners_);
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
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  T G;
  std::vector<T> rho0;
  AnalyticalF1D<T,T>& interactionPotential;
  std::vector<SpatiallyExtendedObject3D*> partners;
};


// This class is being called to generate the shan chen coupling

template<typename T, typename DESCRIPTOR>
class ModifiedShanChenForcedSingleComponentGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  ModifiedShanChenForcedSingleComponentGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_);
  ModifiedShanChenForcedSingleComponentGenerator3D(T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_);
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject3D*> partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  T G;
  std::vector<T> rho0;
  AnalyticalF1D<T,T>& interactionPotential;
};

}

#endif
