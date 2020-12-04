
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

#ifndef MODIFIED_SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_2D_H
#define MODIFIED_SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_2D_H

#include "core/spatiallyExtendedObject2D.h"
#include "core/postProcessing.h"
#include "core/blockLattice2D.h"


namespace olb {

/**
* Multiphysics class for coupling between different lattices.
*/

// =========================================================================//
// ===========Shan Chen coupling without wall interaction===================//
// =========================================================================//

template<typename T, typename DESCRIPTOR>
class ModifiedShanChenForcedSingleComponentPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  ModifiedShanChenForcedSingleComponentPostProcessor2D (
    int x0_, int x1_, int y0_, int y1_, T G_, std::vector<T> rho0_,
    AnalyticalF1D<T,T>& iP_, std::vector<SpatiallyExtendedObject2D*> partners_);
  ModifiedShanChenForcedSingleComponentPostProcessor2D (
    T G_, std::vector<T> rho0_,
    AnalyticalF1D<T,T>& iP_, std::vector<SpatiallyExtendedObject2D*> partners_);
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
  int x0, x1, y0, y1;
  T G;
  std::vector<T> rho0;
  AnalyticalF1D<T,T>& interactionPotential;
  std::vector<SpatiallyExtendedObject2D*> partners;
};


// This class is being called to generate the shan chen coupling

template<typename T, typename DESCRIPTOR>
class ModifiedShanChenForcedSingleComponentGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  ModifiedShanChenForcedSingleComponentGenerator2D(int x0_, int x1_, int y0_, int y1_, T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_);
  ModifiedShanChenForcedSingleComponentGenerator2D(T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_);
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D*> partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;
private:
  T G;
  std::vector<T> rho0;
  AnalyticalF1D<T,T>& interactionPotential;
};

}

#endif
