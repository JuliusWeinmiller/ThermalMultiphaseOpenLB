
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

#ifndef MODIFIED_SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_2D_HH
#define MODIFIED_SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_2D_HH

#include "modifiedShanChenForcedSingleComponentPostProcessor2D.h"
#include "dynamics/interactionPotential.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "core/finiteDifference2D.h"

namespace olb {

////////  ShanChenForcedSingleComponentPostProcessor2D ///////////////////////////////////


template<typename T, typename DESCRIPTOR>
ModifiedShanChenForcedSingleComponentPostProcessor2D <T,DESCRIPTOR>::
ModifiedShanChenForcedSingleComponentPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
    T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_,
    std::vector<SpatiallyExtendedObject2D*> partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), G(G_), rho0(rho0_), interactionPotential(iP_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
ModifiedShanChenForcedSingleComponentPostProcessor2D <T,DESCRIPTOR>::
ModifiedShanChenForcedSingleComponentPostProcessor2D(T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_,
    std::vector<SpatiallyExtendedObject2D*> partners_)
  :  x0(0), x1(0), y0(0), y1(0), G(G_), rho0(rho0_), interactionPotential(iP_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
void ModifiedShanChenForcedSingleComponentPostProcessor2D<T,DESCRIPTOR>::
processSubDomain( BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                  int x0_, int x1_, int y0_, int y1_ )
{
  typedef DESCRIPTOR L;

  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( x0, x1, y0, y1,
                         x0_, x1_, y0_, y1_,
                         newX0, newX1, newY0, newY1 ) ) {
    int nx = newX1-newX0+3; // include a one-cell boundary
    int ny = newY1-newY0+3; // include a one-cell boundary
    int offsetX = newX0-1;
    int offsetY = newY0-1;

    BlockData2D<T,T> rhoField1(nx,ny);
    BlockData2D<T,T> psiField1(nx,ny);

    // Compute density and velocity on every site of first lattice, and store result
    //   in external scalars; envelope cells are included, because they are needed
    //   to compute the interaction potential in what follows.
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY);
        rhoField1.get(iX-offsetX, iY-offsetY) = cell.computeRho()*rho0[0];

        T psi;
        interactionPotential(&psi, &rhoField1.get(iX-offsetX, iY-offsetY) );
        psiField1.get(iX-offsetX, iY-offsetY) = psi;
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
          Cell<T,DESCRIPTOR>& blockCell   = blockLattice.get(iX,iY);

          // Compute local velocity using precalced rho
          T* j = blockCell.template getFieldPointer<descriptors::VELOCITY>();
          lbHelpers<T,DESCRIPTOR>::computeJ(blockCell,j);
          T rhoPreCalc = rhoField1.get(iX-offsetX, iY-offsetY);
          T *blockU = blockCell.template getFieldPointer<descriptors::VELOCITY>(); // contains precomputed value rho*u
          for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
            blockU[iD] = (blockU[iD]*rho0[0]) / rhoPreCalc;
          }

          
          // Computation of the interaction potential
          T rhoBlockContribution[L::d]   = {T(), T()};
          // T psi;
          // interactionPotential(&psi, &rhoField1.get(iX-offsetX, iY-offsetY) );
          T psi = psiField1.get(iX-offsetX, iY-offsetY);
          for (int iPop = 0; iPop < L::q; ++iPop) {
            int nextX = iX + descriptors::c<L>(iPop,0);
            int nextY = iY + descriptors::c<L>(iPop,1);
            // T blockPsi;
            // interactionPotential(&blockPsi, &rhoField1.get(nextX-offsetX, nextY-offsetY) );
            T blockPsi = psiField1.get(nextX-offsetX, nextY-offsetY);
            rhoBlockContribution[0]   += psi * blockPsi * descriptors::c<L>(iPop,0)* descriptors::t<T,L>(iPop);
            rhoBlockContribution[1]   += psi * blockPsi * descriptors::c<L>(iPop,1)* descriptors::t<T,L>(iPop);
          }

          // Computation and storage of the final velocity, consisting
          //   of u and the momentum difference due to interaction
          //   potential plus external force
          T *blockForce = blockCell.template getFieldPointer<descriptors::FORCE>();
          T *externalBlockForce = blockCell.template getFieldPointer<descriptors::EXTERNAL_FORCE>();

          for (int iD = 0; iD < L::d; ++iD) {
            blockForce[iD] = externalBlockForce[iD] - G*rhoBlockContribution[iD]/rhoField1.get(iX-offsetX, iY-offsetY);
          }

          // Storage of the psi value used for collide functions
          //T *psi_field = blockCell.template getFieldPointer<descriptors::PSI_PSEUDO_RHO>();
          //psi_field = psi;
          blockCell.template setField<descriptors::PSI_PSEUDO_RHO>(psi);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ModifiedShanChenForcedSingleComponentPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


/// LatticeCouplingGenerator for NS coupling

template<typename T, typename DESCRIPTOR>
ModifiedShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>::ModifiedShanChenForcedSingleComponentGenerator2D (
  int x0_, int x1_, int y0_, int y1_, T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), G(G_), rho0(rho0_), interactionPotential(iP_)
{ }

template<typename T, typename DESCRIPTOR>
ModifiedShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>::ModifiedShanChenForcedSingleComponentGenerator2D (
  T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(0, 0, 0, 0), G(G_), rho0(rho0_), interactionPotential(iP_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* ModifiedShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject2D*> partners) const
{
  return new ModifiedShanChenForcedSingleComponentPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, G, rho0, interactionPotential, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* ModifiedShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new ModifiedShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>(*this);
}



}  // namespace olb

#endif
