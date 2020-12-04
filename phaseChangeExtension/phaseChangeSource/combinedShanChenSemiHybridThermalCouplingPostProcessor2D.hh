
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/** \file
 * A modification to the navier stokes into advection diffusiton coupling 
 * -- generic implementation.
 */
#ifndef COMBINED_SEMI_HYBRID_THERMAL_COUPLING_POST_PROCESSOR_2D_HH
#define COMBINED_SEMI_HYBRID_THERMAL_COUPLING_POST_PROCESSOR_2D_HH

#include "combinedShanChenSemiHybridThermalCouplingPostProcessor2D.h"
#include "dynamics/latticeDescriptors.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "core/finiteDifference2D.h"


namespace olb{


//=====================================================================================
//==============  NavierStokesAdvectionDiffusionCouplingPostProcessor2D ===========
//=====================================================================================

template<typename T, typename DESCRIPTOR>
CombinedShanChenSemiHybridThermalCouplingPostProcessor2D<T,DESCRIPTOR>::
CombinedShanChenSemiHybridThermalCouplingPostProcessor2D(
    int x0_, int x1_, int y0_, int y1_,
    T G_, T Tomega_,
    T averageRho_, T gravity_, std::vector<T> dir_,
    FluidProperties<T>& fluidProperties_,
    InteractionPotentialExtended<T,T>& interactionPotential_,
    std::vector<SpatiallyExtendedObject2D* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
     G(G_), Tomega(Tomega_),
     averageRho(averageRho_), gravity(gravity_), dir(dir_),
     fluidProperties(fluidProperties_), interactionPotential(interactionPotential_),
     partners(partners_)
{

  constant_k = ( 1. / Tomega - 0.5) / descriptors::invCs2<T, descriptors::D2Q5<>>() ;
  derivConductivity = 
        (fluidProperties.conductivity(fluidProperties.getDensityLiquid() ) - fluidProperties.conductivity(fluidProperties.getDensityVapour() ) )
      / (fluidProperties.getDensityLiquid()                                - fluidProperties.getDensityVapour() );

  tPartner = dynamic_cast<BlockLattice2D<T,descriptors::D2Q5<descriptors::VELOCITY, descriptors::PHI_THERMAL, descriptors::PREV_T_V, descriptors::PREV_PHI>> *>(partners[0]);

  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD]*dir[iD];
  }
  normDir = sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    gravityForcePrefactor[iD] = gravity * dir[iD];
  }

}

template<typename T, typename DESCRIPTOR>
void CombinedShanChenSemiHybridThermalCouplingPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int nx = newX1-newX0+3; // include a one-cell boundary
    int ny = newY1-newY0+3; // include a one-cell boundary
    int offsetX = newX0-1;
    int offsetY = newY0-1;

    BlockData2D<T,T> rhoField(nx,ny);
    BlockData2D<T,T> psiField(nx,ny);
    BlockData2D<T,T> uField(nx,ny, DESCRIPTOR::d);
    BlockData2D<T,T> tempField(nx,ny);

    // Compute density, velocity and temperature on every site of lattices, and store result
    //   in external scalars; envelope cells are included, because they are needed
    //   to compute the interaction potential in what follows.
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {

        // Fluid calculation of density and velocity
        Cell<T,DESCRIPTOR>& fluidCell = blockLattice.get(iX,iY);
        T rho, u[DESCRIPTOR::d];
        fluidCell.computeRhoU(rho, u);
        rhoField.get(iX-offsetX, iY-offsetY) = rho;
        for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
          uField.get(iX-offsetX, iY-offsetY, iD) = u[iD];
        }

        // Temperature calculation of density
        T temperature, temperatureRatio;
        temperature = tPartner->get(iX,iY).computeRho();
        temperatureRatio = fluidProperties.temperatureRatio(temperature);
        tempField.get(iX-offsetX, iY-offsetY) = temperature;

        // Storage of the psi value used in Li's BGK dynamics collide functions
        // Also used in calculations of pseudopotential force
        T psi;
        interactionPotential(&psi, 
                              &rho, 
                              &temperatureRatio );
        psiField.get(iX-offsetX, iY-offsetY) = psi;
        fluidCell.template setField<descriptors::PSI_PSEUDO_RHO>(psi);

      }
    }


    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
  

// Computation of the pseudopotential force
        Cell<T,DESCRIPTOR>& fluidCell   = blockLattice.get(iX,iY);

        T localRho = rhoField.get(iX-offsetX, iY-offsetY);

        T rhoBlockContribution[L::d]   = {T(), T()};
        const T psi = psiField.get(iX-offsetX, iY-offsetY);

        for (int iPop = 0; iPop < L::q; ++iPop) {
          int nextX = iX + descriptors::c<L>(iPop,0);
          int nextY = iY + descriptors::c<L>(iPop,1);

          const T blockPsi = psiField.get(nextX-offsetX, nextY-offsetY);

          rhoBlockContribution[0]   += psi * blockPsi * descriptors::c<L>(iPop,0)* descriptors::t<T,L>(iPop);
          rhoBlockContribution[1]   += psi * blockPsi * descriptors::c<L>(iPop,1)* descriptors::t<T,L>(iPop);
        }

        // Computation and storage of the final velocity, consisting
        //   of u and the momentum difference due to interaction
        //   potential plus external force
        T *fluidForce = fluidCell.template getFieldPointer<descriptors::FORCE>();
        T *externalFluidForce = fluidCell.template getFieldPointer<descriptors::EXTERNAL_FORCE>();

// Bouyancy
        // Multiphase coupling executed first
        // Reset external Fluid force
        // Add MPF to force
        T localBouyancyMagnitude = ( 1 - averageRho/localRho); 
        for (unsigned iD = 0; iD < L::d; ++iD) {
          T tempBouyancyForce = gravityForcePrefactor[iD] * localBouyancyMagnitude; 
          externalFluidForce[iD] = tempBouyancyForce;
          fluidForce[iD] =  tempBouyancyForce - G*rhoBlockContribution[iD]/localRho;
        }



// Computation of the thermal force

        // reset phi_Li
        T phi_li = 0; 
        T localTemp = tempField.get(iX-offsetX, iY-offsetY);
        T localTempRatio = fluidProperties.temperatureRatio(localTemp);

        T localSpecificHeat = fluidProperties.specificHeat(localRho);
        T invVolumetricHeat = 1. / localRho / localSpecificHeat ;
        
        // second derivative of T at iX, iY
        // using finite difference
        // Units = [k / m^2]
        T nablaNablaT =  tempField.get(iX-offsetX-1, iY-offsetY)
                        +tempField.get(iX-offsetX+1, iY-offsetY)
                        +tempField.get(iX-offsetX, iY-offsetY-1)
                        +tempField.get(iX-offsetX, iY-offsetY+1)
                        -4*localTemp;
                        
        // First derivative of T at iX, iY
        // using central difference scheme
        // Units = [k / m]
        T nablaT =  (-tempField.get(iX-offsetX-1, iY-offsetY)
                     +tempField.get(iX-offsetX+1, iY-offsetY)
                     -tempField.get(iX-offsetX, iY-offsetY-1)
                     +tempField.get(iX-offsetX, iY-offsetY+1))*0.5;
        
        // First derivative of rho at iX, iY
        // using central difference scheme
        // Units = [kg / m^4]
        T nablaRho = (-rhoField.get(iX-offsetX-1, iY-offsetY)
                      +rhoField.get(iX-offsetX+1, iY-offsetY)
                      -rhoField.get(iX-offsetX, iY-offsetY-1)
                      +rhoField.get(iX-offsetX, iY-offsetY+1))*0.5;

        // First derivative of U ( nablaU ) at iX, iY
        // using central difference scheme
        // Units = [1 / s]
        T nablaU = ( uField.get(iX-offsetX+1, iY-offsetY    , 0)
                    -uField.get(iX-offsetX-1, iY-offsetY    , 0)
                    +uField.get(iX-offsetX  , iY-offsetY+1  , 1)
                    -uField.get(iX-offsetX  , iY-offsetY-1  , 1))*0.5;


    // First part: custom diffusion
        // invVolumetricHeat * conductivity = Diffusivity [m^2 / s]
        // https://en.wikipedia.org/wiki/Thermal_diffusivity
        

        // Remove LBM diffusion
        phi_li -= constant_k*nablaNablaT;

        // Add custom diffusion with density dependent diffusivity
        phi_li += invVolumetricHeat * fluidProperties.conductivity(localRho) * nablaNablaT;
        phi_li += invVolumetricHeat * derivConductivity * nablaRho * nablaT;

    // Second part: phase change
        // => Energy taken by phase change
        // T* (1 - 1/(rho cv)* dPeos ) * nabla U  [K /s] 
      
        T dPeos; // calculate derivative of Peos wrt. T
        interactionPotential.computeDerivativeConstRho(&dPeos, &localRho, &localTempRatio ) ;

        phi_li += localTemp * nablaU;
        phi_li -= localTempRatio * invVolumetricHeat * dPeos * nablaU;
        
        tPartner->get(iX,iY).template setField<descriptors::PHI_THERMAL>( phi_li );



      // Velocity coupling
        // Needed to transport the AD lattice at same speed as NS lattice
        T *uThermal = tPartner->get(iX,iY).template getFieldPointer<descriptors::VELOCITY>();
        T *uFluid = fluidCell.template getFieldPointer<descriptors::VELOCITY>();
        for (int iD = 0; iD < L::d; ++iD) {
          T u = uField.get(iX-offsetX, iY-offsetY, iD);
          uThermal[iD] =  u;
          uFluid[iD] = u;
        }

      
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void CombinedShanChenSemiHybridThermalCouplingPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}





// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
CombinedShanChenSemiHybridThermalCouplingGenerator2D<T,DESCRIPTOR>::
CombinedShanChenSemiHybridThermalCouplingGenerator2D(
    int x0_, int x1_, int y0_, int y1_,
    T G_, T Tomega_,
    T averageRho_, T gravity_, std::vector<T> dir_,
    FluidProperties<T>& fluidProperties_,
    InteractionPotentialExtended<T,T>& interactionPotential_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    G(G_), Tomega(Tomega_),
    averageRho(averageRho_), gravity(gravity_), dir(dir_),
    fluidProperties(fluidProperties_), interactionPotential(interactionPotential_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* CombinedShanChenSemiHybridThermalCouplingGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject2D* > partners) const
{
  return new CombinedShanChenSemiHybridThermalCouplingPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,
           G, Tomega,
           averageRho, gravity, dir, 
           fluidProperties, interactionPotential, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* CombinedShanChenSemiHybridThermalCouplingGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new CombinedShanChenSemiHybridThermalCouplingGenerator2D<T,DESCRIPTOR>(*this);
}




} // end namespace

#endif
