/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/* thermodynamicConsistency2D.cpp:
 * This simulation gives a single droplet suspended without body forces
 * periodic boundaries are applied in all three axis
 * The results of the simulations give the thermodynamic consistency of the
 * forcing method and parameters used. 
 */


#include "olb2D.h"
#include "olb2D.hh"   // include full template code
#include "phaseChangeExtension.h"
#include "phaseChangeExtension.hh"

#include <math.h> 
#include <sys/stat.h>
#include <sys/types.h>

using namespace olb;

typedef double SimType;
#define NSDESCRIPTOR olb::descriptors::D2Q9< olb::descriptors::VELOCITY,         \
                                              olb::descriptors::FORCE,            \
                                              olb::descriptors::EXTERNAL_FORCE,   \
                                              olb::descriptors::PSI_PSEUDO_RHO>

// Parameters for the simulation setup

// Resolutions
const int Nx = 200;                   // spacial resolution of the model
const SimType physViscosity = 1e-5;
const SimType physDensity = 1;

const SimType physLength = 100e-6; // m
const SimType maxVelocity = 30;


const SimType lx0   = physLength;    // length of channel
const SimType ly0   = physLength;     // height of channel


const SimType sphereRadius = physLength/6.;
const SimType sphereInterfaceWidth = sphereRadius/4;
const SimType sphereOrigin[2] = {lx0/2., ly0/2.};

const SimType rhoVapourMultiplier = 1.0;

const SimType iterMax = 5000;


const int MatBase = 0;
const int MatFluid = 1;

// End constants definition
////

void prepareGeometry( SuperGeometry2D<SimType>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;


  superGeometry.rename( MatBase, MatFluid );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}



// Set up the geometry of the simulation
void prepareLattice( olb::UnitConverter<SimType, NSDESCRIPTOR> const& converter,
                     olb::SuperLattice2D<SimType, NSDESCRIPTOR>& NSlattice,
                     olb::Dynamics<SimType, NSDESCRIPTOR>& bulkDynamics,
                     SuperGeometry2D<SimType>& superGeometry,
                     SimType rhoVapor, SimType rhoLiquid )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << endl;

  const SimType omega = converter.getLatticeRelaxationFrequency();
  clout << "omega = " << omega << endl;

  clout << "Define NS Dynamics ..." << endl;
  NSlattice.defineDynamics( superGeometry, MatBase, &olb::instances::getNoDynamics<SimType, NSDESCRIPTOR>() );
  NSlattice.defineDynamics( superGeometry, MatFluid, &bulkDynamics );     
  


  olb::SmoothIndicatorCircle2D<SimType,SimType> sphere( sphereOrigin, sphereRadius, sphereInterfaceWidth);


  clout << "Define Initial Conditions ..." << endl;
  // Fluid Initial conditions
  

  olb::AnalyticalConst2D<SimType,SimType> zeroVelocity( 0,0 );

  AnalyticalConst2D<SimType,SimType> noiseVap( rhoVapor*0.1 );
  AnalyticalConst2D<SimType,SimType> noiseLiq( rhoLiquid*0.1 );
  AnalyticalRandom2D<SimType,SimType> random;

  olb::AnalyticalConst2D<SimType,SimType> constRhoVapour( rhoVapor * rhoVapourMultiplier );
  olb::AnalyticalConst2D<SimType,SimType> constRhoLiquid( rhoLiquid );

  // Density of the Field is:
  // Globally the density of the vapour + some noise
  // Plus at the spehere the difference in densities between vapour and liquid 
  olb::AnalyticalIdentity2D<SimType,SimType> noiseIndicatorVap( random*noiseVap  );
  olb::AnalyticalIdentity2D<SimType,SimType> noiseIndicatorLiq( random*noiseLiq  );
  olb::AnalyticalIdentity2D<SimType,SimType> rhoBase( constRhoVapour );
  olb::AnalyticalIdentity2D<SimType,SimType> sphereIndicator( sphere  );

  olb::AnalyticalIdentity2D<SimType,SimType> rho3( rhoBase + noiseIndicatorVap + sphereIndicator*( constRhoLiquid - constRhoVapour + noiseIndicatorLiq) );
  

  //Initialize all values of distribution functions to their local equilibrium

  // Bulk fluid
  NSlattice.defineRhoU( superGeometry, MatFluid, rho3, zeroVelocity );
  NSlattice.iniEquilibrium( superGeometry, MatFluid, rho3, zeroVelocity );
  

  // Setting initial psi values
  olb::AnalyticalConst2D<SimType,SimType> initialPsi( 1.0 );
  auto fluidIndicator = superGeometry.getMaterialIndicator({MatFluid});
  NSlattice.defineField<olb::descriptors::PSI_PSEUDO_RHO>( fluidIndicator, initialPsi ) ; 

  // Make the lattice ready for simulation
  NSlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Output to console and files
void getResults( olb::SuperLattice2D<SimType, NSDESCRIPTOR>& NSlattice,
                 olb::SuperGeometry2D<SimType>& superGeometry,
                 olb::UnitConverter<SimType, NSDESCRIPTOR > &converter,
                 olb::Gnuplot<SimType> gplot,
                 int iT,
                 Timer<SimType>& timer,
                 bool forcedSave = false)
{
  olb::OstreamManager clout( std::cout,"getResults" );

  olb::SuperVTMwriter2D<SimType> vtmWriter( "thermalConsistency");
  
  olb::SuperLatticePhysVelocity2D<SimType, NSDESCRIPTOR> velocity( NSlattice, converter );
  olb::SuperLatticeDensity2D<SimType, NSDESCRIPTOR> density( NSlattice );
  olb::SuperLatticeVelocity2D<SimType, NSDESCRIPTOR> latticeVelocity( NSlattice );
  olb::SuperLatticePhysPressure2D<SimType, NSDESCRIPTOR> pressure( NSlattice, converter );

  
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( latticeVelocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( density );

  const int  vtkIter  = 25;
  const int  statIter = 1000;
  const int  ppmIter = 25;
  const int  plotIter = 200;

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    olb::SuperLatticeGeometry2D<SimType, NSDESCRIPTOR> geometry( NSlattice, superGeometry );
    olb::SuperLatticeCuboid2D<SimType, NSDESCRIPTOR> cuboid( NSlattice );
    olb::SuperLatticeRank2D<SimType, NSDESCRIPTOR> rank( NSlattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Adds data to GNU plot

  if (iT%plotIter == 0 || forcedSave) {

    SimType minRho[1], maxRho[1];
    int input1[3];
    int input2[3];

    olb::SuperMin2D<SimType> densityLatticeMin(density, superGeometry, MatFluid);
    olb::SuperMax2D<SimType> densityLatticeMax(density, superGeometry, MatFluid);
    densityLatticeMin(minRho, input1);
    densityLatticeMax(maxRho, input2);

    gplot.setData( SimType (iT), { minRho[0], maxRho[0] } , {"min","max"}, "left", {'l','l'}  );
    gplot.writePNG();
  }

  // Writes the ppm files
  if ( iT%vtkIter==0 || forcedSave ) {
    vtmWriter.write( iT );
  }

  if ( iT%ppmIter==0 || (iT<=50 && (iT%10==0) ) || forcedSave ) {
    olb::BlockReduction2D2D<SimType> planeReduction( density );
    // write output as JPEG
    olb::heatmap::write(planeReduction, iT);
  }

  // Writes output on the console
  if ( iT%statIter==0 || forcedSave ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    NSlattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

}

void writeGeometry( olb::SuperLattice2D<SimType, NSDESCRIPTOR>& NSlattice,
                    olb::SuperGeometry2D<SimType>& superGeometry)
{
  olb::SuperVTMwriter2D<SimType> vtmWriter( "thermalConsistency2D" );
  // Writes the geometry, cuboid no. and rank no. as vti file for visualization
  olb::SuperLatticeGeometry2D<SimType, NSDESCRIPTOR> geometry( NSlattice, superGeometry );
  olb::SuperLatticeCuboid2D<SimType, NSDESCRIPTOR> cuboid( NSlattice );
  olb::SuperLatticeRank2D<SimType, NSDESCRIPTOR> rank( NSlattice );
  vtmWriter.write( geometry );
  vtmWriter.write( cuboid );
  vtmWriter.write( rank );
  //vtmWriter.createMasterFile();
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  bool verbosity = false;
  olb::olbInit( &argc, &argv, verbosity );
  olb::OstreamManager clout( std::cout,"main" );

  // === 2nd Step: Prepare Geometry ===
  olb::Vector<SimType,2> extend( lx0, ly0 );
  olb::Vector<SimType,2> origin;
  olb::IndicatorCuboid2D<SimType> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
  #ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = olb::singleton::mpi().getSize();
  #else
  const int noOfCuboids = 1;
  #endif
  const SimType physDeltaX = physLength/Nx;
  olb::CuboidGeometry2D<SimType> cuboidGeometry( cuboid, physDeltaX, noOfCuboids );
  cuboidGeometry.setPeriodicity( true, true );

  // Instantiation of a loadBalancer
  olb::HeuristicLoadBalancer<SimType> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  olb::SuperGeometry2D<SimType> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( superGeometry );


  // Loop through variables
  SimType latticeRelaxationTime; // = 1.0;
  SimType temperatureRatio;      // = 0.7;
  SimType sigma;                 // = 0.31;
  SimType eosA;                 // = 2./49.;

  for(latticeRelaxationTime = 0.8; latticeRelaxationTime<=0.801; latticeRelaxationTime+=0.2 ){
    for(temperatureRatio = 0.6; temperatureRatio<=0.851; temperatureRatio+=0.05 ){
      for( sigma = 0.33; sigma<=0.351; sigma+=0.01 ){
        for( eosA = 0.5/49.; eosA<=2.01/49.; eosA+=0.5/49. ){
          
          std::ostringstream outputDirStream;
          outputDirStream << std::setprecision(2);
          outputDirStream << "./ThermConsistency2D/";
          mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
          outputDirStream << "tau_" << latticeRelaxationTime << "/";
          mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
          outputDirStream << "sigma_" << sigma << "/";
          mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
          outputDirStream << "tr_" << temperatureRatio << "/";
          mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
          outputDirStream << "a_" << eosA << "/";
          olb::singleton::directories().setOutputDir( outputDirStream.str() );

          olb::UnitConverterFromResolutionAndRelaxationTime<SimType, NSDESCRIPTOR> converter(
            (int) Nx,
            (SimType) latticeRelaxationTime,
            (SimType) physLength,          // charPhysLength: reference length of simulation geometry
            (SimType) maxVelocity,         // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (SimType) physViscosity,       // physViscosity: physical kinematic viscosity in __m^2 / s__
            (SimType) physDensity          // physDensity: physical density in __kg / m^3__
          );

          // Prints the converter log as console output
          converter.print();
          // Writes the converter log in a file
          converter.write("ThermConsistencyConverter");

          // === 3rd Step: Prepare Lattice ===
          olb::SuperLattice2D<SimType, NSDESCRIPTOR> NSlattice(superGeometry);

          // Modified Guo forcing dynamics
          olb::Dynamics<SimType, NSDESCRIPTOR>* ModifiedGuoBulkdynamics;

          ModifiedGuoBulkdynamics = new olb::ForcedModifiedGuoBGKdynamics<SimType, NSDESCRIPTOR>(
            converter.getLatticeRelaxationFrequency(),
            olb::instances::getBulkMomenta<SimType,NSDESCRIPTOR>(), sigma);

          // === Adding coupling to the NS lattice
          // Check that SC force is not overwritten by another coupling
          // boussinesq force in thermal coupling by default overwrites SC phase force

          // Pseudopotential multiphase single component dynamics
          std::vector<SimType> rho0 ={1.0};
          const SimType G     = -1;

          olb::PengRobinsonExt<SimType,SimType> interactionPotential( G, 0.3443, eosA, 2./21., temperatureRatio );
          interactionPotential.write();

          olb::MaxwellConstruction<SimType, decltype(interactionPotential)> maxwellConstruction( interactionPotential );
          SimType rho[2];
          clout << "Performing maxwell construction" << std::endl;
          maxwellConstruction(rho, temperatureRatio); //temperatureRatio
          SimType densityLatticeVapour = rho[1]; 
          SimType densityLatticeLiquid = rho[0];
          clout << "Density liquid: " <<  densityLatticeLiquid  << std::endl;
          clout << "Density vapour: " <<  densityLatticeVapour  << std::endl;

          olb::ModifiedShanChenForcedSingleComponentGenerator2D<SimType, NSDESCRIPTOR> multiphaseCoupling( G, rho0, interactionPotential );

          NSlattice.addLatticeCoupling( multiphaseCoupling, NSlattice );

          // Prep lattice 

          prepareLattice( converter,
                          NSlattice,
                          *ModifiedGuoBulkdynamics,
                          superGeometry,
                          densityLatticeVapour, densityLatticeLiquid );


          // === 4th Step: Main Loop with Timer ===
          clout << "starting simulation... (" << outputDirStream.str() << ")" << endl;
          olb::util::Timer<SimType> timer( iterMax, superGeometry.getStatistics().getNvoxel() );
          timer.start();
          olb::Gnuplot<SimType> gplot("density");

          int interval = 75;        //over the period of X steps
          SimType epsilon = 1e-5;
          olb::util::ValueTracer<SimType> converge( interval, epsilon );

          writeGeometry(NSlattice, superGeometry);
          getResults(NSlattice, superGeometry, converter, gplot, 0, timer );
          for ( int iT = 0; iT <= iterMax; ++iT ) {

            // === 5th Step: Definition of Initial and Boundary Conditions ===
            // NA

            // === 6th Step: Collide and Stream Execution ===
            NSlattice.collideAndStream();
            NSlattice.communicate();
            NSlattice.executeCoupling();

            // Psi is initally undef (default: 0), thus first a round of post processing (executeCoupling) needs to be done before collision is able to execute normally
            // However, it is possible to set Psi to i.e. 1.0 and thus allow for first a collision set to occur


            // === 7th Step: Computation and Output of the Results ===
            getResults(NSlattice, superGeometry, converter, gplot, iT+1, timer );
            if ( std::isnan( NSlattice.getStatistics().getAverageRho() ) ) { 
                clout << "Terminating this run at step " << iT+1 << " => density is NAN" << std::endl;
                break;
            }

            converge.takeValue( NSlattice.getStatistics().getAverageEnergy(), false );
            if ( converge.hasConverged() ) {
              clout << "Terminating this run at step " << iT+1 << "=> density has converged" <<std::endl;
              getResults(NSlattice, superGeometry, converter, gplot, iT+1, timer, true );
              break;
            }


          }

          timer.stop();
          timer.printSummary();
          delete ModifiedGuoBulkdynamics;
        } // end for 4
      } // end for 3
    } // end for 2
  } // end for 1
} // end main
