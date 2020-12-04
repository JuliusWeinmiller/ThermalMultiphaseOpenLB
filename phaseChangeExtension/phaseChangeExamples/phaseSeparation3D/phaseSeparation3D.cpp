/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */


#include "olb3D.h"
#include "olb3D.hh"   // include full template code
#include "phaseChangeExtension.h"
#include "phaseChangeExtension.hh"

#include <math.h> 
#include <sys/stat.h>
#include <sys/types.h>

using namespace olb;

typedef double SimType;
#define NSDESCRIPTOR olb::descriptors::D3Q19< olb::descriptors::VELOCITY,         \
                                              olb::descriptors::FORCE,            \
                                              olb::descriptors::EXTERNAL_FORCE,   \
                                              olb::descriptors::PSI_PSEUDO_RHO>

// Parameters for the simulation setup

// Resolutions
const int Nx = 100;                   // spacial resolution of the model
const SimType physViscosity = 1e-5;
const SimType physDensity = 1;

const SimType physLength = 100e-6; // m
const SimType maxVelocity = 10;


const SimType lx0   = physLength;    // length of channel
const SimType ly0   = physLength;     // height of channel
const SimType lz0   = physLength;     // width of channel

const SimType sphereRadius = physLength/6.;
const SimType sphereInterfaceWidth = sphereRadius/4;
const SimType sphereOrigin[3] = {lx0/2., ly0/2., lz0/2};

const SimType iterMax = 10000;

const int MatBase = 0;
const int MatFluid = 1;

// End constants definition
////

void prepareGeometry( SuperGeometry3D<SimType>& superGeometry )
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
                     olb::SuperLattice3D<SimType, NSDESCRIPTOR>& NSlattice,
                     olb::Dynamics<SimType, NSDESCRIPTOR>& bulkDynamics,
                     SuperGeometry3D<SimType>& superGeometry,
                     SimType rhoVapor, SimType rhoLiquid )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << endl;

  const SimType omega = converter.getLatticeRelaxationFrequency();

  clout << "omega = " << omega << endl;

  clout << "Define NS Dynamics ..." << endl;
  // Material=0 (BaseMaterial) -->do nothing
  NSlattice.defineDynamics( superGeometry, MatBase, &olb::instances::getNoDynamics<SimType, NSDESCRIPTOR>() );
  NSlattice.defineDynamics( superGeometry, MatFluid, &bulkDynamics );     
  
  clout << "Define Initial Conditions ..." << endl;
  // Fluid Initial conditions

  olb::AnalyticalConst3D<SimType,SimType> zeroVelocity( 0,0,0 );

  AnalyticalRandom3D<SimType,SimType> random;
  olb::AnalyticalConst3D<SimType,SimType> constRhoAvg( (rhoVapor + rhoLiquid)/10 );

  olb::AnalyticalIdentity3D<SimType,SimType> noiseIndicatorVap( random*constRhoAvg + constRhoAvg  );

  //Initialize all values of distribution functions to their local equilibrium

  // Bulk fluid
  NSlattice.defineRhoU( superGeometry, MatFluid, noiseIndicatorVap, zeroVelocity );
  NSlattice.iniEquilibrium( superGeometry, MatFluid, noiseIndicatorVap, zeroVelocity );
  

  // Setting initial psi values
  // Psi is initally undef (default: 0), thus first a round of post processing (executeCoupling) needs to be done before collision is able to execute normally
  // However, it is possible to set Psi to i.e. 1.0 and thus allow for first a collision set to occur
  olb::AnalyticalConst3D<SimType,SimType> initialPsi( 1.0 );
  auto fluidIndicator = superGeometry.getMaterialIndicator({MatFluid});
  NSlattice.defineField<olb::descriptors::PSI_PSEUDO_RHO>( fluidIndicator, initialPsi ) ; 

  // Make the lattice ready for simulation
  NSlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Output to console and files
void getResults( olb::SuperLattice3D<SimType, NSDESCRIPTOR>& NSlattice,
                 olb::SuperGeometry3D<SimType>& superGeometry,
                 olb::UnitConverter<SimType, NSDESCRIPTOR > &converter,
                 olb::Gnuplot<SimType> gplot,
                 int iT,
                 Timer<SimType>& timer)
{
  olb::OstreamManager clout( std::cout,"getResults" );

  olb::SuperVTMwriter3D<SimType> vtmWriter( "phaseSep");
  
  olb::SuperLatticePhysVelocity3D<SimType, NSDESCRIPTOR> velocity( NSlattice, converter );
  olb::SuperLatticeDensity3D<SimType, NSDESCRIPTOR> density( NSlattice );
  olb::SuperLatticeVelocity3D<SimType, NSDESCRIPTOR> latticeVelocity( NSlattice );
  olb::SuperLatticePhysPressure3D<SimType, NSDESCRIPTOR> pressure( NSlattice, converter );

  
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( latticeVelocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( density );

  const int  vtkIter  = 200;
  const int  statIter = 100;
  const int  plotIter = 50;

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    olb::SuperLatticeGeometry3D<SimType, NSDESCRIPTOR> geometry( NSlattice, superGeometry );
    olb::SuperLatticeCuboid3D<SimType, NSDESCRIPTOR> cuboid( NSlattice );
    olb::SuperLatticeRank3D<SimType, NSDESCRIPTOR> rank( NSlattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Adds data to GNU plot

  if (iT%plotIter == 0) {

    SimType minRho[1], maxRho[1];
    int input[3];

    olb::SuperMin3D<SimType, SimType> densityLatticeMin(density, superGeometry, MatFluid);
    olb::SuperMax3D<SimType, SimType> densityLatticeMax(density, superGeometry, MatFluid);
    densityLatticeMin(minRho, input);
    densityLatticeMax(maxRho, input);

    gplot.setData( SimType (iT), { minRho[0], maxRho[0] } , {"min","max"}, "left", {'l','l'}  );
    gplot.writePNG();
  }

  // Writes the ppm files
  if ( iT%vtkIter==0 ) {
    vtmWriter.write( iT );

    olb::BlockReduction3D2D<SimType> planeReduction( density, {0, 0, 1} );
    // write output as JPEG
    olb::heatmap::write(planeReduction, iT);
  }

  // Writes output on the console
  if ( iT%statIter==0 && iT>=0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    NSlattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

}

// Only write geomerty files
// Handy for debugging geometry related issues
void writeGeometry( olb::SuperLattice3D<SimType, NSDESCRIPTOR>& NSlattice,
                    olb::SuperGeometry3D<SimType>& superGeometry)
{
  olb::SuperVTMwriter3D<SimType> vtmWriter( "phaseSep" );
  // Writes the geometry, cuboid no. and rank no. as vti file for visualization
  olb::SuperLatticeGeometry3D<SimType, NSDESCRIPTOR> geometry( NSlattice, superGeometry );
  olb::SuperLatticeCuboid3D<SimType, NSDESCRIPTOR> cuboid( NSlattice );
  olb::SuperLatticeRank3D<SimType, NSDESCRIPTOR> rank( NSlattice );
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
  olb::singleton::directories().setOutputDir( "./tmp/" );
  olb::OstreamManager clout( std::cout,"main" );

  // === 2nd Step: Prepare Geometry ===
  olb::Vector<SimType,3> extend( lx0, ly0, lz0 );
  olb::Vector<SimType,3> origin;
  olb::IndicatorCuboid3D<SimType> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
  #ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = olb::singleton::mpi().getSize();
  #else
  const int noOfCuboids = 1;
  #endif
  const SimType physDeltaX = physLength/Nx;
  olb::CuboidGeometry3D<SimType> cuboidGeometry( cuboid, physDeltaX, noOfCuboids );
  cuboidGeometry.setPeriodicity( true, true, true );

  // Instantiation of a loadBalancer
  olb::HeuristicLoadBalancer<SimType> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  olb::SuperGeometry3D<SimType> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( superGeometry );


  // Loop through variables
  SimType latticeRelaxationTime = 0.8; // = 1.0;
  SimType temperatureRatio = 0.7;      // = 0.7;
  SimType sigma = 0.31;                 // = 0.31;
  SimType eosA = 0.5/49;                 // = 2./49.;

  olb::UnitConverterFromResolutionAndRelaxationTime<SimType, NSDESCRIPTOR> converter(
    (int) Nx,
    (SimType) latticeRelaxationTime,
    (SimType) physLength,          // charPhysLength: reference length of simulation geometry
    (SimType) maxVelocity,         // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (SimType) physViscosity,       // physViscosity: physical kinematic viscosity in __m^2 / s__
    (SimType) physDensity          // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  //converter.print();
  // Writes the converter log in a file
  converter.write("PhaseSepConverter");

  // === 3rd Step: Prepare Lattice ===
  olb::SuperLattice3D<SimType, NSDESCRIPTOR> NSlattice(superGeometry);

  // Modified Guo forcing dynamics

  olb::ForcedModifiedGuoBGKdynamics<SimType, NSDESCRIPTOR> ModifiedGuoBulkdynamics (
    converter.getLatticeRelaxationFrequency(),
    olb::instances::getExternalVelocityMomenta<SimType,NSDESCRIPTOR>(), sigma );


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

  olb::ModifiedShanChenForcedSingleComponentGenerator3D<SimType, NSDESCRIPTOR> multiphaseCoupling( G, rho0, interactionPotential );


  NSlattice.addLatticeCoupling( multiphaseCoupling, NSlattice );

  // Prep lattice 

  prepareLattice( converter,
                  NSlattice,
                  ModifiedGuoBulkdynamics,
                  superGeometry,
                  densityLatticeVapour, densityLatticeLiquid );


  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  olb::util::Timer<SimType> timer( iterMax, superGeometry.getStatistics().getNvoxel() );
  timer.start();
  olb::Gnuplot<SimType> gplot("density");

  writeGeometry(NSlattice, superGeometry);
  getResults(NSlattice, superGeometry, converter, gplot, 0, timer );
  for ( int iT = 0; iT <= iterMax; ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // NA

    // === 6th Step: Collide and Stream Execution ===
    NSlattice.collideAndStream();
    NSlattice.communicate();
    NSlattice.executeCoupling();

    // === 7th Step: Computation and Output of the Results ===
    getResults(NSlattice, superGeometry, converter, gplot, iT+1, timer );
    if ( std::isnan( NSlattice.getStatistics().getAverageRho() ) ) { 
        clout << "Terminating this run at step " << iT+1 << " => density is NAN" << std::endl;
        break;
    }

  }

  timer.stop();
  timer.printSummary();


} // end main
