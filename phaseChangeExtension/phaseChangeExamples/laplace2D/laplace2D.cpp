/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */
 
/* laplace.cpp:
 * This simulation gives a single droplet suspended without body forces
 * periodic boundaries are applied in all three axis

 * The results of the simulations includes multiphase pressure
 * from which the laplace law can be validated
 * multiple droplet diameters are simulated automatically
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
#define NSDESCRIPTOR olb::descriptors::D2Q9< olb::descriptors::VELOCITY,        \
                                             olb::descriptors::FORCE,           \
                                             olb::descriptors::EXTERNAL_FORCE,  \
                                             olb::descriptors::PSI_PSEUDO_RHO>          


// Parameters for the simulation setup

// Resolutions
const int Nx = 125;                   // spacial resolution of the model
const SimType physViscosity = 1e-5;
const SimType physDensity = 1;
// const SimType physDensity = 947.13/8.72518;


const SimType physLength = 100e-6; // m
const SimType maxVelocity = 15;

const SimType lx0   = 5*physLength;    // length of channel
const SimType ly0   = 5*physLength;     // height of channel

const SimType sphereInterfaceWidth = physLength/Nx*8;
const SimType sphereOrigin[2] = {lx0/2., ly0/2. };

const SimType rhoVapourMultiplier = 1.1;
const SimType rhoLiquidMultiplier = 1.0;

const SimType randomFactor = 0.01;

// Pseudopotential multiphase single component dynamics
const std::vector<SimType> rho0 ={1.0};
const SimType G     = -1;


const SimType iterMax = 1000000;

const SimType latticeRelaxationTime = 0.8;
const SimType temperatureRatio = 0.7;     

const SimType sigma = 0.35;                 
const SimType eosA = 0.5/49.0; 

const int MatBase = 0;
const int MatFluid = 1;

// const SimType convP = physDensity * dx * dx / dt / dt;

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
                     SimType rhoVapor, SimType rhoLiquid,
                     SimType sphereRadius )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << endl;

  const SimType omega = converter.getLatticeRelaxationFrequency();

  clout << "omega = " << omega << endl;

  clout << "Define NS Dynamics ..." << endl;
  // Material=0 (BaseMaterial) -->do nothing
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

  olb::SuperVTMwriter2D<SimType> vtmWriter( "laplace2D" );
  
  olb::SuperLatticePhysVelocity2D<SimType, NSDESCRIPTOR> velocity( NSlattice, converter );
  olb::SuperLatticeDensity2D<SimType, NSDESCRIPTOR> density( NSlattice );
  olb::SuperLatticeVelocity2D<SimType, NSDESCRIPTOR> latticeVelocity( NSlattice );
  olb::SuperLatticeField2D<SimType, NSDESCRIPTOR, olb::descriptors::PSI_PSEUDO_RHO> pseudopotential( NSlattice );

  // Multiphase pressure calculations
  olb::SuperCalcPower2D<SimType, SimType> psiSquared (pseudopotential, 2 );
  olb::SuperCalcMultiplication2D<SimType, SimType> halfPsiSquared ( 0.5*G, psiSquared);
  olb::SuperCalcPlus2D<SimType, SimType> addedDensity (density, halfPsiSquared);
  SimType cs2 = 1/olb::descriptors::invCs2<SimType,NSDESCRIPTOR>();
  olb::SuperCalcMultiplication2D<SimType, SimType> pressureNonDim ( cs2 , addedDensity);
  pressureNonDim.getName() = "Lattice pressure";
  SimType PhysPressureConversion = converter.getConversionFactorPressure();
  olb::SuperCalcMultiplication2D<SimType, SimType> pressureMPF ( PhysPressureConversion , pressureNonDim);
  pressureMPF.getName() = "Physical pressure";
  
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( latticeVelocity );
  vtmWriter.addFunctor( pressureNonDim );
  vtmWriter.addFunctor( pressureMPF );
  vtmWriter.addFunctor( density );

  const int  vtkIter  = 50000;
  const int  statIter = 5000;
  const int  ppmIter = 5000;
  const int  plotIter = 250;

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

    SimType maxRho[1];
    int input[3];

    // olb::SuperMin2D<SimType> pressureLatticeMin(pressureNonDim, superGeometry, MatFluid);
    olb::SuperMax2D<SimType> pressureLatticeMax(pressureNonDim, superGeometry, MatFluid);
    // pressureLatticeMin(minRho, input);
    pressureLatticeMax(maxRho, input);

    gplot.setData( SimType (iT), { maxRho[0] } , {"max"}, "left", {'l'}  );
    gplot.writePNG();
  }


  // Writes the vtk files
  if ( iT%vtkIter==0 || forcedSave ) {
    vtmWriter.write( iT );
  }

  // Writes the ppm files
  if ( iT%ppmIter==0 || (iT<=50 && (iT%10==0) ) || forcedSave ) {
    olb::BlockReduction2D2D<SimType> planeReduction( density );
    olb::BlockReduction2D2D<SimType> planeReductionP( pressureNonDim );
    // write output as JPEG
    olb::heatmap::plotParam<SimType> params;
    params.name = "pressure";
    olb::heatmap::write(planeReduction, iT);
    olb::heatmap::write(planeReductionP, iT, params);
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
  olb::SuperVTMwriter2D<SimType> vtmWriter( "laplace2D" );
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


  olb::UnitConverterFromResolutionAndRelaxationTime<SimType, NSDESCRIPTOR> initialConverter(
    (int) Nx,
    (SimType) latticeRelaxationTime,
    (SimType) physLength,          // charPhysLength: reference length of simulation geometry
    (SimType) maxVelocity,         // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (SimType) physViscosity,       // physViscosity: physical kinematic viscosity in __m^2 / s__
    (SimType) physDensity          // physDensity: physical density in __kg / m^3__
  );

  SimType pressureOffset = -initialConverter.getPhysPressure( olb::util::pressureFromDensity<SimType,NSDESCRIPTOR>(0) );
  clout << "Pressure offset: " << pressureOffset << std::endl;

  // Loop through variables
  SimType sphereRadiusFraction;
  for(sphereRadiusFraction = 0.1; sphereRadiusFraction<=0.601; sphereRadiusFraction+=0.1 ){

        
    const SimType sphereRadius = sphereRadiusFraction*lx0/2.0;


    std::ostringstream outputDirStream;
    outputDirStream << std::setprecision(2);
    outputDirStream << "./laplace2D_results/";
    mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    outputDirStream << "r_" << sphereRadiusFraction << "/";
    olb::singleton::directories().setOutputDir( outputDirStream.str() );

    olb::UnitConverterFromResolutionAndRelaxationTime<SimType, NSDESCRIPTOR> converter(
      (int) Nx,
      (SimType) latticeRelaxationTime,
      (SimType) physLength,          // charPhysLength: reference length of simulation geometry
      (SimType) maxVelocity,         // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (SimType) physViscosity,       // physViscosity: physical kinematic viscosity in __m^2 / s__
      (SimType) physDensity,         // physDensity: physical density in __kg / m^3__
      (SimType) pressureOffset       // characteristicPressure: physical pressure at 1 rho
    );

    // Prints the converter log as console output
    converter.print();
    // Writes the converter log in a file
    converter.write("LaplaceConverter");

    // === 3rd Step: Prepare Lattice ===
    olb::SuperLattice2D<SimType, NSDESCRIPTOR> NSlattice(superGeometry);

    // Modified Guo forcing dynamics
    olb::ForcedModifiedGuoBGKdynamics<SimType, NSDESCRIPTOR> ModifiedGuoBulkdynamics (
      converter.getLatticeRelaxationFrequency(),
      olb::instances::getExternalVelocityMomenta<SimType,NSDESCRIPTOR>(), sigma );


    // === Adding coupling to the NS lattice

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
                    ModifiedGuoBulkdynamics,
                    superGeometry,
                    densityLatticeVapour, densityLatticeLiquid,
                    sphereRadius );


    // === 4th Step: Main Loop with Timer ===
    clout << "starting simulation..." << endl;
    olb::util::Timer<SimType> timer( iterMax, superGeometry.getStatistics().getNvoxel() );
    timer.start();

    olb::Gnuplot<SimType> gplot("pressure");

    int interval = 150;        //over the period of X steps
    SimType epsilon = 1e-5;
    olb::util::ValueTracer<SimType> converge( interval, epsilon );

    writeGeometry(NSlattice, superGeometry);
    // getResults(NSlattice, superGeometry, converter, gplot, 0, timer);
    for ( int iT = 0; iT <= iterMax; ++iT ) {

      // === 5th Step: Definition of Initial and Boundary Conditions ===
      // NA

      // === 6th Step: Collide and Stream Execution ===
      NSlattice.collideAndStream();
      NSlattice.communicate();
      NSlattice.executeCoupling();

      // === 7th Step: Computation and Output of the Results ===
      getResults(NSlattice, superGeometry, converter, gplot, iT, timer );

      // == Check for crash or convergence 
      if ( std::isnan( NSlattice.getStatistics().getAverageRho() ) ) { 
          clout << "Terminating this run at step " << iT+1 << " => density is NAN" << std::endl;
          getResults(NSlattice, superGeometry, converter, gplot, iT, timer, true );
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

  } // end for loop 1
} // end main
