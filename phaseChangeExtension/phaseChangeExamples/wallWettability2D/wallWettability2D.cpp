/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/* wallWettability.cpp:
 * This simulation will show the effect of wall density on the contact angle of a droplet
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
#define NSDESCRIPTOR olb::descriptors::D2Q9<olb::descriptors::VELOCITY,olb::descriptors::FORCE,olb::descriptors::EXTERNAL_FORCE,olb::descriptors::PSI_PSEUDO_RHO>

// Parameters for the simulation setup

// Resolutions
const int Nx = 150;                   // spacial resolution of the model
const SimType physViscosity = 1e-5;
const SimType physDensity = 1;

const SimType physLength = 100e-6; // m
const SimType maxVelocity = 40;

const SimType conversionLength = physLength/Nx ;

const SimType rhoVaporMultiplier = 1.0;
const SimType rhoLiquidMultiplier = 0.99;


const SimType lx0   = 10*physLength;      // widths of simulation domain
const SimType ly0   = physLength;               // depth  of simulation domain

const SimType dx = physLength/Nx;

const SimType sphereRadius = ly0/6.;
const SimType sphereInterfaceWidth = sphereRadius/4.;
const SimType sphereOrigin[2] = {lx0/2., sphereRadius+sphereInterfaceWidth};

const int iterMax = 100000;
const int rhoRampIter = 5000;
const int iTsmoothUpdate = 10;

// End constants definition
////

const int materialBase = 0;
const int materialFluid = 1;
const int materialLowerWall = 2;
const int materialUpperBC = 3;

void prepareGeometry( SuperGeometry2D<SimType>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // create lower wall and fluid
  superGeometry.rename( materialBase, materialLowerWall );
  superGeometry.rename( materialLowerWall, materialFluid, 0, 1 );

  // create upper BC
  Vector<SimType,2> extend(lx0, 2*conversionLength);
  Vector<SimType,2> origin(0, ly0 - conversionLength*1.5);

  olb::IndicatorCuboid2D<SimType> upperBC( extend, origin );
  superGeometry.rename( materialLowerWall, materialUpperBC, materialFluid, upperBC );

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
                     olb::BounceBack<SimType, NSDESCRIPTOR>& bounceBackUpperBC,
                     olb::BounceBack<SimType, NSDESCRIPTOR>& bounceBackWall,
                     olb::Dynamics<SimType, NSDESCRIPTOR>& bulkDynamics,
                     SuperGeometry2D<SimType>& superGeometry,
                     SimType rhoVapor, SimType rhoLiquid, SimType rhoWall )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << endl;

  const SimType omega = converter.getLatticeRelaxationFrequency();

  clout << "omega = " << omega << endl;

  clout << "Define NS Dynamics ..." << endl;
  NSlattice.defineDynamics( superGeometry, materialBase, &olb::instances::getNoDynamics<SimType, NSDESCRIPTOR>() );
  NSlattice.defineDynamics( superGeometry, materialFluid, &bulkDynamics );
  NSlattice.defineDynamics( superGeometry, materialUpperBC, &bounceBackUpperBC );

  NSlattice.defineDynamics( superGeometry, materialLowerWall, &bounceBackWall ); 

  clout << "Define Initial Conditions ..." << endl;
  // Fluid Initial conditions
  

  olb::AnalyticalConst2D<SimType,SimType> zeroVelocity( 0,0,0 );

  AnalyticalConst2D<SimType,SimType> noiseVap( rhoVapor*0.1 );
  AnalyticalConst2D<SimType,SimType> noiseLiq( rhoLiquid*0.1 );
  AnalyticalRandom2D<SimType,SimType> random;

  SimType estimatedAdditionalAverageVapour = 2*(rhoWall - rhoVapor)*conversionLength/ly0;
  clout << "estimatedAdditionalAverageVapour: " << estimatedAdditionalAverageVapour << std::endl;
  // this estimate assumes that 2 full fluid nodes equal to the wall density are needed to compensate for the wall density gradient
  // This additional density is then spread over the entire y domain

  olb::AnalyticalConst2D<SimType,SimType> constRhoVapour( (rhoVapor+estimatedAdditionalAverageVapour) );
  olb::AnalyticalConst2D<SimType,SimType> constRhoUpperBC( rhoVapor );
  olb::AnalyticalConst2D<SimType,SimType> constRhoLiquid( rhoLiquid*rhoLiquidMultiplier );
  olb::AnalyticalConst2D<SimType,SimType> constRhoWall( rhoVapor );

  // Density of the Field is:
  // Globally the density of the vapour + some noise
  // Plus at the spehere the difference in densities between vapour and liquid 
  olb::AnalyticalIdentity2D<SimType,SimType> noiseIndicatorVap( random*noiseVap  );
  olb::AnalyticalIdentity2D<SimType,SimType> noiseIndicatorLiq( random*noiseLiq  );
  olb::AnalyticalIdentity2D<SimType,SimType> rhoBase( constRhoVapour );

  olb::SmoothIndicatorCircle2D<SimType,SimType> sphere( sphereOrigin, sphereRadius, sphereInterfaceWidth);
  olb::AnalyticalIdentity2D<SimType,SimType> sphereIndicator( sphere  );

  olb::AnalyticalIdentity2D<SimType,SimType> rho3( rhoBase + noiseIndicatorVap + sphereIndicator*( constRhoLiquid - constRhoVapour + noiseIndicatorLiq) );
  

  //Initialize all values of distribution functions to their local equilibrium

  // Bulk fluid
  NSlattice.defineRhoU( superGeometry, materialFluid, rho3, zeroVelocity );
  NSlattice.iniEquilibrium( superGeometry, materialFluid, rho3, zeroVelocity );

  // Wall
  NSlattice.defineRhoU( superGeometry, materialLowerWall, constRhoWall, zeroVelocity );
  NSlattice.iniEquilibrium( superGeometry, materialLowerWall, constRhoWall, zeroVelocity );


  // upper BC
  NSlattice.defineRhoU( superGeometry, materialUpperBC, constRhoUpperBC, zeroVelocity );
  NSlattice.iniEquilibrium( superGeometry, materialUpperBC, constRhoUpperBC, zeroVelocity );
  
  clout << "Define Psi ..." << endl;
  // Setting initial psi values
  olb::AnalyticalConst2D<SimType,SimType> initialPsi( 1.0 );
  auto fluidIndicator = superGeometry.getMaterialIndicator({materialFluid});
  NSlattice.defineField<olb::descriptors::PSI_PSEUDO_RHO>( fluidIndicator, initialPsi ) ; 

  // Make the lattice ready for simulation
  NSlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Slowly guides wall density to wanted density
void setBoundaryValues( olb::UnitConverter<SimType, NSDESCRIPTOR> &converter,
                        olb::SuperLattice2D<SimType, NSDESCRIPTOR>& NSlattice,
                        int iT,
                        SuperGeometry2D<SimType>& superGeometry,
                        olb::BounceBackVariableRho<SimType, NSDESCRIPTOR> &bounceBackWall,
                        SimType wallRhoMax, SimType rhoVapour )
{
  olb::OstreamManager clout( std::cout,"setBoundaryValues" );

if ( iT%iTsmoothUpdate==0 ) {

    SimType frac[1]= {};
    
    // if below iTmaxStart, get percentage ramp up, else 1
    if (iT<= rhoRampIter){
      // Smooth start curve, polynomial
      olb::PolynomialStartScale<SimType,int> startScale( rhoRampIter, SimType( 1 ) );
      int iTvec[1]= {iT};
      startScale( frac,iTvec );
    }
    else {
      frac[0] = 1;
    }
    bounceBackWall.setDensity( frac[0]*(wallRhoMax-rhoVapour) + rhoVapour );
  }
}

// Output to console and files
void getResults( olb::SuperLattice2D<SimType, NSDESCRIPTOR>& NSlattice,
                 olb::SuperGeometry2D<SimType>& superGeometry,
                 olb::UnitConverter<SimType, NSDESCRIPTOR > &converter,
                 olb::Gnuplot<SimType> gplot,
                 int iT,
                 Timer<SimType>& timer,
                 bool forceSave=false)
{
  olb::OstreamManager clout( std::cout,"getResults" );

  olb::SuperVTMwriter2D<SimType> vtmWriter( "wallWettability2D");
  
  olb::SuperLatticePhysVelocity2D<SimType, NSDESCRIPTOR> velocity( NSlattice, converter );
  olb::SuperLatticeDensity2D<SimType, NSDESCRIPTOR> density( NSlattice );
  olb::SuperLatticeVelocity2D<SimType, NSDESCRIPTOR> latticeVelocity( NSlattice );
  olb::SuperLatticePhysPressure2D<SimType, NSDESCRIPTOR> pressure( NSlattice, converter );

  
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( latticeVelocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( density );

  const int  vtkIter  = 20000;
  const int  heatmapIter = 5000;
  const int  statIter = 5000;

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
  // Removed - not needed

  // Writes the vtk files
  if ( iT%vtkIter==0 || (iT<=50 && (iT%5 == 0) ) || forceSave ) {
    vtmWriter.write( iT );
  }

  // Writes the heatmap files
  if ( iT%heatmapIter==0 || (iT<=50 && (iT%10 == 0) ) || forceSave) {
    olb::BlockReduction2D2D<SimType> planeReduction( density );
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

void writeGeometry( olb::SuperLattice2D<SimType, NSDESCRIPTOR>& NSlattice,
                    olb::SuperGeometry2D<SimType>& superGeometry)
{
  olb::SuperVTMwriter2D<SimType> vtmWriter( "wallWettability" );
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
  cuboidGeometry.setPeriodicity( true, false );

  // Instantiation of a loadBalancer
  olb::HeuristicLoadBalancer<SimType> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  olb::SuperGeometry2D<SimType> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( superGeometry );


  // Loop through variables
  SimType latticeRelaxationTime; 
  SimType temperatureRatio;      
  SimType wallRho;               

  // Removed from loop
  SimType sigma = 0.35;                 
  SimType eosA = 0.5/49.0;                 

  for(latticeRelaxationTime = 0.8; latticeRelaxationTime<=0.801; latticeRelaxationTime+=0.2 ){
    for(temperatureRatio = 0.7; temperatureRatio<=0.701; temperatureRatio+=0.1 ){
      for( wallRho = 0.50; wallRho<=7.51; wallRho+=0.5 ){

          std::ostringstream outputDirStream;
          outputDirStream << std::setprecision(2);
          outputDirStream << "./wallWettabilityResults/";
          mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
          outputDirStream << "tau_" << latticeRelaxationTime << "/";
          mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
          outputDirStream << "tr_" << temperatureRatio << "/";
          mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
          outputDirStream << "Wrho_" << wallRho << "/";
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
          //converter.print();
          // Writes the converter log in a file
          converter.write("wallWettabilityConverter");

          // === 3rd Step: Prepare Lattice ===
          olb::SuperLattice2D<SimType, NSDESCRIPTOR> NSlattice(superGeometry);
          writeGeometry(NSlattice, superGeometry);

          // Modified Guo forcing dynamics
          olb::ForcedModifiedGuoBGKdynamics<SimType, NSDESCRIPTOR> ModifiedGuoBulkdynamics (
            converter.getLatticeRelaxationFrequency(),
            olb::instances::getExternalVelocityMomenta<SimType,NSDESCRIPTOR>(), sigma );


          // Pseudopotential multiphase single component dynamics
          std::vector<SimType> rho0 ={1.0};
          const SimType G     = -1;

          olb::PengRobinsonExt<SimType,SimType> interactionPotential( G, 0.3443, eosA, 2./21., temperatureRatio );
          interactionPotential.write();

          olb::MaxwellConstruction<SimType, decltype(interactionPotential)> maxwellConstruction( interactionPotential, 1e-1, 1e20  );
          SimType rho[2];
          clout << "Performing maxwell construction" << std::endl;
          maxwellConstruction(rho, temperatureRatio); //temperatureRatio
          clout << "Density liquid: " <<  rho[0]  << std::endl;
          clout << "Density vapour: " <<  rho[1]  << std::endl;


          olb::LatticeCouplingGenerator2D<SimType, NSDESCRIPTOR>* multiphaseCoupling;

          multiphaseCoupling = new olb::ModifiedShanChenForcedSingleComponentGenerator2D<SimType, NSDESCRIPTOR>( G, rho0, interactionPotential );


          //auto multiphaseCouplingIndicator = superGeometry.getMaterialIndicator({materialFluid, materialLowerWall });
          NSlattice.addLatticeCoupling( *multiphaseCoupling, NSlattice );


          // A bounce-back node with fictitious density
          //   which is experienced by the partner fluid
          olb::BounceBackVariableRho<SimType, NSDESCRIPTOR> bounceBackUpperBC(rho[1]);

          // A bounce-back node with fictitious density
          //   which is experienced by the partner fluid
          olb::BounceBackVariableRho<SimType, NSDESCRIPTOR> bounceBackWall(rho[1]);
          
          // Prep lattice 
          prepareLattice( converter,
                          NSlattice,
                          bounceBackUpperBC,
                          bounceBackWall,
                          ModifiedGuoBulkdynamics,
                          superGeometry,
                          rho[1], rho[0], wallRho );


          // === 4th Step: Main Loop with Timer ===
          clout << "starting simulation..." << endl;
          olb::util::Timer<SimType> timer( iterMax, superGeometry.getStatistics().getNvoxel() );
          timer.start();
          olb::Gnuplot<SimType> gplot("density");

          // Set up convergence check
          int interval = 100;        //over the period of X steps
          SimType epsilon = 1e-7;
          olb::util::ValueTracer<SimType> converge( interval, epsilon );
          
          // Write inital setup
          getResults(NSlattice, superGeometry, converter, gplot, 0, timer );
          
          for ( int iT = 0; iT <= iterMax; ++iT ) {

            // === 5th Step: Definition of Initial and Boundary Conditions ===
            setBoundaryValues(converter, NSlattice, iT, superGeometry, bounceBackWall, wallRho, rho[1]);

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
              getResults(NSlattice, superGeometry, converter, gplot, iT+1, timer, true );
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

      } // end for 3
    } // end for 2
  } // end for 1
} // end main
