
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/* poolBoiling2D.cpp
 * This simulation showcases a implementation of gravity
 * in addition to phase changing multiphase flow
 * This method has some issues with the metastability of
 * the non-ideal equation of state. This can be alleviated
 * using high disturbance fluctuations. Another factor
 * that helps in pool boiling is using a higher gravity.
 */


#include "olb2D.h"
#include "olb2D.hh"   // include full template code
#include "phaseChangeExtension.h"
#include "phaseChangeExtension.hh"

using namespace olb;
#define ThermalMRT

typedef double SimType;
#define NSDESCRIPTOR olb::descriptors::D2Q9<olb::descriptors::VELOCITY,        \
                                             olb::descriptors::FORCE,           \
                                             olb::descriptors::EXTERNAL_FORCE,  \
                                             olb::descriptors::PSI_PSEUDO_RHO>          

#ifdef ThermalMRT
#define TDESCRIPTOR olb::descriptors::D2Q5<olb::descriptors::tag::MRT,      \
                                           olb::descriptors::VELOCITY,      \
                                           olb::descriptors::PHI_THERMAL,   \
                                           olb::descriptors::PREV_T_V,      \
                                           olb::descriptors::PREV_PHI>
#else
#define TDESCRIPTOR olb::descriptors::D2Q5<olb::descriptors::VELOCITY,      \
                                           olb::descriptors::PHI_THERMAL,   \
                                           olb::descriptors::PREV_T_V,      \
                                           olb::descriptors::PREV_PHI>
#endif

// Parameters for the simulation setup

// Resolutions
const int Nx = 75;            // spacial resolution of the model
const SimType Nt = 7;         // temporal resolution of the model 

const SimType physLength = 100e-6; // m
const SimType lx0   = 8*physLength;     // length of channel
const SimType ly0   = 3*physLength;     // height of channel
const SimType conversionLength = physLength/Nx ;

const SimType physViscosity = 1e-5; //  m2 / s

const SimType maxVelocity = 20;
const SimType charTime = physLength/maxVelocity;


const SimType endTimeInit = 10.0*charTime;        // Time to allow fluid to come to rest
const SimType maxPhysT = 1000.0*charTime;         // max. simulation time in s, SI unit
const SimType maxPhysT_ms = maxPhysT*1000;        // max. simulation time in ms

const SimType dx = physLength/ Nx;
const SimType dt = charTime/Nx/Nt; // L/V / Nx/Nt

// Tau from visc
const SimType invCs2  = olb::descriptors::invCs2<SimType,NSDESCRIPTOR>() ;
const SimType tauCalc  = physViscosity/dx/dx*dt * invCs2 + 0.5;

const SimType noiseFraction = 0.05;
const SimType vapourDensityMultiplier = 1.05;

const SimType wallRho = 5.0;

const SimType initLiquidHeightFraction = 0.4;
const SimType liquidInterfaceWidth = physLength/Nx*8;  // Y nodes across

const SimType gravityFactor = 10000;
const SimType grav = -9.81*gravityFactor/dx*dt*dt;

// == Define simulations constants
const SimType sigma = 0.33;
const SimType eosA  = 1.5/49.0; 
std::vector<SimType> rho0 ={1.0};
const SimType G     = -1;


////
// Thermal constants
const SimType TempCrit = 647; // [K] critical temperature of water in Kelvin

const SimType initTempRatio = 0.70;
const SimType newTempRatio  = 0.90;
const SimType charLowTemp = 0.70*TempCrit;
const SimType charHighTemp = 0.90*TempCrit;

const SimType tempNoiseFraction = 2.5;


const SimType specificGasConstantVapour = 461.52;  
const SimType propertyTemperatureRatio = 0.6;
const SimType physWaterTemp = TempCrit * propertyTemperatureRatio - 273.15;
const SimType physDensityLiquid = 947.13;
const SimType specificHeatCapacityVolumeLiquid = 3.8e3;    // J / kg K REAL
const SimType thermalConductivityLiquid = 682e-3; // W / m K REAL
const SimType physDensityVapour = 0.9643;
const SimType specificHeatCapacityVolumeVapour = 2.3e3 - specificGasConstantVapour; // Cv = Cp - Rsp REAL
const SimType thermalConductivityVapour =  31e-3;         // W / m K  REAL        

const SimType liquidDensityRatio = physDensityLiquid/8.72518;


const SimType specificHeatCapacityVolume = 3e3;    // J / kg K
const SimType thermalConductivity = 1; // W / m K       

// End constants definition
////

const int materialBase = 0;
const int materialFluid = 1;
const int materialLowerWall = 2;
const int materialUpperWall = 3;

void prepareGeometry( SuperGeometry2D<SimType>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // create lower wall and fluid
  superGeometry.rename( materialBase, materialLowerWall );
  superGeometry.rename( materialLowerWall, materialFluid, 0, 1 );

  // create upper BC
  Vector<SimType,2> extend(lx0, ly0);
  Vector<SimType,2> origin(0, ly0 - conversionLength*1.5);

  olb::IndicatorCuboid2D<SimType> upperBC( extend, origin );
  superGeometry.rename( materialLowerWall, materialUpperWall, materialFluid, upperBC );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}



// Set up the geometry of the simulation
void prepareLattice( olb::ThermalUnitConverter<SimType, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     olb::SuperLattice2D<SimType, NSDESCRIPTOR>& NSlattice,
                     olb::SuperLattice2D<SimType, TDESCRIPTOR>& ADlattice,
                     olb::Dynamics<SimType, NSDESCRIPTOR>& bulkDynamics,
                     olb::Dynamics<SimType, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     olb::BounceBack<SimType, NSDESCRIPTOR>& bounceBackUpperWall,
                     olb::BounceBack<SimType, NSDESCRIPTOR>& bounceBackLowerWall,
                     olb::sOnLatticeBoundaryCondition2D<SimType,TDESCRIPTOR>& TboundaryCondition,
                     olb::sOnLatticeBoundaryCondition2D<SimType,NSDESCRIPTOR>& NSboundaryCondition,
                     olb::SuperGeometry2D<SimType>& superGeometry,
                     SimType TempInitial, SimType rhoVapor, SimType rhoLiquid)
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << endl;

  const SimType omega = converter.getLatticeRelaxationFrequency();
  const SimType Tomega  = converter.getLatticeThermalRelaxationFrequency();

  clout << "NS omega = " << omega << endl;
  clout << "T  omega = " << Tomega << endl;

  auto fluidIndicator = superGeometry.getMaterialIndicator( {materialFluid});
  auto bulkIndicator = superGeometry.getMaterialIndicator({materialFluid, materialLowerWall, materialUpperWall});
  auto wallIndicator = superGeometry.getMaterialIndicator({materialLowerWall, materialUpperWall});


  clout << "Define NS Dynamics ..." << endl;
  NSlattice.defineDynamics( superGeometry, materialBase, &olb::instances::getNoDynamics<SimType, NSDESCRIPTOR>() );
  NSlattice.defineDynamics( fluidIndicator, &bulkDynamics );
  NSlattice.defineDynamics( superGeometry, materialUpperWall, &bounceBackUpperWall );
  NSlattice.defineDynamics( superGeometry, materialLowerWall, &bounceBackLowerWall );

  clout << "Define AD Dynamics ..." << endl;
  ADlattice.defineDynamics( superGeometry, materialBase, &olb::instances::getNoDynamics<SimType, TDESCRIPTOR>());

  ADlattice.defineDynamics( fluidIndicator, &advectionDiffusionBulkDynamics );
  ADlattice.defineDynamics( wallIndicator,  &advectionDiffusionBulkDynamics );
  TboundaryCondition.addTemperatureBoundary(superGeometry, materialUpperWall, Tomega);
  TboundaryCondition.addTemperatureBoundary(superGeometry, materialLowerWall, Tomega);


  clout << "Define Initial Conditions ..." << endl;
  // Fluid Initial conditions
  olb::AnalyticalConst2D<SimType,SimType> zeroVelocity( 0,0 );

  // random density around average density
  olb::AnalyticalConst2D<SimType,SimType> noiseVap( rhoVapor*noiseFraction );
  olb::AnalyticalConst2D<SimType,SimType> noiseLiq( rhoLiquid*noiseFraction );
  olb::AnalyticalConst2D<SimType,SimType> constRhoVapour( rhoVapor*vapourDensityMultiplier );
  olb::AnalyticalConst2D<SimType,SimType> constRhoLiquid( rhoLiquid );
  olb::AnalyticalConst2D<SimType,SimType> constRhoWall( wallRho );
  olb::AnalyticalRandom2D<SimType,SimType> random; // between 0 and 1

  SimType offsetEstimate = 1 - (wallRho - rhoVapor)/(rhoLiquid-rhoVapor) ;

  SimType cubeOrigin[2] = {lx0/2., (ly0*(initLiquidHeightFraction)/2.0 + liquidInterfaceWidth * offsetEstimate ) };
  olb::SmoothIndicatorCuboid2D<SimType,SimType> cube(cubeOrigin, 2*lx0, ly0*initLiquidHeightFraction, liquidInterfaceWidth);

  olb::AnalyticalIdentity2D<SimType,SimType> fluidDensityLower( cube );
  olb::AnalyticalIdentity2D<SimType,SimType> fluidDensityLowerRandom( constRhoVapour + noiseVap + cube*( constRhoLiquid - constRhoVapour + noiseLiq) );


  // Temperature Inital Conditions
  olb::AnalyticalConst2D<SimType,SimType> T_init(  converter.getLatticeTemperature(TempInitial));
  olb::AnalyticalIdentity2D<SimType,SimType> T_0( T_init );

  //Initialize all values of distribution functions to their local equilibrium

  // Setting initial psi values
  olb::AnalyticalConst2D<SimType,SimType> initialPsi( 1.0 );
  NSlattice.defineField<olb::descriptors::PSI_PSEUDO_RHO>( superGeometry, materialFluid, initialPsi ) ; 

  // Fluid init
  NSlattice.defineRhoU( superGeometry, materialFluid, fluidDensityLowerRandom, zeroVelocity );
  NSlattice.iniEquilibrium( superGeometry, materialFluid, fluidDensityLowerRandom, zeroVelocity );

  NSlattice.defineRhoU( superGeometry, materialUpperWall, constRhoVapour, zeroVelocity );
  NSlattice.iniEquilibrium( superGeometry, materialUpperWall, constRhoVapour, zeroVelocity );

  NSlattice.defineRhoU( superGeometry, materialLowerWall, constRhoWall, zeroVelocity );
  NSlattice.iniEquilibrium( superGeometry, materialLowerWall, constRhoWall, zeroVelocity );

  // Fluid walls do not need any init

  // Thermal init
  ADlattice.defineRho(wallIndicator, T_0);
  ADlattice.iniEquilibrium(wallIndicator, T_0, zeroVelocity);
  ADlattice.defineRho(fluidIndicator, T_0);
  ADlattice.iniEquilibrium(fluidIndicator, T_0, zeroVelocity);


  // Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( olb::ThermalUnitConverter<SimType, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                        olb::SuperLattice2D<SimType, NSDESCRIPTOR>& NSlattice,
                        olb::SuperLattice2D<SimType, TDESCRIPTOR>& ADlattice,
                        int iT,
                        SuperGeometry2D<SimType>& superGeometry )
{

  olb::OstreamManager clout( std::cout,"setBoundaryValues" );

}

// Output to console and files
void getResults( olb::SuperLattice2D<SimType, NSDESCRIPTOR>& NSlattice,
                 olb::SuperLattice2D<SimType, TDESCRIPTOR>& ADlattice,
                 olb::SuperGeometry2D<SimType>& superGeometry,
                 olb::ThermalUnitConverter<SimType, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                 int iT,
                 Timer<SimType>& timer,
                 bool forcedSave = false)
{
  olb::OstreamManager clout( std::cout,"getResults" );

  olb::SuperVTMwriter2D<SimType> vtmWriterNS( "fluidLatticeResults" );
  olb::SuperLatticePhysVelocity2D<SimType, NSDESCRIPTOR> velocity( NSlattice, converter );
  olb::SuperLatticeDensity2D<SimType, NSDESCRIPTOR> density( NSlattice );
  olb::SuperLatticeVelocity2D<SimType, NSDESCRIPTOR> latticeVelocity( NSlattice );
  olb::SuperLatticePhysPressure2D<SimType, NSDESCRIPTOR> pressure( NSlattice, converter );

  olb::SuperVTMwriter2D<SimType> vtmWriterThermal( "thermalLatticeResults" );
  olb::SuperLatticeDensity2D<SimType, TDESCRIPTOR> latticeTemperature( ADlattice );
  olb::SuperLatticePhysTemperature2D<SimType, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
  
  vtmWriterNS.addFunctor( velocity );
  vtmWriterNS.addFunctor( latticeVelocity );
  vtmWriterNS.addFunctor( pressure );
  vtmWriterNS.addFunctor( density );
  vtmWriterThermal.addFunctor( latticeTemperature );
  vtmWriterThermal.addFunctor( temperature );

  const int vtkIter  = converter.getLatticeTime( 20000*charTime/Nx/Nt );
  const int statIter = converter.getLatticeTime( 10000*charTime/Nx/Nt );
  const int ppmIter = converter.getLatticeTime( 5000*charTime/Nx/Nt ); 

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    olb::SuperLatticeGeometry2D<SimType, NSDESCRIPTOR> geometry( NSlattice, superGeometry );
    olb::SuperLatticeCuboid2D<SimType, NSDESCRIPTOR> cuboid( NSlattice );
    olb::SuperLatticeRank2D<SimType, NSDESCRIPTOR> rank( NSlattice );
    vtmWriterNS.write( geometry );
    vtmWriterNS.write( cuboid );
    vtmWriterNS.write( rank );
    vtmWriterNS.createMasterFile();
    vtmWriterThermal.createMasterFile();
  }

  // Writes the vtk files
  if ( iT%vtkIter==0 || (iT<=100 && (iT%10==0) ) || forcedSave ) {
    vtmWriterNS.write( iT );
    vtmWriterThermal.write( iT );
  }

  // Writes the ppm files
  if ( iT%ppmIter==0 || (iT<=50 && (iT%5==0) ) || forcedSave ) {
    olb::BlockReduction2D2D<SimType> planeReductionRho( density );
    olb::BlockReduction2D2D<SimType> planeReductionT( temperature );
    planeReductionT.getName() = "Phys Temperature";
    // write output as JPEG
    olb::heatmap::plotParam<SimType> thermalPlotParam;
    thermalPlotParam.name = "physTemperature";
    thermalPlotParam.minValue = initTempRatio*TempCrit-10;
    thermalPlotParam.maxValue = newTempRatio*TempCrit+10;
    olb::heatmap::write(planeReductionRho, iT);
    olb::heatmap::write(planeReductionT, iT, thermalPlotParam);
  }

  // Writes output on the console
 if ( (iT%statIter==0 && iT>=0) || forcedSave ) {
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
  olb::SuperVTMwriter2D<SimType> vtmWriter( "poolBoiling" );
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
  olb::singleton::directories().setOutputDir( "./tmp/" );
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
  olb::CuboidGeometry2D<SimType> cuboidGeometry( cuboid, conversionLength, noOfCuboids );
  cuboidGeometry.setPeriodicity( true, false );

  // Instantiation of a loadBalancer
  olb::HeuristicLoadBalancer<SimType> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  olb::SuperGeometry2D<SimType> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( superGeometry );

// === 3rd Step: Prepare Lattice ===
  olb::SuperLattice2D<SimType, TDESCRIPTOR> ADlattice(superGeometry);
  olb::SuperLattice2D<SimType, NSDESCRIPTOR> NSlattice(superGeometry);
  writeGeometry(NSlattice, superGeometry);
  
  // Pseudopotential multiphase single component dynamics
  SimType molarMass = 18.0152;
  SimType Rdimless = 8314.0/molarMass /(dx*dx/dt/dt) * TempCrit;
 
  olb::PengRobinsonExt<SimType,SimType> interactionPotential( G, 0.3443, eosA, 2./21., initTempRatio, Rdimless );
  interactionPotential.write();

  // Density from maxwell construction
  SimType rho[2];
  olb::MaxwellConstruction<SimType, decltype(interactionPotential)> maxwellConstruction( interactionPotential );
  maxwellConstruction(rho, initTempRatio); //temperatureRatio
  SimType densityLatticeVapour = rho[1]; 
  SimType densityLatticeLiquid = rho[0];

  olb::ThermalUnitConverter<SimType, NSDESCRIPTOR, TDESCRIPTOR> converter(
    (SimType)   dx,                // physDeltaX: spacing between two lattice cells in __m__
    (SimType)   dt,                // physDeltaT = charLatticeVelocity / charPhysVelocity * dx
    (SimType)   physLength,        // charPhysLength: reference length of simulation geometry
    (SimType)   maxVelocity,       // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (SimType)   physViscosity,     // physViscosity: physical kinematic viscosity in __m^2 / s__
    (SimType)   liquidDensityRatio,               // physDensity: physical density in __kg / m^3__
    (SimType)   thermalConductivity,          // physThermalConductivity
    (SimType)   specificHeatCapacityVolume, // physSpecificHeatCapacity
    (SimType)   0.0,                    //physThermalExpansionCoefficient
    (SimType)   charLowTemp,
    (SimType)   charHighTemp
  );


  // Prints the converter log as console output
  converter.print();
  converter.write("thermalUnitConverter");
  clout << "Density liquid: " <<  densityLatticeLiquid  << std::endl;
  clout << "Density vapour: " <<  densityLatticeVapour  << std::endl;

  olb::FluidProperties<SimType> fluidProperties( 
                    converter.getLatticeDensity(physDensityVapour),
                    converter.getLatticeDensity(physDensityLiquid),
                    converter.getLatticeSpecificHeatCapacity( specificHeatCapacityVolumeVapour ),
                    converter.getLatticeSpecificHeatCapacity( specificHeatCapacityVolumeLiquid ),
                    converter.getLatticeThermalConductivity( thermalConductivityVapour ),
                    converter.getLatticeThermalConductivity( thermalConductivityLiquid ),
                    converter.getCharPhysLowTemperature(), converter.getCharPhysHighTemperature(), TempCrit ); 

  fluidProperties.print();
  fluidProperties.write();


  // Modified Guo forcing dynamics
  olb::ForcedModifiedGuoBGKdynamics<SimType, NSDESCRIPTOR> fluidSimulationDynamics (
   converter.getLatticeRelaxationFrequency(),
   olb::instances::getBulkMomenta<SimType,NSDESCRIPTOR>(), sigma 
   );

  olb::Dynamics<SimType, TDESCRIPTOR>* thermalSimulationDynamics;
  // Modified Advection Diffusion dynamics
  #ifdef ThermalMRT
  thermalSimulationDynamics = new olb::PseudopotentialAdvectionDiffusionMRTdynamics<SimType, TDESCRIPTOR>(
    converter.getLatticeThermalRelaxationFrequency(),
    olb::instances::getBulkMomenta<SimType,TDESCRIPTOR>()
  );
  #else
  thermalSimulationDynamics = new olb::PseudopotentialAdvectionDiffusionBGKdynamics<SimType, TDESCRIPTOR>(
    converter.getLatticeThermalRelaxationFrequency(),
    olb::instances::getBulkMomenta<SimType,TDESCRIPTOR>()
  );
  #endif

  // boundary condition
  olb::sOnLatticeBoundaryCondition2D<SimType,TDESCRIPTOR> TboundaryCondition(ADlattice);
  olb::createAdvectionDiffusionBoundaryCondition2D<SimType,TDESCRIPTOR>(TboundaryCondition);  // OnLatticeAdvectionDiffusionBoundaryCondition3D

  // A bounce-back node with fictitious density
  //   which is experienced by the partner fluid
  olb::BounceBack<SimType, NSDESCRIPTOR> bounceBackUpperWall(densityLatticeVapour);
  olb::BounceBack<SimType, NSDESCRIPTOR> bounceBackLowerWall(wallRho);

  olb::sOnLatticeBoundaryCondition2D<SimType,NSDESCRIPTOR> NSboundaryCondition(NSlattice);


  prepareLattice( converter,
                  NSlattice, ADlattice,
                  fluidSimulationDynamics, *thermalSimulationDynamics,
                  bounceBackUpperWall, bounceBackLowerWall,
                  TboundaryCondition, NSboundaryCondition,
                  superGeometry,
                  initTempRatio*TempCrit,
                  densityLatticeVapour, densityLatticeLiquid);


  SimType gravity =  -9.81  / converter.getConversionFactorLength() * converter.getConversionFactorTime() * converter.getConversionFactorTime() ;
  std::vector<SimType> dir{0.0, 1.0};
  SimType averageRho = (initLiquidHeightFraction*densityLatticeLiquid + (1-initLiquidHeightFraction)*densityLatticeVapour*vapourDensityMultiplier)*(1+noiseFraction);
  clout << "Average Density estimate: " <<  averageRho  << std::endl;


    // === 4th Step: Prepare Lattice Coupling
    olb::LatticeCouplingGenerator2D<SimType, NSDESCRIPTOR> *thermalMultiphaseCoupling;
    #ifdef ThermalMRT
      thermalMultiphaseCoupling = new olb::CombinedShanChenThermalMRTCouplingGenerator2D<SimType, NSDESCRIPTOR>(
          0, converter.getLatticeLength(lx0), 
          0, converter.getLatticeLength(ly0),
          G, converter.getLatticeThermalRelaxationFrequency(),
          averageRho, gravity, {0,1},
          fluidProperties, interactionPotential);
    #else
      thermalMultiphaseCoupling = new olb::CombinedShanChenThermalCouplingGenerator2D<SimType, NSDESCRIPTOR>(
          0, converter.getLatticeLength(lx0), 
          0, converter.getLatticeLength(ly0),
          G, converter.getLatticeThermalRelaxationFrequency(),
          averageRho, gravity, {0,1},
          fluidProperties, interactionPotential);
    #endif

  // auto fluidIndicator = superGeometry.getMaterialIndicator({materialFluid}); 
  
  NSlattice.addLatticeCoupling( *thermalMultiphaseCoupling, ADlattice);

// === 5th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  olb::util::Timer<SimType> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();
  int latticeInitTime = converter.getLatticeTime( endTimeInit );

  // Perform the cycle once to get correct NS.average rho
  getResults(NSlattice, ADlattice, superGeometry, converter, 0, timer );

  olb::AnalyticalConst2D<SimType,SimType> zeroVelocity( 0,0 );
  olb::AnalyticalConst2D<SimType,SimType> T_init( converter.getLatticeTemperature(TempCrit * initTempRatio) );
  
  for ( int iT = 0; iT < latticeInitTime; ++iT ){

    // === 6th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, NSlattice, ADlattice, iT, superGeometry);

    // === 7th Step: Collide and Stream Execution ===
  
    NSlattice.collideAndStream();
    ADlattice.collideAndStream();
    NSlattice.communicate();
    ADlattice.communicate();
    NSlattice.executeCoupling();

    // === Reset temperature field
    ADlattice.defineRho(superGeometry, materialFluid, T_init);
    ADlattice.iniEquilibrium(superGeometry, materialFluid, T_init, zeroVelocity);

    // === 8th Step: Computation and Output of the Results ===
    getResults(NSlattice, ADlattice, superGeometry, converter, iT+1, timer );
    if (iT%5==0){ //check every 5 steps
      if ( std::isnan( NSlattice.getStatistics().getAverageRho() ) ) { 
          clout << "Terminating this run at step " << iT+1 << " => density is NAN" << std::endl;
          getResults(NSlattice, ADlattice, superGeometry, converter, iT+1, timer, true );
          break;
      }
    }
  }
  
  olb::AnalyticalConst2D<SimType,SimType> T_new(   converter.getLatticeTemperature(TempCrit * newTempRatio));//
  olb::AnalyticalRandom2D<SimType,SimType> random; // between 0 and 1
  olb::AnalyticalConst2D<SimType,SimType> randomMagnitude (2*tempNoiseFraction) ;
  olb::AnalyticalConst2D<SimType,SimType> half(0.5);
  olb::AnalyticalIdentity2D<SimType,SimType> T_0(  T_new + randomMagnitude*(random*random-half*half)  );

  ADlattice.defineRho(superGeometry, materialLowerWall, T_0);

  olb::SuperLatticePhysTemperature2D<SimType, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
  olb::SuperAverage2D<SimType, SimType> localAverage(temperature, superGeometry, materialLowerWall );

  SimType rhoAvgOut[1]; 
  int input[3];
  localAverage(rhoAvgOut, input);

  clout << "Avg wall temp: " <<  rhoAvgOut[0] << std::endl;
  clout << "New wall temp: " <<  TempCrit * newTempRatio << std::endl;


  for ( int iT = latticeInitTime; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {

    // === 6th Step: Definition of Initial and Boundary Conditions ===
    //setBoundaryValues( converter, NSlattice, ADlattice, iT - latticeInitTime, superGeometry, newTempRatio);

    // === 7th Step: Collide and Stream Execution ===
  
    NSlattice.collideAndStream();
    ADlattice.collideAndStream();
    NSlattice.communicate();
    ADlattice.communicate();
    NSlattice.executeCoupling();

    // === 8th Step: Computation and Output of the Results ===
    getResults(NSlattice, ADlattice, superGeometry, converter, iT+1, timer );
    if (iT%50==0){ //check every 50 steps
      if ( std::isnan( NSlattice.getStatistics().getAverageRho() ) ) { 
          clout << "Terminating this run at step " << iT+1 << " => density is NAN" << std::endl;
          getResults(NSlattice, ADlattice, superGeometry, converter, iT+1, timer, true );
          break;
      }

      // olb::AnalyticalIdentity2D<SimType,SimType> T_0(  T_new + randomMagnitude*(random-half)  );
      // ADlattice.defineRho(superGeometry, materialLowerWall, T_0);
    }
  }

  delete thermalSimulationDynamics;
  timer.stop();
  timer.printSummary();
}
