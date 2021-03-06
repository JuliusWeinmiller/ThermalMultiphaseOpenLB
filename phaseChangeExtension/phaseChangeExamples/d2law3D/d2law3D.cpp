
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/* D2Law.cpp:
 * This simulation will combine both AD and NS lattice for a multiphase problem
 * A droplet is spended without body forces
 * periodic boundaries are applied in all three axis
 * The termpature is increased and evaporation rate recorded
 * The results of the simulations show that the simulations follow the 
 * well known D2 evaporative law. 
 */


#include "olb3D.h"
#include "olb3D.hh"   // include full template code
#include "phaseChangeExtension.h"
#include "phaseChangeExtension.hh"

using namespace olb;
// #define ThermalMRT
// WARNING
// For 3D simulations, thermalMRT does not work
// This is because of odd definiton of the D3Q7 velocity set
// The MRT dynamics would need to use the index 1,2,4 and not 1,2,3 as it is currently doing
// see: https://www.openlb.net/forum/topic/d3q7-mrt-velocity-directions/

typedef double SimType;
#define NSDESCRIPTOR olb::descriptors::D3Q19<olb::descriptors::VELOCITY,        \
                                             olb::descriptors::FORCE,           \
                                             olb::descriptors::EXTERNAL_FORCE,  \
                                             olb::descriptors::PSI_PSEUDO_RHO>          

#ifdef ThermalMRT
#define TDESCRIPTOR olb::descriptors::D3Q7<olb::descriptors::tag::MRT,      \
                                           olb::descriptors::VELOCITY,      \
                                           olb::descriptors::PHI_THERMAL,   \
                                           olb::descriptors::PREV_T_V,      \
                                           olb::descriptors::PREV_PHI>
#else
#define TDESCRIPTOR olb::descriptors::D3Q7<olb::descriptors::VELOCITY,      \
                                           olb::descriptors::PHI_THERMAL,   \
                                           olb::descriptors::PREV_T_V,      \
                                           olb::descriptors::PREV_PHI>
#endif


// Parameters for the simulation setup

// Resolutions
const int Nx = 75;       // spacial resolution of the model
const SimType Nt = 7;         // temporal resolution of the model 

const SimType physLength = 100e-6; // m
const SimType lx0   = 2*physLength;    // length of channel
const SimType ly0   = 2*physLength;     // height of channel
const SimType lz0   = 2*physLength;     // width of channel
const SimType conversionLength = physLength/Nx ;

const SimType sphereRadius = physLength/5.;
const SimType sphereInterfaceWidth = sphereRadius/4;
const SimType sphereOrigin[3] = {lx0/2., ly0/2., lz0/2};

const SimType rhoVapourMultiplier = 1.0;
const SimType rhoLiquidMultiplier = 0.99;

const SimType physViscosity = 2e-5; //  m2 / s

const SimType maxVelocity = 20;
const SimType charTime = physLength/maxVelocity;


const SimType endTimeInit = 0.1*charTime;        // Time to allow droplet to come to rest
const SimType endTimeEvap = 50.0*charTime;     // max. simulation time in s, SI unit
const SimType maxPhysT_ms = endTimeEvap*1000;  //

const SimType dx = physLength/ Nx;
const SimType dt = charTime/Nx/Nt;



////
// Thermal constants
const SimType TempCrit = 647; // [K] critical temperature of water in Kelvin

const SimType specificGasConstantVapour = 461.52;  
const SimType propertyTemperatureRatio = 0.6;
const SimType physWaterTemp = TempCrit * propertyTemperatureRatio - 273.15;
// http://www.thermopedia.com/content/1150/
// Table 8 for thermal conductivity
// Table 4 for specific heat at constant pressure
// Table 1 for density
const SimType physDensityLiquid = 947.13;
const SimType specificHeatCapacityVolumeLiquid = 3.8e3;    // J / kg K REAL
const SimType thermalConductivityLiquid = 682e-3; // W / m K REAL
const SimType physDensityVapour = 0.9643;
const SimType specificHeatCapacityVolumeVapour = 2.3e3 - specificGasConstantVapour; // Cv = Cp - Rsp REAL
const SimType thermalConductivityVapour =  31e-3;         // W / m K  REAL

const SimType liquidDensityRatio = physDensityLiquid/8.72518;


const SimType specificHeatCapacityVolume = 4e3;    // J / kg K
const SimType thermalConductivity = 4; // W / m K      


const int materialBase = 0;
const int materialFluid = 1;
const int materialBorder = 2;

void prepareGeometry( SuperGeometry3D<SimType>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( materialBase, materialBorder );

  Vector<SimType,3> extend( lx0, ly0, lz0 );
  Vector<SimType,3> origin;
  
  //clout << "Prepare Internal Box ..." << std::endl;
  origin = {conversionLength, conversionLength, conversionLength};
  extend = {lx0-2*conversionLength, ly0-2*conversionLength, lz0-2*conversionLength}; 
  olb::IndicatorCuboid3D<SimType> internalBox( extend, origin );
  superGeometry.rename( materialBorder, materialFluid, internalBox );

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
                     olb::SuperLattice3D<SimType, NSDESCRIPTOR>& NSlattice,
                     olb::SuperLattice3D<SimType, TDESCRIPTOR>& ADlattice,
                     olb::Dynamics<SimType, NSDESCRIPTOR>& bulkDynamics,
                     olb::Dynamics<SimType, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     olb::sOnLatticeBoundaryCondition3D<SimType,TDESCRIPTOR>& TboundaryCondition,
                     SuperGeometry3D<SimType>& superGeometry,
                     SimType TempInitial,
                     SimType rhoVapor, SimType rhoLiquid )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << endl;

  const SimType omega = converter.getLatticeRelaxationFrequency();
  const SimType Tomega  = converter.getLatticeThermalRelaxationFrequency();

  clout << "NS omega = " << omega << endl;
  clout << "T  omega = " << Tomega << endl;

  auto bulkIndicator = superGeometry.getMaterialIndicator({materialFluid, materialBorder});

  clout << "Define NS Dynamics ..." << endl;
  NSlattice.defineDynamics( superGeometry, materialBase, &olb::instances::getNoDynamics<SimType, NSDESCRIPTOR>() );
  NSlattice.defineDynamics( bulkIndicator, &bulkDynamics );   

  clout << "Define AD Dynamics ..." << endl;
  ADlattice.defineDynamics( superGeometry, materialBase, &olb::instances::getNoDynamics<SimType, TDESCRIPTOR>());
  ADlattice.defineDynamics( bulkIndicator, &advectionDiffusionBulkDynamics );

  TboundaryCondition.addTemperatureBoundary(superGeometry, materialBorder, Tomega);


  clout << "Define Initial Conditions ..." << endl;
  // Fluid Initial conditions
  olb::AnalyticalConst3D<SimType,SimType> zeroVelocity( 0,0,0 );

  AnalyticalConst3D<SimType,SimType> noiseVap( rhoVapor*0.1 );
  AnalyticalConst3D<SimType,SimType> noiseLiq( rhoLiquid*0.1 );
  AnalyticalRandom3D<SimType,SimType> random;

  olb::AnalyticalConst3D<SimType,SimType> constRhoVapour( rhoVapor  * rhoVapourMultiplier );
  olb::AnalyticalConst3D<SimType,SimType> constRhoLiquid( rhoLiquid * rhoLiquidMultiplier );


  // Temperature Inital Conditions
  olb::AnalyticalConst3D<SimType,SimType> T_init(   converter.getLatticeTemperature(TempInitial));//

  // Random density distribution across inlet
  olb::AnalyticalIdentity3D<SimType,SimType> noiseIndicatorVap( random*noiseVap  );
  olb::AnalyticalIdentity3D<SimType,SimType> noiseIndicatorLiq( random*noiseLiq  );
  olb::AnalyticalIdentity3D<SimType,SimType> rhoBase( constRhoVapour );
  olb::SmoothIndicatorSphere3D<SimType,SimType> sphere( sphereOrigin, sphereRadius, sphereInterfaceWidth);

  olb::AnalyticalIdentity3D<SimType,SimType> sphereIndicator( sphere  );

  olb::AnalyticalIdentity3D<SimType,SimType> rho3( rhoBase + noiseIndicatorVap + sphereIndicator*( constRhoLiquid - constRhoVapour + noiseIndicatorLiq) );

  //Initialize all values of distribution functions to their local equilibrium

  // Setting initial psi values
  olb::AnalyticalConst3D<SimType,SimType> initialPsi( 1.0 );
  NSlattice.defineField<olb::descriptors::PSI_PSEUDO_RHO>( bulkIndicator, initialPsi ) ; 


  // Fluid init
  NSlattice.defineRhoU( bulkIndicator, rho3, zeroVelocity );
  NSlattice.iniEquilibrium( bulkIndicator, rho3, zeroVelocity );

  // Thermal init
  ADlattice.defineRho(bulkIndicator, T_init);
  ADlattice.iniEquilibrium(bulkIndicator, T_init, zeroVelocity);


  // Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( olb::ThermalUnitConverter<SimType, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                        olb::SuperLattice3D<SimType, NSDESCRIPTOR>& NSlattice,
                        olb::SuperLattice3D<SimType, TDESCRIPTOR>& ADlattice,
                        int iT,
                        SuperGeometry3D<SimType>& superGeometry,
                        SimType newTempRatio )
{
  olb::OstreamManager clout( std::cout,"setBoundaryValues" );
}

// Output to console and files
void getResults( olb::SuperLattice3D<SimType, NSDESCRIPTOR>& NSlattice,
                 olb::SuperLattice3D<SimType, TDESCRIPTOR>& ADlattice,
                 olb::SuperGeometry3D<SimType>& superGeometry,
                 olb::ThermalUnitConverter<SimType, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                 int iT,
                 Timer<SimType>& timer,
                 bool forcedSave=false)
{
  olb::OstreamManager clout( std::cout,"getResults" );

  olb::SuperVTMwriter3D<SimType> vtmWriterNS( "fluidLatticeResults" );
  olb::SuperLatticePhysVelocity3D<SimType, NSDESCRIPTOR> velocity( NSlattice, converter );
  olb::SuperLatticeDensity3D<SimType, NSDESCRIPTOR> density( NSlattice );
  olb::SuperLatticeVelocity3D<SimType, NSDESCRIPTOR> latticeVelocity( NSlattice );
  olb::SuperLatticePhysPressure3D<SimType, NSDESCRIPTOR> pressure( NSlattice, converter );

  olb::SuperVTMwriter3D<SimType> vtmWriterThermal( "thermalLatticeResults" );
  olb::SuperLatticeDensity3D<SimType, TDESCRIPTOR> latticeTemperature( ADlattice );
  olb::SuperLatticePhysTemperature3D<SimType, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
  
  vtmWriterNS.addFunctor( velocity );
  vtmWriterNS.addFunctor( latticeVelocity );
  vtmWriterNS.addFunctor( pressure );
  vtmWriterNS.addFunctor( density );
  vtmWriterThermal.addFunctor( latticeTemperature );
  vtmWriterThermal.addFunctor( temperature );

  const int vtkIter  = converter.getLatticeTime( 250*charTime/Nx/Nt );
  const int statIter = converter.getLatticeTime( 250*charTime/Nx/Nt );
  const int ppmIter = converter.getLatticeTime( 250*charTime/Nx/Nt ); 

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    olb::SuperLatticeGeometry3D<SimType, NSDESCRIPTOR> geometry( NSlattice, superGeometry );
    olb::SuperLatticeCuboid3D<SimType, NSDESCRIPTOR> cuboid( NSlattice );
    olb::SuperLatticeRank3D<SimType, NSDESCRIPTOR> rank( NSlattice );
    vtmWriterNS.write( geometry );
    vtmWriterNS.write( cuboid );
    vtmWriterNS.write( rank );
    vtmWriterNS.createMasterFile();
    vtmWriterThermal.createMasterFile();
  }

  // Writes the vtk files
  if ( iT%vtkIter==0 || forcedSave) {
    vtmWriterNS.write( iT );
    vtmWriterThermal.write( iT );
  }

  // Writes the ppm files
  if ( iT%ppmIter==0 || (iT<=100 && (iT%10==0) ) || forcedSave) {
    olb::BlockReduction3D2D<SimType> planeReductionRho( density, {0, 1, 0} );
    olb::BlockReduction3D2D<SimType> planeReductionT( temperature, {0, 1, 0} );
    olb::heatmap::plotParam<SimType> thermalPlotParam;
    thermalPlotParam.name = "physTemperature";
    thermalPlotParam.minValue = converter.getCharPhysLowTemperature()-10;
    thermalPlotParam.maxValue = converter.getCharPhysHighTemperature()+10;
    // write output as JPEG
    olb::heatmap::write(planeReductionRho, iT);
    olb::heatmap::write(planeReductionT, iT, thermalPlotParam);
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

void writeGeometry( olb::SuperLattice3D<SimType, NSDESCRIPTOR>& NSlattice,
                    olb::SuperGeometry3D<SimType>& superGeometry)
{
  olb::SuperVTMwriter3D<SimType> vtmWriter( "D2Results" );
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
  olb::OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);


  // == Define simulations constants
  const SimType sigma = 0.31;
  const SimType eosA  = 0.5/49.0;
  std::vector<SimType> rho0 ={1.0};
  const SimType G     = -1;

  SimType molarMass = 18.0152;


  SimType temperatureRatio = 0.7;
  SimType newTempRatio =     0.8;

  SimType physThermalConductivity;
  for( physThermalConductivity = 100.0e-3; physThermalConductivity<=500.01e-3; physThermalConductivity+=200.0e-3 ){

    std::ostringstream outputDirStream;
    outputDirStream << std::setprecision(3);
    outputDirStream << "./results3D/";
    mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    outputDirStream << "lambda_" << physThermalConductivity << "/";
    olb::singleton::directories().setOutputDir( outputDirStream.str() );

    

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
    olb::CuboidGeometry3D<SimType> cuboidGeometry( cuboid, conversionLength, noOfCuboids );
    cuboidGeometry.setPeriodicity( true, true, true );

    // Instantiation of a loadBalancer
    olb::HeuristicLoadBalancer<SimType> loadBalancer( cuboidGeometry );

    // Instantiation of a superGeometry
    olb::SuperGeometry3D<SimType> superGeometry( cuboidGeometry, loadBalancer, 2 );

    prepareGeometry( superGeometry );

  // === 3rd Step: Prepare Lattice ===
    olb::SuperLattice3D<SimType, TDESCRIPTOR> ADlattice(superGeometry);
    olb::SuperLattice3D<SimType, NSDESCRIPTOR> NSlattice(superGeometry);
    writeGeometry(NSlattice, superGeometry);
    
    // Pseudopotential multiphase single component dynamics
    // A dimensionless R is needed to recover correct latent heat
    SimType Rdimless = 8314.0/molarMass /(dx*dx/dt/dt) * TempCrit;

    // Pseudopotential multiphase single component dynamics
    olb::PengRobinsonExt<SimType,SimType> interactionPotential( G, 0.3443, eosA, 2./21., temperatureRatio, Rdimless );
    interactionPotential.write();

    SimType rho[2];
    olb::MaxwellConstruction<SimType, decltype(interactionPotential)> maxwellConstruction( interactionPotential );
    maxwellConstruction(rho, temperatureRatio); //temperatureRatio
    SimType densityLatticeVapour = rho[1]; 
    SimType densityLatticeLiquid = rho[0];
    clout << "Density liquid: " <<  densityLatticeLiquid  << std::endl;
    clout << "Density vapour: " <<  densityLatticeVapour  << std::endl;

    // == Define converter
    olb::ThermalUnitConverter<SimType, NSDESCRIPTOR, TDESCRIPTOR> converter(
      (SimType)   dx,                // physDeltaX: spacing between two lattice cells in __m__
      (SimType)   dt,                // physDeltaT = charLatticeVelocity / charPhysVelocity * dx
      (SimType)   physLength,        // charPhysLength: reference length of simulation geometry
      (SimType)   maxVelocity,       // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (SimType)   physViscosity,     // physViscosity: physical kinematic viscosity in __m^2 / s__
      (SimType)   liquidDensityRatio,          // physDensity: physical density in __kg / m^3__
      (SimType)   thermalConductivity,           // physThermalConductivity
      (SimType)   specificHeatCapacityVolume,    // specificHeatCapacity
      (SimType)   0.0,                           // physThermalExpansionCoefficient
      (SimType)   temperatureRatio*TempCrit,                 // charPhysLowTemperature
      (SimType)   newTempRatio*TempCrit           // charPhysHighTemperature
    );
      // Prints the converter log as console output
    converter.print();
    converter.write();

    olb::FluidProperties<SimType> fluidProperties( 
                      densityLatticeVapour, densityLatticeLiquid,
                      converter.getLatticeSpecificHeatCapacity( specificHeatCapacityVolume ),
                      converter.getLatticeSpecificHeatCapacity( specificHeatCapacityVolume ),
                      converter.getLatticeThermalConductivity( physThermalConductivity ),
                      converter.getLatticeThermalConductivity( physThermalConductivity ),
                      converter.getCharPhysLowTemperature(), converter.getCharPhysHighTemperature(), TempCrit );

    fluidProperties.write();

    // Modified Guo forcing dynamics
    olb::ForcedModifiedGuoBGKdynamics<SimType, NSDESCRIPTOR> fluidSimulationDynamics (
    converter.getLatticeRelaxationFrequency(),
    olb::instances::getBulkMomenta<SimType,NSDESCRIPTOR>(), sigma 
    );

    // Modified Advection Diffusion dynamics
    olb::Dynamics<SimType, TDESCRIPTOR>* thermalSimulationDynamics;
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
    olb::sOnLatticeBoundaryCondition3D<SimType,TDESCRIPTOR> TboundaryCondition(ADlattice);
    olb::createAdvectionDiffusionBoundaryCondition3D<SimType,TDESCRIPTOR>(TboundaryCondition);  // OnLatticeAdvectionDiffusionBoundaryCondition3D


    prepareLattice( converter,
                    NSlattice, ADlattice,
                    fluidSimulationDynamics, *thermalSimulationDynamics,
                    TboundaryCondition,
                    superGeometry,
                    TempCrit * temperatureRatio,
                    densityLatticeVapour, densityLatticeLiquid );

  // === 4th Step: Prepare Lattice Coupling
    olb::LatticeCouplingGenerator3D<SimType, NSDESCRIPTOR> *thermalMultiphaseCoupling;
    #ifdef ThermalMRT
      thermalMultiphaseCoupling = new olb::CombinedShanChenThermalMRTCouplingGenerator3D<SimType, NSDESCRIPTOR>(
          0, converter.getLatticeLength(lx0), 
          0, converter.getLatticeLength(ly0),
          0, converter.getLatticeLength(lz0), 
          G, converter.getLatticeThermalRelaxationFrequency(),
          0, 0, {0,0,1},
          fluidProperties, interactionPotential);
    #else
      thermalMultiphaseCoupling = new olb::CombinedShanChenThermalCouplingGenerator3D<SimType, NSDESCRIPTOR>(
          0, converter.getLatticeLength(lx0), 
          0, converter.getLatticeLength(ly0),
          0, converter.getLatticeLength(lz0), 
          G, converter.getLatticeThermalRelaxationFrequency(),
          0, 0, {0,0,1},
          fluidProperties, interactionPotential);
    #endif

    NSlattice.addLatticeCoupling(*thermalMultiphaseCoupling, ADlattice);

  // === 5th Step: Main Loop with Timer ===
    clout << "starting simulation..." << endl;
    olb::util::Timer<SimType> timer( converter.getLatticeTime( endTimeEvap ), superGeometry.getStatistics().getNvoxel() );
    timer.start();
    int latticeInitTime = converter.getLatticeTime( endTimeInit );
    clout << "Initial timestep end: " << latticeInitTime << std::endl;

    //setBoundaryValues( converter, NSlattice, ADlattice, 0, superGeometry, temperatureRatio);
    getResults(NSlattice, ADlattice, superGeometry, converter, 0, timer );

    olb::AnalyticalConst3D<SimType,SimType> zeroVelocity( 0,0,0 );
    olb::AnalyticalConst3D<SimType,SimType> T_const( converter.getLatticeTemperature(TempCrit * temperatureRatio));
    for ( int iT = 0; iT < latticeInitTime; ++iT ){

      // === 6th Step: Definition of Initial and Boundary Conditions ===

      // === 7th Step: Collide and Stream Execution ===
    
      NSlattice.collideAndStream();
      ADlattice.collideAndStream();
      NSlattice.communicate();
      ADlattice.communicate();
      NSlattice.executeCoupling();

      // === Reset temperature field
      ADlattice.defineRho(superGeometry, materialFluid, T_const);
      ADlattice.iniEquilibrium(superGeometry, materialFluid, T_const, zeroVelocity);

      // === 8th Step: Computation and Output of the Results ===
      getResults(NSlattice, ADlattice, superGeometry, converter, iT+1, timer );
      if ( std::isnan( NSlattice.getStatistics().getAverageRho() ) ) { 
                clout << "Terminating this run at step " << iT+1 << " => density is NAN" << std::endl;
                getResults(NSlattice, ADlattice, superGeometry, converter, iT+1, timer, true );
                break;
      }
    }

    olb::AnalyticalConst3D<SimType,SimType> T_new( converter.getLatticeTemperature(TempCrit * newTempRatio));//
    ADlattice.defineRho(superGeometry, materialBorder, T_new);

    for ( int iT = latticeInitTime; iT < converter.getLatticeTime( endTimeEvap ); ++iT ) {

      // === 6th Step: Definition of Initial and Boundary Conditions ===

      // === 7th Step: Collide and Stream Execution ===
    
      NSlattice.collideAndStream();
      ADlattice.collideAndStream();
      NSlattice.communicate();
      ADlattice.communicate();
      NSlattice.executeCoupling();

      // === 8th Step: Computation and Output of the Results ===
      getResults(NSlattice, ADlattice, superGeometry, converter, iT+1, timer );
      if ( std::isnan( NSlattice.getStatistics().getAverageRho() ) ) { 
                clout << "Terminating this run at step " << iT+1 << " => density is NAN" << std::endl;
                getResults(NSlattice, ADlattice, superGeometry, converter, iT+1, timer, true );
                break;
      }
    }
    delete thermalSimulationDynamics; // free previous alloc to avoid memory leak
    delete thermalMultiphaseCoupling;
    timer.stop();
    timer.printSummary();
  } // end for loop
}
