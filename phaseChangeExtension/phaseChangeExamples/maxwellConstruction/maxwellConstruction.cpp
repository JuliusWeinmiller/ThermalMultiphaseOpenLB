
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/* maxwellConstruction.cpp:
 * This file will give the predicted densities and print them to file
 */

// Either the olb2D or olb3D is needed for setup functionalities
#include "olb3D.h"
#include "olb3D.hh" 

#include "phaseChangeExtension.h"
#include "phaseChangeExtension.hh"


using namespace olb;


typedef double SimType;


int main( int argc, char* argv[] )
{

    // === 1st Step: Initialization ===
    olb::OstreamManager clout( std::cout,"main" );

    // Loop through variables
    SimType temperatureRatio;      

    const SimType constDerivRho[1] = {0.1};

    std::vector<SimType> temperatureRatios;   
    std::vector<SimType> liquidDensities; 
    std::vector<SimType> vapourDensities;
    std::vector<SimType> pressureDerivative;

    std::ostringstream outputDirStream;
    outputDirStream << std::setprecision(2);
    outputDirStream << "./tmp/";
    mkdir( outputDirStream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    // olb::singleton::directories().setOutputDir( outputDirStream.str() );

    for(temperatureRatio = 0.49; temperatureRatio<=0.851; temperatureRatio+=0.01 ){
    
        // Pseudopotential multiphase single component dynamics
        std::vector<SimType> rho0 ={1.0};
        const SimType G     = -1;

        // olb::PengRobinsonExt<SimType,SimType> interactionPotential( G, 0.3443, 0.5/49., 2./21., temperatureRatio );
        // olb::MaxwellConstruction<SimType, decltype(interactionPotential)> maxwellConstruction( interactionPotential);                      // For PR use default interval

        olb::CarnahanStarlingExt<SimType,SimType> interactionPotential( G, 0.5, 4, temperatureRatio );
        olb::MaxwellConstruction<SimType, decltype(interactionPotential)> maxwellConstruction( interactionPotential, 1.01, 1e22 );         // For CS use interval 1.01, 1e22


        
        maxwellConstruction.setAccuracy(1e-5);
        SimType rho[2];
        clout << "Performing maxwell construction with Tr = "<< temperatureRatio << std::endl;
        maxwellConstruction(rho, temperatureRatio);
        SimType densityLatticeVapour = rho[1]; 
        SimType densityLatticeLiquid = rho[0];
        clout << "Density liquid: " <<  densityLatticeLiquid  << std::endl;
        clout << "Density vapour: " <<  densityLatticeVapour  << std::endl;

        SimType calcPderiv[1];
        SimType passTr[1] = {temperatureRatio};
        interactionPotential.computeDerivativeConstRho(calcPderiv, constDerivRho, passTr );

        temperatureRatios.push_back(temperatureRatio);
        liquidDensities.push_back(densityLatticeLiquid); 
        vapourDensities.push_back(densityLatticeVapour);
        pressureDerivative.push_back(calcPderiv[0]);

    } // end for temperatureRatio

    std::string dataFile = singleton::directories().getLogOutDir() + "predictedDensities.dat";
    if (singleton::mpi().isMainProcessor())
    {
        std::ofstream fout;
        fout.open(dataFile.c_str(), std::ios::trunc);

        for(uint iter=0; iter<temperatureRatios.size(); iter++){
            fout << temperatureRatios[iter] << ", " << liquidDensities[iter] << ", " << vapourDensities[iter] << ", " << pressureDerivative[iter] << std::endl;
        }

        fout.close();
    }

} // end main
