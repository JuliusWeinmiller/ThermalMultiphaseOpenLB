/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

#ifndef MISC_FEATURES_H
#define MISC_FEATURES_H

namespace olb{

template <typename T, typename DESCRIPTOR, typename ThermalLattice>
void ThermalUnitConverter<T, DESCRIPTOR, ThermalLattice>::write(std::string const& title) const
{
  std::string dataFile = singleton::directories().getLogOutDir() + title + ".dat";

  if (singleton::mpi().isMainProcessor())
  {
    std::ofstream fout;
    fout.open(dataFile.c_str(), std::ios::trunc);

    fout << "----------------- UnitConverter information -----------------" << std::endl;
    fout << "-- Parameters:" << std::endl;
    fout << "Resolution:                                 N=                              " << this->getResolution() << std::endl;
    fout << "Lattice velocity:                           latticeU=                       " << this->getCharLatticeVelocity() << std::endl;
    fout << "Lattice relaxation frequency:               omega=                          " << this->getLatticeRelaxationFrequency() << std::endl;
    fout << "Lattice relaxation time:                    tau=                            " << this->getLatticeRelaxationTime() << std::endl;
    fout << "Thermal Lattice relaxation frequency:       omega_AD=                       " << this->getLatticeThermalRelaxationFrequency() << std::endl;
    fout << "Thermal Lattice relaxation time:            tau_AD=                         " << this->getLatticeThermalRelaxationTime() << std::endl;
    fout << "Characteristical length(m):                 charL=                          " << this->getCharPhysLength() << std::endl;
    fout << "Characteristical speed(m/s):                charU=                          " << this->getCharPhysVelocity() << std::endl;
    fout << "Phys. kinematic viscosity(m^2/s):           charNu=                         " << this->getPhysViscosity() << std::endl;
    fout << "Phys. density(kg/m^d):                      charRho=                        " << this->getPhysDensity() << std::endl;
    fout << "Characteristical pressure(N/m^2):           charPressure=                   " << this->getCharPhysPressure() << std::endl;
    fout << "Reynolds number:                            reynoldsNumber=                 " << this->getReynoldsNumber() << std::endl;

    fout << "-------------------------------------------------------------" << std::endl;

    fout << "----------------- ThermalUnitConverter information -----------------" << std::endl;
    fout << "-- Parameters:" << std::endl;
    fout << "Phys. Delta X(m):                           physDeltaX=                     " << this->getPhysDeltaX() << std::endl;
    fout << "Phys. Delta T(s):                           physDeltaT=                     " << this->getPhysDeltaT() << std::endl;
    fout << "Phys. Thermal Conductivity(W/m/K):          physThermalCondcticity=         " << getThermalConductivity() << std::endl;
    fout << "Phys. specific Heat Capacity(J/kg/K):       physSpecificHeatCapacity=       " << getPhysSpecificHeatCapacity() << std::endl;
    fout << "Phys. Low Temperature(K):                   physLowTemperature=             " << getCharPhysLowTemperature() << std::endl;
    fout << "Phys. High Temperature(K):                  physHighTemperature=            " << getCharPhysHighTemperature() << std::endl;
    fout << "Prandtl number:                             prandtlNumber=                  " << getPrandtlNumber() << std::endl;


    fout << "-------------------------------------------------------------" << std::endl;

    fout << "----------------- Conversion factors:-----------------" << std::endl;
    fout << "Voxel length(m):                            physDeltaX=                     " << this->getConversionFactorLength() << std::endl;
    fout << "Time step(s):                               physDeltaT=                     " << this->getConversionFactorTime() << std::endl;
    fout << "Velocity factor(m/s):                       physVelocity=                   " << this->getConversionFactorVelocity() << std::endl;
    fout << "Density factor(kg/m^3):                     physDensity=                    " << this->getConversionFactorDensity() <<  std::endl;
    fout << "Mass factor(kg):                            physMass=                       " << this->getConversionFactorMass() << std::endl;
    fout << "Viscosity factor(m^2/s):                    physViscosity=                  " << this->getConversionFactorViscosity() << std::endl;
    fout << "Force factor(N):                            physForce=                      " << this->getConversionFactorForce() << std::endl;
    fout << "Pressure factor(N/m^2):                     physPressure=                   " << this->getConversionFactorPressure() << std::endl;

    fout << "-------------------------------------------------------------" << std::endl;

    fout << "----------------- ThermalConversion factors:-----------------" << std::endl;
    fout << "Temperature(K):                             temperature=                    " << getConversionFactorTemperature() << std::endl;
    fout << "Thermal Diffusivity(m^2/s):                 physThermalDiffusivity=         " << getConversionFactorThermalDiffusivity() << std::endl;
    fout << "Specific Heat Capacity(J/kg):               physSpecificHeatCapacity=       " << getConversionFactorSpecificHeatCapacity() << std::endl;
    fout << "Thermal Conductivity(W/m/K):                physThermalCondcticity=         " << getConversionFactorThermalConductivity() <<  std::endl;
    fout << "Heat Flux(W):                               physHeatFlux=                   " << getConversionFactorHeatFlux() << std::endl;

    fout << "-------------------------------------------------------------" << std::endl;

    fout.close();
  }
}

/// Implementation of "bounce-back" dynamics
/** This is a very popular way to implement no-slip boundary conditions,
 * because the dynamics are independent of the orientation of the boundary.
 * 
 * The setDensity function is added such that the fixed density can be changed.
 * 
 * By default the varibale _rho of BB is private. This needs to be changed to protected.
 * 
 */
template<typename T, typename DESCRIPTOR>
class BounceBackVariableRho : public BounceBack<T,DESCRIPTOR> {
public:
  BounceBackVariableRho(T rho): BounceBack<T,DESCRIPTOR>(rho){};

  void setDensity( T newRho ){
    this->_rho = newRho;
  };

};

}// end namespace
#endif