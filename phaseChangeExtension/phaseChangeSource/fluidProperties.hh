
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

/** \file
 * A compact method to pass vapour & liquid properties to other functions
 * -- implementation file.
 */
#ifndef FLUID_PROPERTIES_HH
#define FLUID_PROPERTIES_HH

#include "fluidProperties.h"

namespace olb{

namespace fpf {
template <typename T>
T clamp(T x, T a, T b) {
  if (x < a) {
    return a;
  } else if (x > b) {
    return b;
  } else {
    return x;
  }
}
}

template <typename T>
FluidProperties<T>::FluidProperties(
    T densityVapour, T densityLiquid,
    T specificHeatVapour, T specificHeatLiquid,
    T conductivityVapour, T conductivityLiquid,
    T tempLower, T tempHigher, T tempCrit ):
    _densityVapour(densityVapour), _densityLiquid(densityLiquid), 
    _specificHeatVapour(specificHeatVapour), _specificHeatLiquid(specificHeatLiquid),
    _conductivityVapour(conductivityVapour), _conductivityLiquid(conductivityLiquid),
    _tempLowerRatio(tempLower/tempCrit), _tempRatioDifference((tempHigher-tempLower)/tempCrit),
    _tempCrit(tempCrit)
    { 
        assert( specificHeatLiquid >= specificHeatVapour);
        assert( tempHigher >= tempLower);
        assert( densityLiquid >= densityVapour);

        verify();
    }

template <typename T>
T FluidProperties<T>::scaleDensity(const T localDensity){
    return (fpf::clamp(localDensity, _densityVapour, _densityLiquid) - _densityVapour)/(_densityLiquid - _densityVapour);
}

template <typename T>
T FluidProperties<T>::density(const T localDensity){
    return fpf::clamp(localDensity, _densityVapour, _densityLiquid);
}

template <typename T>
T FluidProperties<T>::specificHeat(const T localDensity){
    return scaleProperty(scaleDensity(localDensity), _specificHeatVapour, _specificHeatLiquid );
}

template <typename T>
T FluidProperties<T>::conductivity(const T localDensity){
    return scaleProperty(scaleDensity(localDensity), _conductivityVapour, _conductivityLiquid );
}

template <typename T>
T FluidProperties<T>::scaleProperty(const T localScale, const T propertyVapour, const T propertyLiquid){
    return propertyVapour*(1-localScale) + propertyLiquid*localScale;
}

template <typename T>
T FluidProperties<T>::temperatureRatio(const T latticeTemp){
    return _tempLowerRatio + (latticeTemp - 0.5)*(_tempRatioDifference);
}

template <typename T>
T FluidProperties<T>::getDensityVapour(){return _densityVapour;}

template <typename T>
T FluidProperties<T>::getDensityLiquid(){return _densityLiquid;}

template <typename T>
T FluidProperties<T>::getTempConvertionRatio(){return _tempRatioDifference;}


template <typename T>
void FluidProperties<T>::verify(){

    T rhoV = getDensityVapour();
    T rhoL = getDensityLiquid();
    T cvV = specificHeat(getDensityVapour());
    T cvL = specificHeat(getDensityLiquid());
    T kV = conductivity(getDensityVapour());
    T kL = conductivity(getDensityLiquid());

    assert( olb::util::nearZero(_densityVapour      - rhoV  ) );
    assert( olb::util::nearZero(_densityLiquid      - rhoL  ) );
    assert( olb::util::nearZero(_specificHeatVapour - cvV  ) );
    assert( olb::util::nearZero(_specificHeatLiquid - cvL  ) );
    assert( olb::util::nearZero(_conductivityVapour - kV  ) );
    assert( olb::util::nearZero(_conductivityLiquid - kL  ) );

};

template <typename T>
void FluidProperties<T>::print(){
    OstreamManager clout( std::cout,"ScalingProperties" );

    clout<< "Simulation fluid properties: " <<std::endl;
    clout<< "Vapour Density:                              " << _densityVapour      << std::endl;
    clout<< "Vapour Specific heat at constant volume:     " << _specificHeatVapour << std::endl;
    clout<< "Vapour Thermal conductivity:                 " << _conductivityVapour << std::endl;
    clout<< "Vapour Thermal diffusivity (calculated):     " << _conductivityVapour/(_specificHeatVapour*_densityVapour) << std::endl;
    clout<< "Liquid Density:                              " << _densityLiquid      << std::endl;
    clout<< "Liquid Specific heat at constant volume:     " << _specificHeatLiquid << std::endl;
    clout<< "Liquid Thermal conductivity:                 " << _conductivityLiquid << std::endl;
    clout<< "Liquid Thermal diffusivity (calculated):     " << _conductivityLiquid/(_specificHeatLiquid*_densityLiquid) << std::endl;

};

template <typename T>
void FluidProperties<T>::write(std::string const& title){
    
    std::string dataFile = singleton::directories().getLogOutDir() + title + ".dat";

    if (singleton::mpi().isMainProcessor())
    {
        std::ofstream fout;
        fout.open(dataFile.c_str(), std::ios::trunc);

        fout<< "Simulation fluid properties: " <<std::endl;
        fout<< "Vapour Density:                              " << _densityVapour      << std::endl;
        fout<< "Vapour Specific heat at constant volume:     " << _specificHeatVapour << std::endl;
        fout<< "Vapour Thermal conductivity:                 " << _conductivityVapour << std::endl;
        fout<< "Vapour Thermal diffusivity (calculated):     " << _conductivityVapour/(_specificHeatVapour*_densityVapour) << std::endl;
        fout<< "Liquid Density:                              " << _densityLiquid      << std::endl;
        fout<< "Liquid Specific heat at constant volume:     " << _specificHeatLiquid << std::endl;
        fout<< "Liquid Thermal conductivity:                 " << _conductivityLiquid << std::endl;
        fout<< "Liquid Thermal diffusivity (calculated):     " << _conductivityLiquid/(_specificHeatLiquid*_densityLiquid) << std::endl;
        fout.close();
    }

};

}// end namespace
#endif



  


