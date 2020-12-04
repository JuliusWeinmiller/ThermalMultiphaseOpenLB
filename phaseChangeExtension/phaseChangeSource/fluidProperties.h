
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
 * -- header file.
 */
#ifndef FLUID_PROPERTIES_H
#define FLUID_PROPERTIES_H

namespace olb {

template <typename T>
class FluidProperties{
    public:
        FluidProperties(T densityVapour, T densityLiquid,
                          T specificHeatVapour, T specificHeatLiquid,
                          T conductivityVapour, T conductivityLiquid,
                          T tempLower, T tempHigher, T tempCritical);


        // SpecificHeat( T ) -> T
        // Returns the specific heat depending on the density
        T specificHeat(const T latticeDensity);

        // conductivity( T ) -> T
        // Returns the conductivity depending on the density
        T conductivity(const T latticeDensity);

        // density( T ) -> T
        // Returns the conductivity depending on the density
        T density(const T latticeDensity);

        // temperatureRatio( T ) -> T
        // Returns the temperature ratio depending on the lattice temperature
        // Should be fewer computations than using converter to convert from lattice to physical to ratio 
        T temperatureRatio(const T latticeTemperature);


        T getDensityVapour();
        T getDensityLiquid();


        T getTempConvertionRatio();

        void verify();
        void print();
        void write(std::string const& title = "FluidProperties");

    protected:
        static T scaleProperty(const T localScale, const T propertyVapour, const T propertyLiquid);        
        
        // ScaleDensity( T ) -> T
        // Returns 1 for liquid and 0 for vapour 
        T scaleDensity(const T localDensity);

    private:
        const T _densityVapour;
        const T _densityLiquid;
        const T _specificHeatVapour;
        const T _specificHeatLiquid;
        const T _conductivityVapour;
        const T _conductivityLiquid;
        const T _tempLowerRatio;
        const T _tempRatioDifference;
        const T _tempCrit;
};

} // end namespace

#endif
