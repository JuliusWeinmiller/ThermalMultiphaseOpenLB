
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

#include "interactionPotentialExtended.hh"
#include "modifiedShanChenForcedSingleComponentPostProcessor3D.hh"
#include "modifiedShanChenForcedSingleComponentPostProcessor2D.hh"
#include "pseudopotentialAdvectionDiffusionDynamics.hh"
#include "modifiedGuoBGKDynamics.hh"
#include "maxwellConstruction.hh"
#include "fluidProperties.hh"
#include "combinedShanChenThermalCouplingPostProcessor3D.hh"
#include "combinedShanChenThermalCouplingPostProcessor2D.hh"
#include "combinedShanChenThermalMRTCouplingPostProcessor3D.hh"
#include "combinedShanChenThermalMRTCouplingPostProcessor2D.hh"
#include "combinedShanChenSemiHybridThermalCouplingPostProcessor3D.hh"
#include "combinedShanChenSemiHybridThermalCouplingPostProcessor2D.hh"
#include "combinedShanChenSemiHybridThermalMRTCouplingPostProcessor3D.hh"
#include "combinedShanChenSemiHybridThermalMRTCouplingPostProcessor2D.hh"

#include "superAverage2D.hh"
#include "blockAverage2D.hh"
#include "blockIntegralF2D.hh"
