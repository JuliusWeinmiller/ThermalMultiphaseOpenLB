
/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

#include "descriptorFieldExtended.h"
#include "interactionPotentialExtended.h"
#include "modifiedShanChenForcedSingleComponentPostProcessor3D.h"
#include "modifiedShanChenForcedSingleComponentPostProcessor2D.h"
#include "pseudopotentialAdvectionDiffusionDynamics.h"
#include "modifiedGuoBGKDynamics.h"
#include "maxwellConstruction.h"
#include "fluidProperties.h"
#include "combinedShanChenThermalCouplingPostProcessor3D.h"
#include "combinedShanChenThermalCouplingPostProcessor2D.h"
#include "combinedShanChenThermalMRTCouplingPostProcessor3D.h"
#include "combinedShanChenThermalMRTCouplingPostProcessor2D.h"
#include "combinedShanChenSemiHybridThermalCouplingPostProcessor3D.h"
#include "combinedShanChenSemiHybridThermalCouplingPostProcessor2D.h"
#include "combinedShanChenSemiHybridThermalMRTCouplingPostProcessor3D.h"
#include "combinedShanChenSemiHybridThermalMRTCouplingPostProcessor2D.h"

#include "superAverage2D.h"
#include "blockAverage2D.h"
#include "blockIntegralF2D.h"

#include "miscFeatures.h"
