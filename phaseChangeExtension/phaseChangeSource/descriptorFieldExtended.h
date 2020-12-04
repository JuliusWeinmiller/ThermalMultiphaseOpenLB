/*  OpenLB Extension:
 *  Phase-change multiphase flow  
 *  Written by Julius Weinmiller for his thesis
 *  https://www.openlb.net/forum/users/julius-weinmiller/
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 */

#ifndef DESCRIPTOR_FIELD_EXTENDED_H
#define DESCRIPTOR_FIELD_EXTENDED_H


#include "dynamics/descriptorField.h"

namespace olb {
namespace descriptors {

struct PSI_PSEUDO_RHO      : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct PHI_THERMAL         : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct PREV_T_V            : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct PREV_PHI            : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };

} // end descriptors namespace
} // end olb namespace

#endif