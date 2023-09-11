// ****************************************************************************
//
//              SiMuScale - Multi-scale simulation framework
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package
// E-mail: simuscale-contact@lists.gforge.inria.fr
// Original Authors : Samuel Bernard, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

#ifndef SIMUSCALE_PARAMS_H__
#define SIMUSCALE_PARAMS_H__

// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>

#define NBR_SIGNALS 87//60 //33 //56
#define INTERCELLSIGNALS(...) enum class InterCellSignal { __VA_ARGS__ }; \
        inline const char *  InterCellSignalName() { return #__VA_ARGS__; }

#define NBR_CELLTYPES 7
#define CELLTYPES(...) typedef enum { __VA_ARGS__ } CellType; \
        inline const char *  CellTypeName() { return #__VA_ARGS__; }

/** Intercellular Signals  
 *  Add signals as needed, using format
 *  NEW_SIGNAL_NAME, \
 *  respecting the whitespaces after the comma 
 *  Increment NBR_SIGNALS
 **/
INTERCELLSIGNALS(\
CYCLE, \
CYCLE_X, \
CYCLE_Z, \
CLOCK, \
DEATH, \
NICHE, \
SPIKY_1, \
LYMPHOCYTE_TYPE, \
LYMPHOCYTE_P1, \
LYMPHOCYTE_P2, \
LYMPHOCYTE_P3, \
LYMPHOCYTE_P4, \
LYMPHOCYTE_P5, \
LYMPHOCYTE_P6, \
LYMPHOCYTE_P7, \
LYMPHOCYTE_P8, \
LYMPHOCYTE_P9, \
LYMPHOCYTE_P1_1, \
LYMPHOCYTE_P2_1, \
LYMPHOCYTE_P3_1, \
LYMPHOCYTE_P4_1, \
LYMPHOCYTE_P5_1, \
LYMPHOCYTE_P6_1, \
LYMPHOCYTE_P7_1, \
LYMPHOCYTE_P8_1, \
LYMPHOCYTE_P9_1, \
LYMPHOCYTE_P1_2, \
LYMPHOCYTE_P2_2, \
LYMPHOCYTE_P3_2, \
LYMPHOCYTE_P4_2, \
LYMPHOCYTE_P5_2, \
LYMPHOCYTE_P6_2, \
LYMPHOCYTE_P7_2, \
LYMPHOCYTE_P8_2, \
LYMPHOCYTE_P9_2, \
LYMPHOCYTE_P1_3, \
LYMPHOCYTE_P2_3, \
LYMPHOCYTE_P3_3, \
LYMPHOCYTE_P4_3, \
LYMPHOCYTE_P5_3, \
LYMPHOCYTE_P6_3, \
LYMPHOCYTE_P7_3, \
LYMPHOCYTE_P8_3, \
LYMPHOCYTE_P9_3, \
REMEMBER_DIVISION, \
LYMPHOCYTE_mRNA1, \
LYMPHOCYTE_mRNA2, \
LYMPHOCYTE_mRNA3, \
LYMPHOCYTE_mRNA4, \
LYMPHOCYTE_mRNA5, \
LYMPHOCYTE_mRNA6, \
LYMPHOCYTE_mRNA7, \
LYMPHOCYTE_mRNA8, \
LYMPHOCYTE_mRNA9, \
LYMPHOCYTE_mRNA1_1, \
LYMPHOCYTE_mRNA2_1, \
LYMPHOCYTE_mRNA3_1, \
LYMPHOCYTE_mRNA4_1, \
LYMPHOCYTE_mRNA5_1, \
LYMPHOCYTE_mRNA6_1, \
LYMPHOCYTE_mRNA7_1, \
LYMPHOCYTE_mRNA8_1, \
LYMPHOCYTE_mRNA9_1, \
LYMPHOCYTE_mRNA1_2, \
LYMPHOCYTE_mRNA2_2, \
LYMPHOCYTE_mRNA3_2, \
LYMPHOCYTE_mRNA4_2, \
LYMPHOCYTE_mRNA5_2, \
LYMPHOCYTE_mRNA6_2, \
LYMPHOCYTE_mRNA7_2, \
LYMPHOCYTE_mRNA8_2, \
LYMPHOCYTE_mRNA9_2, \
LYMPHOCYTE_mRNA1_3, \
LYMPHOCYTE_mRNA2_3, \
LYMPHOCYTE_mRNA3_3, \
LYMPHOCYTE_mRNA4_3, \
LYMPHOCYTE_mRNA5_3, \
LYMPHOCYTE_mRNA6_3, \
LYMPHOCYTE_mRNA7_3, \
LYMPHOCYTE_mRNA8_3, \
LYMPHOCYTE_mRNA9_3, \
APC_CONTACT, \
APC_DURATION, \
T_DURATION,	    \
LYMPHOCYTE_CONTACT, \
VOLUME, \
DIFFERENTIATION_STATE \
);

using CellFormalism = uint32_t;

/** Cell Types  
 *  Add cell types as needed, using format
 *  NEW_CELL_TYPE, \
 *  respecting the whitespaces after the comma 
 *  Increment NBR_CELLTYPES
 **/
CELLTYPES(\
STEM, \
CANCER, \
DIFFA, \
DIFFB, \
LYMPHOCYTE, \
APC, \
NICHE \
)

#endif // SIMUSCALE_PARAMS_H__
