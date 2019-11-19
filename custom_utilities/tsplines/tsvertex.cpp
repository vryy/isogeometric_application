//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Apr 2015 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_utilities/tsplines/tsvertex.h"

namespace Kratos
{

const int TsVertex::UNDEFINED_JOINT = -2;
const int TsVertex::BORDER_JOINT = -1;
const int TsVertex::NORMAL_JOINT = 0;
const int TsVertex::T_JOINT_LEFT = 1;
const int TsVertex::T_JOINT_RIGHT = 2;
const int TsVertex::T_JOINT_UP = 3;
const int TsVertex::T_JOINT_DOWN = 4;
const int TsVertex::L_JOINT_LEFT_DOWN = 5;
const int TsVertex::L_JOINT_LEFT_UP = 6;
const int TsVertex::L_JOINT_RIGHT_DOWN = 7;
const int TsVertex::L_JOINT_RIGHT_UP = 8;
const int TsVertex::I_JOINT_LEFT_RIGHT = 9;
const int TsVertex::I_JOINT_UP_DOWN = 10;

}// namespace Kratos.
