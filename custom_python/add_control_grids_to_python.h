/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        KratosIsogeometricApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 11, 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CONTROL_GRIDS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CONTROL_GRIDS_TO_PYTHON_H_INCLUDED

#include "custom_utilities/control_grid.h"

namespace Kratos
{

namespace Python
{

///////////////////////////////////////////////////////

template<typename TDataType>
TDataType ControlGrid_GetItem(ControlGrid<TDataType>& rDummy, int index)
{
    return rDummy.GetData(index);
}

template<typename TDataType>
void ControlGrid_SetItem(ControlGrid<TDataType>& rDummy, int index, const TDataType& value)
{
    rDummy.SetData(index, value);
}

////////////////////////////////////////////////////////
void  IsogeometricApplication_AddControlGridsToPython();
////////////////////////////////////////////////////////

}  // namespace Python.

} // Namespace Kratos

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CONTROL_GRIDS_TO_PYTHON_H_INCLUDED
