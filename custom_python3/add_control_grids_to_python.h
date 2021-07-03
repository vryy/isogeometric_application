/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 11, 2017 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CONTROL_GRIDS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CONTROL_GRIDS_TO_PYTHON_H_INCLUDED

#include "includes/define_python.h"
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
void  IsogeometricApplication_AddControlGridsToPython(pybind11::module& m);
////////////////////////////////////////////////////////

}  // namespace Python.

} // Namespace Kratos

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CONTROL_GRIDS_TO_PYTHON_H_INCLUDED
