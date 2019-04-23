/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Feb 18, 2019 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_POINT_BASED_BSPLINES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_POINT_BASED_BSPLINES_TO_PYTHON_H_INCLUDED

#include <pybind11/pybind11.h>

namespace Kratos
{

namespace Python
{

void  IsogeometricApplication_AddPBBSplinesToPython(pybind11::module& m);

}  // namespace Python.

} // Namespace Kratos

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_POINT_BASED_BSPLINES_TO_PYTHON_H_INCLUDED
