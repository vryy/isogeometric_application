/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 18, 2013 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED

#include "includes/define_python.h"

namespace Kratos
{

namespace Python
{

void  IsogeometricApplication_AddBackendUtilitiesToPython(pybind11::module& m);
void  IsogeometricApplication_AddFrontendUtilitiesToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined
