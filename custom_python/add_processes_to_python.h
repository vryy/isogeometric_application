/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 Jan 2015 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_PROCESSES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_PROCESSES_TO_PYTHON_H_INCLUDED



// System includes


// External includes
#include <pybind11/pybind11.h>

// Project includes


namespace Kratos
{

namespace Python
{

void  IsogeometricApplication_AddProcessesToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_PROCESSES_TO_PYTHON_H_INCLUDED  defined
