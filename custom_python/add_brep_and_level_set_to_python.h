// see brep_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Feb 2021 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_BREP_AND_LEVEL_SET_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_BREP_AND_LEVEL_SET_TO_PYTHON_H_INCLUDED


// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  IsogeometricApplication_AddBRepAndLevelSetToPython();

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_BREP_AND_LEVEL_SET_TO_PYTHON_H_INCLUDED  defined
