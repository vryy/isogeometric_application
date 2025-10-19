/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        KratosIsogeometricApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 26 Jan 2021 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

#include "custom_python/add_brep_and_level_set_to_python.h"
#include "isogeometric_brep_application.h"

namespace Kratos
{

namespace Python
{

BOOST_PYTHON_MODULE(KratosIsogeometricBRepApplication)
{

    using namespace boost::python;

    class_<KratosIsogeometricBRepApplication,
           KratosIsogeometricBRepApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable > ("KratosIsogeometricBRepApplication")
           ;

    IsogeometricApplication_AddBRepAndLevelSetToPython();

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
