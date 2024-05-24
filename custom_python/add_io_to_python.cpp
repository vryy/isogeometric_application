/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 22 Oct 2015 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part_io.h"
#include "custom_io/bezier_model_part_io.h"
#include "custom_io/isogeometric_model_part_io.h"
#include "custom_python/add_io_to_python.h"

namespace Kratos
{

namespace Python
{

void  IsogeometricApplication_AddIOToPython()
{
    using namespace boost::python;

    class_<IsogeometricModelPartIO, IsogeometricModelPartIO::Pointer, bases<IO>,  boost::noncopyable>(
        "IsogeometricModelPartIO", init<std::string const&>())
    ;

    class_<BezierModelPartIO, BezierModelPartIO::Pointer, bases<ModelPartIO>,  boost::noncopyable>(
        "BezierModelPartIO", init<std::string const&>())
//        .def(init<std::string const&, const Flags>())
    ;
}

}  // namespace Python.

} // Namespace Kratos
