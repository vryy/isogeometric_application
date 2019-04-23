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
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "includes/model_part_io.h"
#include "custom_io/bezier_model_part_io.h"
#include "custom_io/isogeometric_model_part_io.h"
#include "add_io_to_python.h"

namespace Kratos
{

namespace Python
{

void  IsogeometricApplication_AddIOToPython(pybind11::module& m)
{
    pybind11::class_<IsogeometricModelPartIO, IsogeometricModelPartIO::Pointer, IO>(m, "IsogeometricModelPartIO")
    .def(pybind11::init<std::string const&>())
    ;

    pybind11::class_<BezierModelPartIO, BezierModelPartIO::Pointer, ModelPartIO>(m, "BezierModelPartIO")
    .def(pybind11::init<std::string const&>())
//   .def(pybind11::init<std::string const&, const Flags>())
    ;
}

}  // namespace Python.

} // Namespace Kratos

