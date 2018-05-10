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


// System includes
#include <string>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_utilities/trans/transformation.h"
#include "custom_utilities/trans/translation.h"
#include "custom_utilities/trans/rotation.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

//////////////////////////////////////////////////

void IsogeometricApplication_AddTransformationToPython()
{
    class_<Transformation<double>, Transformation<double>::Pointer, boost::noncopyable>
    ("Transformation", init<>())
    .def("AppendTransformation", &Transformation<double>::AppendTransformation)
    // .def(boost::python::operators<boost::python::op_mul>());
    .def(self_ns::str(self))
    ;

    class_<Translation<double>, Translation<double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("Translation", init<const double&, const double&, const double&>())
    .def(self_ns::str(self))
    ;

    class_<Rotation<0, double>, Rotation<0, double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("RotationX", init<const double&>())
    .def(self_ns::str(self))
    ;

    class_<Rotation<1, double>, Rotation<1, double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("RotationY", init<const double&>())
    .def(self_ns::str(self))
    ;

    class_<Rotation<2, double>, Rotation<2, double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("RotationZ", init<const double&>())
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

