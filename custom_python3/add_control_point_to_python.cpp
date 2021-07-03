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

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "custom_utilities/control_point.h"
#include "custom_python3/add_control_point_to_python.h"



namespace Kratos
{

namespace Python
{

using namespace pybind11;

////////////////////////////////////////

inline double ControlPoint_GetX(ControlPoint<double>& rDummy)
{
    return rDummy.X();
}

inline double ControlPoint_GetWX(ControlPoint<double>& rDummy)
{
    return rDummy.WX();
}

inline void ControlPoint_SetWX(ControlPoint<double>& rDummy, const double& newWX)
{
    rDummy.WX() = newWX;
}

inline double ControlPoint_GetY(ControlPoint<double>& rDummy)
{
    return rDummy.Y();
}

inline double ControlPoint_GetWY(ControlPoint<double>& rDummy)
{
    return rDummy.WY();
}

inline void ControlPoint_SetWY(ControlPoint<double>& rDummy, const double& newWY)
{
    rDummy.WY() = newWY;
}

inline double ControlPoint_GetWZ(ControlPoint<double>& rDummy)
{
    return rDummy.WZ();
}

inline double ControlPoint_GetZ(ControlPoint<double>& rDummy)
{
    return rDummy.Z();
}

inline void ControlPoint_SetWZ(ControlPoint<double>& rDummy, const double& newWZ)
{
    rDummy.WZ() = newWZ;
}

inline double ControlPoint_GetW(ControlPoint<double>& rDummy)
{
    return rDummy.W();
}

inline void ControlPoint_SetW(ControlPoint<double>& rDummy, const double& newW)
{
    rDummy.W() = newW;
}

inline void ControlPoint_ApplyTransformation(ControlPoint<double>& rDummy, const Transformation<double>& trans)
{
    rDummy.ApplyTransformation(trans);
}

////////////////////////////////////////

void IsogeometricApplication_AddControlPointToPython(pybind11::module& m)
{
    class_<ControlPoint<double>, ControlPoint<double>::Pointer>
    (m, "ControlPoint")
    .def(init<>())
    .def(init<const double&, const double&, const double&, const double&>())
    .def_property("WX", ControlPoint_GetWX, ControlPoint_SetWX)
    .def_property("WY", ControlPoint_GetWY, ControlPoint_SetWY)
    .def_property("WZ", ControlPoint_GetWZ, ControlPoint_SetWZ)
    .def_property("W", ControlPoint_GetW, ControlPoint_SetW)
    .def("ApplyTransformation", &ControlPoint_ApplyTransformation)
    .def("X", &ControlPoint_GetX)
    .def("Y", &ControlPoint_GetY)
    .def("Z", &ControlPoint_GetZ)
    .def("__str__", &PrintObject<ControlPoint<double> >)
    ;

    class_<Variable<ControlPoint<double> >, VariableData>(m, "ControlPointVariable")
    .def("__str__", &PrintObject<Variable<ControlPoint<double> > >)
    ;
}

}  // namespace Python.

} // Namespace Kratos

