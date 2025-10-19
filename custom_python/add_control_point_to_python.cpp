/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        KratosIsogeometricApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 11, 2017 $
//
//

// System includes
#include <string>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "custom_utilities/control_point.h"
#include "custom_python/add_control_point_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

inline double ControlPoint_GetX(ControlPoint<double>& rDummy)
{
    return rDummy.X();
}

inline double ControlPoint_GetWX(ControlPoint<double>& rDummy)
{
    return rDummy.WX();
}

inline void ControlPoint_SetWX(ControlPoint<double>& rDummy, double newWX)
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

inline void ControlPoint_SetWY(ControlPoint<double>& rDummy, double newWY)
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

inline void ControlPoint_SetWZ(ControlPoint<double>& rDummy, double newWZ)
{
    rDummy.WZ() = newWZ;
}

inline double ControlPoint_GetW(ControlPoint<double>& rDummy)
{
    return rDummy.W();
}

inline void ControlPoint_SetW(ControlPoint<double>& rDummy, double newW)
{
    rDummy.W() = newW;
}

inline void ControlPoint_ApplyTransformation(ControlPoint<double>& rDummy, const Transformation<double>& trans)
{
    rDummy.ApplyTransformation(trans);
}

////////////////////////////////////////

template<typename TDataType, typename TWeightType>
void IsogeometricApplication_AddControlPointToPythonImpl(const std::string& Prefix)
{
    typedef ControlPoint<TDataType, TWeightType> ControlPointType;

    class_<ControlPointType, typename ControlPointType::Pointer>
    ((Prefix+"ControlPoint").c_str(), init<>())
    .def(init<TDataType, TDataType, TDataType, TWeightType>())
    .add_property("WX", ControlPoint_GetWX, ControlPoint_SetWX)
    .add_property("WY", ControlPoint_GetWY, ControlPoint_SetWY)
    .add_property("WZ", ControlPoint_GetWZ, ControlPoint_SetWZ)
    .add_property("W", ControlPoint_GetW, ControlPoint_SetW)
    .def("ApplyTransformation", &ControlPoint_ApplyTransformation)
    .def("X", &ControlPoint_GetX)
    .def("Y", &ControlPoint_GetY)
    .def("Z", &ControlPoint_GetZ)
    .def(self_ns::str(self))
    ;

    class_<Variable<ControlPointType>, bases<VariableData>, boost::noncopyable >( (Prefix+"ControlPointVariable").c_str(), no_init )
    .def( self_ns::str( self ) )
    ;
}

void IsogeometricApplication_AddControlPointToPython()
{

    IsogeometricApplication_AddControlPointToPythonImpl<KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE>("");
    IsogeometricApplication_AddControlPointToPythonImpl<KRATOS_COMPLEX_TYPE, KRATOS_DOUBLE_TYPE>("Complex");

}

}  // namespace Python.

} // Namespace Kratos
