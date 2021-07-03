/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Nov 30, 2017 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_POINT_BASED_CONTROL_GRID_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_POINT_BASED_CONTROL_GRID_TO_PYTHON_H_INCLUDED



// System includes
#include <string>
#include <sstream>

// External includes


// Project includes
#include "includes/define_python.h"
#include "custom_utilities/point_based_control_grid.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

template<class TVariableType, class TFESpaceType>
struct IsogeometricApplication_AddPointBasedControlGrid_Helper
{
    static void Execute();
};

template<class TFESpaceType>
struct IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<double>, TFESpaceType>
{
    typedef Variable<double> VariableType;
    typedef ControlGrid<typename VariableType::Type> ControlGridType;
    typedef PointBasedControlGrid<VariableType, TFESpaceType> PointBasedControlGridType;

    static void Execute(pybind11::module& m)
    {
        std::stringstream ss;
        ss << TFESpaceType::StaticType() << "PointBasedDoubleControlGrid";
        class_<PointBasedControlGridType, typename PointBasedControlGridType::Pointer, ControlGridType>
        (m, ss.str().c_str())
        .def(init<const VariableType&, typename TFESpaceType::Pointer>())
        .def("__str__", &PrintObject<PointBasedControlGridType>)
        ;
    }
};

template<class TFESpaceType>
struct IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<array_1d<double, 3> >, TFESpaceType>
{
    typedef Variable<array_1d<double, 3> > VariableType;
    typedef ControlGrid<typename VariableType::Type> ControlGridType;
    typedef PointBasedControlGrid<VariableType, TFESpaceType> PointBasedControlGridType;

    static void Execute(pybind11::module& m)
    {
        std::stringstream ss;
        ss << TFESpaceType::StaticType() << "PointBasedArray1DControlGrid";
        class_<PointBasedControlGridType, typename PointBasedControlGridType::Pointer, ControlGridType>
        (m, ss.str().c_str())
        .def(init<const VariableType&, typename TFESpaceType::Pointer>())
        .def("__str__", &PrintObject<PointBasedControlGridType>)
        ;
    }
};

template<class TFESpaceType>
struct IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<Vector>, TFESpaceType>
{
    typedef Variable<Vector> VariableType;
    typedef ControlGrid<typename VariableType::Type> ControlGridType;
    typedef PointBasedControlGrid<VariableType, TFESpaceType> PointBasedControlGridType;

    static void Execute(pybind11::module& m)
    {
        std::stringstream ss;
        ss << TFESpaceType::StaticType() << "PointBasedVectorControlGrid";
        class_<PointBasedControlGridType, typename PointBasedControlGridType::Pointer, ControlGridType>
        (m, ss.str().c_str())
        .def(init<const VariableType&, typename TFESpaceType::Pointer>())
        .def("__str__", &PrintObject<PointBasedControlGridType>)
        ;
    }
};

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_POINT_BASED_CONTROL_GRID_TO_PYTHON_H_INCLUDED  defined

