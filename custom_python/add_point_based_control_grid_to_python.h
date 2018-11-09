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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "custom_utilities/point_based_control_grid.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

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

    static void Execute()
    {
        std::stringstream ss;
        ss << TFESpaceType::StaticType() << "PointBasedDoubleControlGrid";
        class_<PointBasedControlGridType, typename PointBasedControlGridType::Pointer, bases<ControlGridType>, boost::noncopyable>
        (ss.str().c_str(), init<const VariableType&, typename TFESpaceType::Pointer>())
        .def(self_ns::str(self))
        ;
    }
};

template<class TFESpaceType>
struct IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<array_1d<double, 3> >, TFESpaceType>
{
    typedef Variable<array_1d<double, 3> > VariableType;
    typedef ControlGrid<typename VariableType::Type> ControlGridType;
    typedef PointBasedControlGrid<VariableType, TFESpaceType> PointBasedControlGridType;

    static void Execute()
    {
        std::stringstream ss;
        ss << TFESpaceType::StaticType() << "PointBasedArray1DControlGrid";
        class_<PointBasedControlGridType, typename PointBasedControlGridType::Pointer, bases<ControlGridType>, boost::noncopyable>
        (ss.str().c_str(), init<const VariableType&, typename TFESpaceType::Pointer>())
        .def(self_ns::str(self))
        ;
    }
};

template<class TFESpaceType>
struct IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<Vector>, TFESpaceType>
{
    typedef Variable<Vector> VariableType;
    typedef ControlGrid<typename VariableType::Type> ControlGridType;
    typedef PointBasedControlGrid<VariableType, TFESpaceType> PointBasedControlGridType;

    static void Execute()
    {
        std::stringstream ss;
        ss << TFESpaceType::StaticType() << "PointBasedVectorControlGrid";
        class_<PointBasedControlGridType, typename PointBasedControlGridType::Pointer, bases<ControlGridType>, boost::noncopyable>
        (ss.str().c_str(), init<const VariableType&, typename TFESpaceType::Pointer>())
        .def(self_ns::str(self))
        ;
    }
};

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_POINT_BASED_CONTROL_GRID_TO_PYTHON_H_INCLUDED  defined

