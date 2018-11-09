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
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_python/add_control_grids_to_python.h"



namespace Kratos
{

namespace Python
{

using namespace boost::python;

void IsogeometricApplication_AddControlPoint()
{
    class_<ControlPoint<double>, ControlPoint<double>::Pointer>
    ("ControlPoint", init<>())
    .def(init<const double&, const double&, const double&, const double&>())
    .add_property("WX", ControlPoint_GetWX, ControlPoint_SetWX)
    .add_property("WY", ControlPoint_GetWY, ControlPoint_SetWY)
    .add_property("WZ", ControlPoint_GetWZ, ControlPoint_SetWZ)
    .add_property("W", ControlPoint_GetW, ControlPoint_SetW)
    .def("ApplyTransformation", &ControlPoint_ApplyTransformation)
    .def(self_ns::str(self))
    ;

    class_<Variable<ControlPoint<double> >, bases<VariableData>, boost::noncopyable >( "ControlPointVariable", no_init )
    .def( self_ns::str( self ) )
    ;
}

void IsogeometricApplication_AddControlGrids()
{
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<ControlGrid<ControlPoint<double> >, ControlGrid<ControlPoint<double> >::Pointer, boost::noncopyable>
    ("ControlPointControlGrid", init<>())
    .def("Size", &ControlGrid<ControlPoint<double> >::Size)
    .def("size", &ControlGrid<ControlPoint<double> >::Size)
    .def("__setitem__", &ControlGrid_SetItem<ControlPoint<double> >)
    .def("__getitem__", &ControlGrid_GetItem<ControlPoint<double> >)
    .def(self_ns::str(self))
    ;

    class_<ControlGrid<double>, ControlGrid<double>::Pointer, boost::noncopyable>
    ("DoubleControlGrid", init<>())
    .def("Size", &ControlGrid<double>::Size)
    .def("size", &ControlGrid<double>::Size)
    .def("__setitem__", &ControlGrid_SetItem<double>)
    .def("__getitem__", &ControlGrid_GetItem<double>)
    .def(self_ns::str(self))
    ;

    class_<ControlGrid<array_1d<double, 3> >, ControlGrid<array_1d<double, 3> >::Pointer, boost::noncopyable>
    ("Array1DControlGrid", init<>())
    .def("Size", &ControlGrid<array_1d<double, 3> >::Size)
    .def("size", &ControlGrid<array_1d<double, 3> >::Size)
    .def("__setitem__", &ControlGrid_SetItem<array_1d<double, 3> >)
    .def("__getitem__", &ControlGrid_GetItem<array_1d<double, 3> >)
    .def(self_ns::str(self))
    ;

    class_<ControlGrid<Vector>, ControlGrid<Vector>::Pointer, boost::noncopyable>
    ("VectorControlGrid", init<>())
    .def("Size", &ControlGrid<Vector>::Size)
    .def("size", &ControlGrid<Vector>::Size)
    .def("__setitem__", &ControlGrid_SetItem<Vector>)
    .def("__getitem__", &ControlGrid_GetItem<Vector>)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<UnstructuredControlGrid<ControlPoint<double> >, UnstructuredControlGrid<ControlPoint<double> >::Pointer, bases<ControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("UnstructuredControlPointGrid", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<UnstructuredControlGrid<double>, UnstructuredControlGrid<double>::Pointer, bases<ControlGrid<double> >, boost::noncopyable>
    ("UnstructuredDoubleControlGrid", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<UnstructuredControlGrid<array_1d<double, 3> >, UnstructuredControlGrid<array_1d<double, 3> >::Pointer, bases<ControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("UnstructuredArray1DControlGrid", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<UnstructuredControlGrid<Vector>, UnstructuredControlGrid<Vector>::Pointer, bases<ControlGrid<Vector> >, boost::noncopyable>
    ("UnstructuredVectorControlGrid", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////
}

void IsogeometricApplication_AddControlGridsToPython()
{

    /////////////////////////////////////////////////////////////////
    ///////////////////CONTROL POINT/////////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddControlPoint();

    /////////////////////////////////////////////////////////////////
    ///////////////////////CONTROL GRIDS/////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddControlGrids();

    /////////////////////////////////////////////////////////////////

    class_<ControlGridLibrary, ControlGridLibrary::Pointer, boost::noncopyable>
    ("ControlGridLibrary", init<>())
    .def("CreateLinearControlPointGrid", &ControlGridLibrary_CreateLinearControlPointGrid)
    .def("CreateRectangularControlPointGrid", &ControlGridLibrary_CreateRectangularControlPointGrid1)
    .def("CreateRectangularControlPointGrid", &ControlGridLibrary_CreateRectangularControlPointGrid2)
    .def("CreateCubicControlPointGrid", &ControlGridLibrary_CreateCubicControlPointGrid1)
    .def("CreateCubicControlPointGrid", &ControlGridLibrary_CreateCubicControlPointGrid2)
    .def("CreateLinearZeroDoubleControlGrid", &ControlGridLibrary_CreateLinearZeroControlGridWithVariable<Variable<double> >)
    .def("CreateRectangularZeroDoubleControlGrid", &ControlGridLibrary_CreateRectangularZeroControlGridWithVariable<Variable<double> >)
    .def("CreateCubicZeroDoubleControlGrid", &ControlGridLibrary_CreateCubicZeroControlGridWithVariable<Variable<double> >)
    .def("CreateLinearZeroArray1DControlGrid", &ControlGridLibrary_CreateLinearZeroControlGridWithVariable<Variable<array_1d<double, 3> > >)
    .def("CreateRectangularZeroArray1DControlGrid", &ControlGridLibrary_CreateRectangularZeroControlGridWithVariable<Variable<array_1d<double, 3> > >)
    .def("CreateCubicZeroArray1DControlGrid", &ControlGridLibrary_CreateCubicZeroControlGridWithVariable<Variable<array_1d<double, 3> > >)
    ;

}

}  // namespace Python.

} // Namespace Kratos

