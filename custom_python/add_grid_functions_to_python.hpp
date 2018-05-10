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
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/grid_function.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

template<int TDim, typename TDataType>
typename FESpace<TDim>::Pointer GridFunction_GetFESpace(GridFunction<TDim, TDataType>& rDummy)
{
    return rDummy.pFESpace();
}

template<int TDim, typename TDataType>
void GridFunction_SetFESpace(GridFunction<TDim, TDataType>& rDummy, typename FESpace<TDim>::Pointer pNewFESpace)
{
    rDummy.SetFESpace(pNewFESpace);
}

template<int TDim, typename TDataType>
typename ControlGrid<TDataType>::Pointer GridFunction_GetControlGrid(GridFunction<TDim, TDataType>& rDummy)
{
    return rDummy.pControlGrid();
}

template<int TDim, typename TDataType>
void GridFunction_SetControlGrid(GridFunction<TDim, TDataType>& rDummy, typename ControlGrid<TDataType>::Pointer pNewControlGrid)
{
    rDummy.SetControlGrid(pNewControlGrid);
}

template<int TDim, typename TDataType>
TDataType GridFunction_GetValue(GridFunction<TDim, TDataType>& rDummy, const boost::python::list& xi)
{
    std::vector<double> xi_vec;
    typedef boost::python::stl_input_iterator<double> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& v, std::make_pair(iterator_value_type(xi), iterator_value_type() ) )
    {
        xi_vec.push_back(v);
    }

    return rDummy.GetValue(xi_vec);
}

///////////////////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddGridFunctionsToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "ControlPointGridFunction" << TDim << "D";
    class_<GridFunction<TDim, ControlPoint<double> >, typename GridFunction<TDim, ControlPoint<double> >::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<ControlPoint<double> >::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, ControlPoint<double> >, GridFunction_SetFESpace<TDim, ControlPoint<double> >)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, ControlPoint<double> >, GridFunction_SetControlGrid<TDim, ControlPoint<double> >)
    .def("GetValue", &GridFunction_GetValue<TDim, ControlPoint<double> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "DoubleGridFunction" << TDim << "D";
    class_<GridFunction<TDim, double>, typename GridFunction<TDim, double>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<double>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, double>, GridFunction_SetFESpace<TDim, double>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, double>, GridFunction_SetControlGrid<TDim, double>)
    .def("GetValue", &GridFunction_GetValue<TDim, double>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "Array1DGridFunction" << TDim << "D";
    class_<GridFunction<TDim, array_1d<double, 3> >, typename GridFunction<TDim, array_1d<double, 3> >::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<array_1d<double, 3> >::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, array_1d<double, 3> >, GridFunction_SetFESpace<TDim, array_1d<double, 3> >)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, array_1d<double, 3> >, GridFunction_SetControlGrid<TDim, array_1d<double, 3> >)
    .def("GetValue", &GridFunction_GetValue<TDim, array_1d<double, 3> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "VectorGridFunction" << TDim << "D";
    class_<GridFunction<TDim, Vector>, typename GridFunction<TDim, Vector>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<Vector>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, Vector>, GridFunction_SetFESpace<TDim, Vector>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, Vector>, GridFunction_SetControlGrid<TDim, Vector>)
    .def("GetValue", &GridFunction_GetValue<TDim, Vector>)
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

