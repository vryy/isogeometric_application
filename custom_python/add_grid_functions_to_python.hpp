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

template<class TGridFrunctionType>
typename TGridFrunctionType::FESpaceType::Pointer GridFunction_GetFESpace(TGridFrunctionType& rDummy)
{
    return rDummy.pFESpace();
}

template<class TGridFrunctionType>
void GridFunction_SetFESpace(TGridFrunctionType& rDummy, typename TGridFrunctionType::FESpaceType::Pointer pNewFESpace)
{
    rDummy.SetFESpace(pNewFESpace);
}

template<class TGridFrunctionType>
typename TGridFrunctionType::ControlGridType::Pointer GridFunction_GetControlGrid(TGridFrunctionType& rDummy)
{
    return rDummy.pControlGrid();
}

template<class TGridFrunctionType>
void GridFunction_SetControlGrid(TGridFrunctionType& rDummy, typename TGridFrunctionType::ControlGridType::Pointer pNewControlGrid)
{
    rDummy.SetControlGrid(pNewControlGrid);
}

template<class TGridFrunctionType>
typename TGridFrunctionType::DataType GridFunction_GetValue(TGridFrunctionType& rDummy, const boost::python::list& xi)
{
    std::vector<double> xi_vec;
    typedef boost::python::stl_input_iterator<double> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& v, std::make_pair(iterator_value_type(xi), iterator_value_type() ) )
    {
        xi_vec.push_back(v);
    }

    return rDummy.GetValue(xi_vec);
}

template<class TGridFrunctionType>
boost::python::list GridFunction_GetDerivative(TGridFrunctionType& rDummy, const boost::python::list& xi)
{
    std::vector<double> xi_vec;
    typedef boost::python::stl_input_iterator<double> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& v, std::make_pair(iterator_value_type(xi), iterator_value_type() ) )
    {
        xi_vec.push_back(v);
    }

    std::vector<typename TGridFrunctionType::DataType> derivatives;
    rDummy.GetDerivative(derivatives, xi_vec);

    boost::python::list results;
    for (std::size_t i = 0; i < derivatives.size(); ++i)
        results.append(derivatives[i]);

    return results;
}

///////////////////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddGridFunctionsToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "ControlPointGridFunction" << TDim << "D";
    typedef GridFunction<TDim, ControlPoint<double> > ControlPointGridFunctionType;
    class_<ControlPointGridFunctionType, typename ControlPointGridFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<ControlPoint<double> >::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<ControlPointGridFunctionType>, GridFunction_SetFESpace<ControlPointGridFunctionType>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<ControlPointGridFunctionType>, GridFunction_SetControlGrid<ControlPointGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue<ControlPointGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative<ControlPointGridFunctionType>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "DoubleGridFunction" << TDim << "D";
    typedef GridFunction<TDim, double> DoubleGridFunctionType;
    class_<DoubleGridFunctionType, typename DoubleGridFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<double>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<DoubleGridFunctionType>, GridFunction_SetFESpace<DoubleGridFunctionType>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<DoubleGridFunctionType>, GridFunction_SetControlGrid<DoubleGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue<DoubleGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative<DoubleGridFunctionType>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "Array1DGridFunction" << TDim << "D";
    typedef GridFunction<TDim, array_1d<double, 3> > Array1DGridFunctionType;
    class_<Array1DGridFunctionType, typename Array1DGridFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<array_1d<double, 3> >::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<Array1DGridFunctionType>, GridFunction_SetFESpace<Array1DGridFunctionType>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<Array1DGridFunctionType>, GridFunction_SetControlGrid<Array1DGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue<Array1DGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative<Array1DGridFunctionType>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "VectorGridFunction" << TDim << "D";
    typedef GridFunction<TDim, Vector> VectorGridFunctionType;
    class_<VectorGridFunctionType, typename VectorGridFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<Vector>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<VectorGridFunctionType>, GridFunction_SetFESpace<VectorGridFunctionType>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<VectorGridFunctionType>, GridFunction_SetControlGrid<VectorGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue<VectorGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative<VectorGridFunctionType>)
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

